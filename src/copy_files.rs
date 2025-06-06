extern crate clap;
extern crate serde_json;
extern crate anyhow;
extern crate serde;
extern crate rayon;
extern crate indicatif;
extern crate num_cpus;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use clap::Parser;
use serde_json::{Value, json};
use anyhow::{Result, Context};
use rayon::prelude::*;
use std::sync::Mutex;
use indicatif::{ProgressBar, ProgressStyle};

/// A CLI tool to copy files from softlinks listed in a JSON file to a new location
#[derive(Parser, Debug)]
#[command(name = "copy_pbWDL")]
#[command(version = "0.1")]
#[command(about = "Copy files from softlinks listed in a JSON file to a new location", long_about = "This tool reads a JSON file containing softlinks and copies the files to a specified output directory. It can also move the haplotagged.bam file if specified.")]
struct Args {
    /// Input JSON file path
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory where files will be copied
    #[arg(short, long)]
    output: PathBuf,

    /// Output JSON file path (optional, defaults to output_dir/outputs.json)
    #[arg(short = 'j', long)]
    json: Option<PathBuf>,

    /// Dry run - only show what would be copied
    #[arg(short, long)]
    dry_run: bool,

    /// Number of threads to use for copying (defaults to number of CPU cores)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Move the haplotagged.bam file instead of copying it
    #[arg(short = 'm', long)]
    move_haplotagged_bam: bool,
}

#[derive(Debug)]
struct FileOperation {
    key: String,
    source: PathBuf,
    destination: PathBuf,
}

fn is_likely_file_path(value: &str) -> bool {
    // Check if the string looks like a file path
    // Returns true if:
    // 1. Starts with "/" (absolute path)
    // 2. Contains a "/" (has directory structure)
    // 3. Has a file extension (contains "." in the last component)
    if value.starts_with('/') {
        return true;
    }
    
    if value.contains('/') {
        return true;
    }
    
    // Check for file extension in the last component
    if let Some(last_component) = value.split('/').last() {
        if last_component.contains('.') {
            return true;
        }
    }
    
    false
}

fn normalize_path(path: &Path) -> Result<PathBuf> {
    let mut result = PathBuf::new();
    
    for component in path.components() {
        match component {
            std::path::Component::ParentDir => {
                result.pop();
            }
            std::path::Component::CurDir => {
                // do nothing
            }
            _ => {
                result.push(component);
            }
        }
    }
    
    Ok(result)
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Get absolute path for output directory
    let absolute_output_dir = if args.output.is_relative() {
        let current_dir = std::env::current_dir()
            .context("Failed to get current directory")?;
        let full_path = current_dir.join(&args.output);
        normalize_path(&full_path)?
    } else {
        normalize_path(&args.output)?
    };

    // Read the input JSON file
    let json_content = fs::read_to_string(&args.input)
        .with_context(|| format!("Failed to read input file: {:?}", args.input))?;
    
    let json_data: Value = serde_json::from_str(&json_content)
        .context("Failed to parse JSON")?;

    // Ensure output directory exists
    if !args.dry_run {
        fs::create_dir_all(&absolute_output_dir)
            .with_context(|| format!("Failed to create output directory: {:?}", absolute_output_dir))?;
    }

    // Process each file in the JSON
    let mut file_operations = Vec::new();
    let mut new_json = HashMap::new();

    if let Value::Object(map) = json_data {
        for (key, value) in map {
            match value {
                Value::String(ref string_value) => {
                    // Check if this string looks like a file path
                    if is_likely_file_path(&string_value) {
                        let source_path = PathBuf::from(&string_value);
                        
                        // Check if the path exists before trying to canonicalize
                        if !source_path.exists() {
                            eprintln!("Warning: Path does not exist for {}: {:?}", key, source_path);
                            new_json.insert(key.clone(), value.clone());
                            continue;
                        }
                        
                        // Get the real path after following softlinks
                        let real_path = match fs::canonicalize(&source_path) {
                            Ok(path) => path,
                            Err(e) => {
                                eprintln!("Warning: Failed to resolve path for {}: {:?} - {}", key, source_path, e);
                                new_json.insert(key.clone(), value.clone());
                                continue;
                            }
                        };
                        
                        // Construct destination path
                        let file_name = real_path.file_name()
                            .ok_or_else(|| anyhow::anyhow!("Failed to get filename from: {:?}", real_path))?;
                        
                        let destination_path = absolute_output_dir.join(file_name);
                        
                        // Create file operation with the original destination path for display
                        let display_destination = args.output.join(file_name);
                        let operation = FileOperation {
                            key: key.clone(),
                            source: real_path,
                            destination: display_destination,
                        };
                        
                        file_operations.push(operation);
                        
                        // Add to new JSON with the absolute path
                        new_json.insert(key, json!(destination_path.to_string_lossy().to_string()));
                    } else {
                        // Not a file path - preserve as-is
                        if args.dry_run {
                            println!("Preserving non-path value for {}: {}", key, string_value);
                        }
                        new_json.insert(key.clone(), value.clone());
                    }
                },
                _ => {
                    // Preserve non-string values as-is
                    new_json.insert(key.clone(), value.clone());
                }
            }
        }
    } else {
        anyhow::bail!("JSON root must be an object");
    }

    // Execute file operations
    if args.dry_run {
        println!("Dry run - the following operations would be performed:");
        for op in &file_operations {
            let operation_type = if args.move_haplotagged_bam && op.key.contains("haplotagged_bam") {
                "Move"
            } else {
                "Copy"
            };
            println!("  {} {:?} to {:?}", operation_type, op.source, op.destination);
        }
    } else {
        println!("Copying files...");
        
        // Set up the thread pool
        let num_threads = args.threads.unwrap_or_else(|| num_cpus::get());
        println!("Using {} threads for parallel copying", num_threads);
        
        // Create a thread pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .context("Failed to create thread pool")?;
        
        // Create progress bar
        let progress_bar = ProgressBar::new(file_operations.len() as u64);
        progress_bar.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} files ({percent}%) {msg}")?
                .progress_chars("=>-")
        );
        
        // Parallel copy operations
        let errors: Mutex<Vec<String>> = Mutex::new(Vec::new());
        pool.install(|| {
            file_operations.par_iter().for_each(|op| {
                // Determine if we should move or copy this file
                let is_haplotagged_bam = op.key.contains("haplotagged_bam");
                let should_move = args.move_haplotagged_bam && is_haplotagged_bam;
                
                // Update progress bar message with current file and operation
                let operation_msg = if should_move { "Moving" } else { "Copying" };
                progress_bar.set_message(format!("{} {}", operation_msg, op.key));
                
                // Use the absolute path for the actual copying/moving
                let absolute_destination = absolute_output_dir.join(op.source.file_name().unwrap());
                
                if should_move {
                    // Move the file
                    if let Err(e) = fs::rename(&op.source, &absolute_destination) {
                        let error_msg = format!("Failed to move {:?} to {:?}: {}", op.source, absolute_destination, e);
                        errors.lock().unwrap().push(error_msg);
                    }
                } else {
                    // Copy the file
                    if let Err(e) = fs::copy(&op.source, &absolute_destination) {
                        let error_msg = format!("Failed to copy {:?} to {:?}: {}", op.source, absolute_destination, e);
                        errors.lock().unwrap().push(error_msg);
                    }
                }
                
                progress_bar.inc(1);
            });
        });
        
        progress_bar.finish_with_message("Done!");
        
        // Check for any errors
        let errors = errors.into_inner().unwrap();
        if !errors.is_empty() {
            for error in &errors {
                eprintln!("Error: {}", error);
            }
            anyhow::bail!("{} files failed to copy", errors.len());
        }
    }

    // Determine the output JSON path and make it absolute
    let json_output_path = match args.json {
        Some(path) if path.is_relative() => {
            let current_dir = std::env::current_dir()
                .context("Failed to get current directory")?;
            normalize_path(&current_dir.join(path))?
        },
        Some(path) => normalize_path(&path)?,
        None => absolute_output_dir.join("outputs.json"),
    };
    
    if args.dry_run {
        println!("\nWould write new JSON to: {:?}", json_output_path);
        println!("JSON content (preview):");
        println!("{}", serde_json::to_string_pretty(&new_json)?);
    } else {
        let json_output = serde_json::to_string_pretty(&new_json)?;
        fs::write(&json_output_path, json_output)
            .with_context(|| format!("Failed to write output JSON to: {:?}", json_output_path))?;
        println!("\nWrote new JSON to: {:?}", json_output_path);
    }

    println!("\nOperation completed successfully!");
    if !file_operations.is_empty() {
        println!("Copied {} files", file_operations.len());
    } else {
        println!("No files were copied (no valid file paths found)");
    }
    
    Ok(())
}