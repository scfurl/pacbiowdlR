extern crate clap;
extern crate serde_json;
extern crate anyhow;
extern crate serde;
extern crate rayon;
extern crate indicatif;
extern crate num_cpus;

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use clap::Parser;
use serde_json::{Value, json};
use anyhow::{Result, Context};
use rayon::prelude::*;
use std::sync::Mutex;
use indicatif::{ProgressBar, ProgressStyle};

/// A CLI tool to copy files from softlinks listed in a JSON file to a new location
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
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

fn main() -> Result<()> {
    let args = Args::parse();

    // Read the input JSON file
    let json_content = fs::read_to_string(&args.input)
        .with_context(|| format!("Failed to read input file: {:?}", args.input))?;
    
    let json_data: Value = serde_json::from_str(&json_content)
        .context("Failed to parse JSON")?;

    // Ensure output directory exists
    if !args.dry_run {
        fs::create_dir_all(&args.output)
            .with_context(|| format!("Failed to create output directory: {:?}", args.output))?;
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
                        
                        // Construct destination path and make it absolute
                        let file_name = real_path.file_name()
                            .ok_or_else(|| anyhow::anyhow!("Failed to get filename from: {:?}", real_path))?;
                        
                        let destination_path = args.output.join(file_name);
                        
                        // Convert to absolute path for the JSON output
                        let absolute_destination_path = if destination_path.is_relative() {
                            std::env::current_dir()
                                .context("Failed to get current directory")?
                                .join(&destination_path)
                        } else {
                            destination_path.clone()
                        };
                        
                        // Create file operation
                        let operation = FileOperation {
                            key: key.clone(),
                            source: real_path,
                            destination: destination_path,
                        };
                        
                        file_operations.push(operation);
                        
                        // Add to new JSON with updated absolute path
                        new_json.insert(key, json!(absolute_destination_path.to_string_lossy().to_string()));
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
            println!("  Copy {:?} to {:?}", op.source, op.destination);
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
                // Update progress bar message with current file
                progress_bar.set_message(format!("Copying {}", op.key));
                
                if let Err(e) = fs::copy(&op.source, &op.destination) {
                    let error_msg = format!("Failed to copy {:?} to {:?}: {}", op.source, op.destination, e);
                    errors.lock().unwrap().push(error_msg);
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

    // Write new JSON file
    let json_output_path = args.json.clone().unwrap_or_else(|| args.output.join("outputs.json"));
    
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