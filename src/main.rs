extern crate clap;
extern crate serde_json;
extern crate anyhow;


use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use clap::Parser;
use serde_json::Value;
use anyhow::{Result, Context};

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
}

#[derive(Debug)]
struct FileOperation {
    key: String,
    source: PathBuf,
    destination: PathBuf,
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
            if let Value::String(file_path) = value {
                let source_path = PathBuf::from(file_path);
                
                // Get the real path after following softlinks
                let real_path = fs::canonicalize(&source_path)
                    .with_context(|| format!("Failed to resolve path: {:?}", source_path))?;
                
                // Construct destination path
                let file_name = real_path.file_name()
                    .ok_or_else(|| anyhow::anyhow!("Failed to get filename from: {:?}", real_path))?;
                
                let destination_path = args.output.join(file_name);
                
                // Create file operation
                let operation = FileOperation {
                    key: key.clone(),
                    source: real_path,
                    destination: destination_path.clone(),
                };
                
                file_operations.push(operation);
                
                // Add to new JSON with updated path
                new_json.insert(key, destination_path.to_string_lossy().to_string());
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
        for op in &file_operations {
            println!("  Copying {} from {:?} to {:?}", op.key, op.source, op.destination);
            fs::copy(&op.source, &op.destination)
                .with_context(|| format!("Failed to copy {:?} to {:?}", op.source, op.destination))?;
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
    Ok(())
}