extern crate rayon;
extern crate flate2;

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use rayon::prelude::*;
use std::collections::BTreeMap;
use flate2::read::GzDecoder;
use std::process::{Command, Stdio};

#[derive(Debug, Clone)]
struct BedEntry {
    chrom: String,
    start: u64,
    end: u64,
    value: f64,
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    
    if args.len() < 3 || args.len() > 4 {
        eprintln!("Usage: {} <input_file> <output_file> [threads]", args[0]);
        eprintln!("  input_file: can be plain text or gzipped (.gz extension)");
        eprintln!("  threads: optional number of threads to use (default: use all available)");
        std::process::exit(1);
    }
    
    let input_file = &args[1];
    let output_file = &args[2];
    
    // Check if input file exists
    if !Path::new(input_file).exists() {
        eprintln!("Error: Input file '{}' does not exist", input_file);
        std::process::exit(1);
    }
    
    // Configure thread pool if specified
    if let Some(threads_str) = args.get(3) {
        if let Ok(num_threads) = threads_str.parse::<usize>() {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .expect("Failed to build thread pool");
            println!("Using {} threads", num_threads);
        } else {
            eprintln!("Invalid thread count: {}", threads_str);
            std::process::exit(1);
        }
    } else {
        println!("Using {} threads (all available)", rayon::current_num_threads());
    }
    
    // Detect if input is gzipped
    println!("Reading input file: {}", input_file);
    let is_gzipped = Path::new(input_file).extension().map_or(false, |ext| ext == "gz");
    if is_gzipped {
        println!("Detected gzipped input file");
    }
    
    // Try the gunzip+awk approach first
    println!("Processing file with native commands...");
    
    let gunzip_cmd = if is_gzipped {
        if cfg!(target_os = "macos") {
            "gzcat"  // macOS uses gzcat
        } else {
            "zcat"   // Linux typically uses zcat
        }
    } else {
        "cat"
    };
    
    let output_cmd = Command::new("bash")
        .arg("-c")
        .arg(format!(
            "{} \"{}\" | awk -F'\\t' '
            # Skip header line
            NR > 1 {{
                # Only take rows where column 6 is -1.0 
                # or any other non -1.0 value we want to keep
                if ($6 == \"-1.0\") {{
                    value = 0.0;
                }} else {{
                    value = $6;
                }}
                # Output the 4 columns we want
                print $1\"\\t\"$2\"\\t\"$3\"\\t\"value;
            }}' > \"{}\"", 
            gunzip_cmd, input_file, output_file
        ))
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()?;
    
    let output = output_cmd.wait_with_output()?;
    
    if output.status.success() {
        println!("Successfully processed file using native commands");
        println!("Output written to: {}", output_file);
        return Ok(());
    } else {
        // Print the error if the command failed
        if !output.stderr.is_empty() {
            eprintln!("Error from command: {}", String::from_utf8_lossy(&output.stderr));
        }
        
        println!("Native command failed, trying Rust implementation...");
    }
    
    // Fallback to Rust implementation if native commands failed
    let mut entries = Vec::new();
    let mut line_count = 0;
    let mut is_header = true;  // Assume first line is header
    
    // Create appropriate reader based on file type
    if is_gzipped {
        let file = File::open(input_file)?;
        let gz = GzDecoder::new(file);
        let reader = BufReader::with_capacity(1024 * 1024, gz); // 1MB buffer
        
        for line in reader.lines() {
            let line = match line {
                Ok(line) => line,
                Err(e) => {
                    eprintln!("Error reading line: {}", e);
                    continue;
                }
            };
            
            line_count += 1;
            
            if is_header {
                is_header = false;
                println!("Skipping header line: {}", line);
                continue;
            }
            
            if line_count % 1_000_000 == 0 {
                println!("Read {} million lines...", line_count / 1_000_000);
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 6 {
                let chrom = fields[0].to_string();
                let start = match fields[1].parse::<u64>() {
                    Ok(val) => val,
                    Err(_) => {
                        eprintln!("Invalid start position in line {}: {}", line_count, fields[1]);
                        continue;
                    }
                };
                let end = match fields[2].parse::<u64>() {
                    Ok(val) => val,
                    Err(_) => {
                        eprintln!("Invalid end position in line {}: {}", line_count, fields[2]);
                        continue;
                    }
                };
                let value = match fields[5].parse::<f64>() {
                    Ok(val) => val,
                    Err(_) => {
                        eprintln!("Invalid value in line {}: {}", line_count, fields[5]);
                        continue;
                    }
                };
                
                entries.push(BedEntry {
                    chrom,
                    start,
                    end,
                    value,
                });
            }
        }
    } else {
        let file = File::open(input_file)?;
        let reader = BufReader::new(file);
        
        for line in reader.lines() {
            let line = line?;
            line_count += 1;
            
            if is_header {
                is_header = false;
                println!("Skipping header line: {}", line);
                continue;
            }
            
            if line_count % 1_000_000 == 0 {
                println!("Read {} million lines...", line_count / 1_000_000);
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 6 {
                entries.push(BedEntry {
                    chrom: fields[0].to_string(),
                    start: fields[1].parse::<u64>().unwrap_or(0),
                    end: fields[2].parse::<u64>().unwrap_or(0),
                    value: fields[5].parse::<f64>().unwrap_or(0.0),
                });
            }
        }
    }
    
    println!("Read {} lines total (excluding header)", line_count - 1);
    println!("Parsed {} valid BED entries", entries.len());
    
    if entries.is_empty() {
        eprintln!("Error: No valid entries parsed from the file");
        std::process::exit(1);
    }
    
    // Group entries by chromosome for parallel processing
    let mut chrom_groups: BTreeMap<String, Vec<BedEntry>> = BTreeMap::new();
    for entry in &entries {
        chrom_groups.entry(entry.chrom.clone()).or_default().push(entry.clone());
    }
    
    // Process each chromosome in parallel
    let merged_groups: Vec<(String, Vec<BedEntry>)> = chrom_groups
        .into_par_iter()
        .map(|(chrom, mut entries)| {
            // Sort entries by start position
            entries.sort_by_key(|e| e.start);
            
            let mut merged = Vec::new();
            if entries.is_empty() {
                return (chrom, merged);
            }
            
            let mut current = entries[0].clone();
            
            for next in entries.iter().skip(1) {
                if current.end == next.start 
                   && current.value == -1.0 
                   && next.value == -1.0 {
                    current.end = next.end;
                } else {
                    merged.push(current);
                    current = next.clone();
                }
            }
            merged.push(current);
            
            (chrom, merged)
        })
        .collect();
    
    // Write output sequentially (file writes can't be easily parallelized)
    let mut output = File::create(output_file)?;
    let mut total_merged = 0;
    
    // Sort chromosomes for consistent output
    let mut sorted_groups = merged_groups;
    sorted_groups.sort_by(|a, b| a.0.cmp(&b.0));
    
    for (_, entries) in sorted_groups {
        for entry in entries {
            let output_value = if entry.value == -1.0 { 0.0 } else { entry.value };
            writeln!(output, "{}\t{}\t{}\t{}", 
                     entry.chrom, entry.start, entry.end, output_value)?;
            total_merged += 1;
        }
    }
    
    println!("Merged {} entries into {} entries", entries.len(), total_merged);
    println!("Output written to: {}", output_file);
    
    Ok(())
}