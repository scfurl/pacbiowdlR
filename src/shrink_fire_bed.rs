extern crate rayon;
extern crate flate2;

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use rayon::prelude::*;
use std::collections::BTreeMap;
use flate2::read::GzDecoder;

#[derive(Debug, Clone)]
struct BedEntry {
    chrom: String,
    start: u64,
    end: u64,
    value: f64,
}

impl BedEntry {
    fn merge_with(&mut self, other: &BedEntry) {
        self.end = other.end;
    }
}

/// Reads lines from a file, automatically detecting gzip compression based on file extension
fn read_lines<P: AsRef<Path>>(path: P) -> io::Result<Vec<String>> {
    let file = File::open(&path)?;
    let file_path = path.as_ref();
    
    if let Some(extension) = file_path.extension() {
        if extension == "gz" {
            // Handle gzipped file
            let decoder = GzDecoder::new(file);
            let reader = BufReader::new(decoder);
            return reader.lines().collect();
        }
    }
    
    // Handle regular file
    let reader = BufReader::new(file);
    reader.lines().collect()
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
    
    // Detect if input is gzipped and read lines accordingly
    println!("Reading input file: {}", input_file);
    let is_gzipped = Path::new(input_file).extension().map_or(false, |ext| ext == "gz");
    if is_gzipped {
        println!("Detected gzipped input file");
    }
    
    // Read the input file lines
    let lines = read_lines(input_file)?;
    println!("Read {} lines from input file", lines.len());
    
    // Parse entries in parallel
    let entries: Vec<BedEntry> = lines
        .par_iter()
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 6 {
                Some(BedEntry {
                    chrom: fields[0].to_string(),
                    start: fields[1].parse::<u64>().unwrap_or(0),
                    end: fields[2].parse::<u64>().unwrap_or(0),
                    value: fields[5].parse::<f64>().unwrap_or(0.0),
                })
            } else {
                None
            }
        })
        .collect();
    
    println!("Parsed {} valid BED entries", entries.len());
    
    // Group entries by chromosome for parallel processing
    let mut chrom_groups: BTreeMap<String, Vec<BedEntry>> = BTreeMap::new();
    for entry in entries {
        chrom_groups.entry(entry.chrom.clone()).or_default().push(entry);
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
                    current.merge_with(next);
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
    let total_original = lines.len();
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
    
    println!("Merged {} entries into {} entries", total_original, total_merged);
    println!("Output written to: {}", output_file);
    
    Ok(())
}