extern crate rayon;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use rayon::prelude::*;
use std::collections::BTreeMap;

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

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    
    if args.len() < 3 || args.len() > 4 {
        eprintln!("Usage: {} <input_file> <output_file> [threads]", args[0]);
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
    
    // Read the entire file into memory for parallel processing
    let file = File::open(input_file)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;
    
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
    Ok(())
}