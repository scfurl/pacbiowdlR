
extern crate rayon;
extern crate flate2;
extern crate clap;
extern crate rust_htslib;
extern crate thiserror;
extern crate log;
extern crate env_logger;


// src/main.rs

use clap::Parser;
use log::{info, warn};
use rayon::prelude::*;
use rust_htslib::tbx::{Reader, Read};
use rust_htslib::errors::Error as TbxError;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};
use std::path::PathBuf;
use std::sync::Arc;
use thiserror::Error;

/// CLI for merging adjacent -1.0 BED intervals from a Tabix-indexed BED file
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Tabix-indexed bgzipped BED file (.bed.gz)
    #[clap(value_parser)]
    input: PathBuf,

    /// Output file to write merged BED entries
    #[clap(value_parser)]
    output: PathBuf,

    /// Chrom sizes file (tab-delimited: chrom, length)
    #[clap(long, value_parser)]
    chrom_sizes: PathBuf,

    /// Number of threads to use (default: all available)
    #[clap(short, long, value_parser)]
    threads: Option<usize>,
}

/// Unified CLI error type
#[derive(Error, Debug)]
enum CliError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("Tabix error: {0}")]
    Tbx(#[from] TbxError),
    #[error("Parse error: {0}")]
    Parse(String),
}

/// Basic BED entry struct
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct BedEntry {
    chrom: String,
    start: u64,
    end: u64,
    value: f64,
}

impl BedEntry {
    fn new(chrom: String, start: u64, end: u64, value: f64) -> Self {
        Self { chrom, start, end, value }
    }
}

impl From<std::num::ParseIntError> for CliError {
    fn from(e: std::num::ParseIntError) -> Self { CliError::Parse(e.to_string()) }
}
impl From<std::num::ParseFloatError> for CliError {
    fn from(e: std::num::ParseFloatError) -> Self { CliError::Parse(e.to_string()) }
}

fn main() -> Result<(), CliError> {
    // Initialize logger at INFO level
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    // Load chrom sizes (required)
let file = File::open(&args.chrom_sizes)?;
let reader = BufReader::new(file);
let sizes_map: HashMap<String, u64> = reader.lines()
    .map(|l| {
        let line = l.map_err(CliError::Io)?;
        let mut parts = line.split_whitespace();
        let chrom = parts.next()
            .ok_or_else(|| CliError::Parse("Missing chrom in sizes file".into()))?;
        let len = parts.next()
            .ok_or_else(|| CliError::Parse(format!("Missing length for {}", chrom)))?
            .parse::<u64>()?;
        Ok::<(std::string::String, u64), CliError>((chrom.to_string(), len))
    })
    .collect::<Result<_, _>>()?;
let sizes = Arc::new(sizes_map);

    // Check for accompanying .tbi index
    let idx = format!("{}.tbi", args.input.display());
    if !std::path::Path::new(&idx).exists() {
        return Err(CliError::Io(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Tabix index not found: {}", idx),
        )));
    }

    // Configure Rayon threads
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to set thread count");
    }
    info!("Using {} threads", rayon::current_num_threads());

    // Open reader and list contigs
    let tbx = Reader::from_path(&args.input)?;
    let contigs: Vec<String> = tbx.seqnames().into_iter().map(String::from).collect();
    info!("Found {} contigs", contigs.len());

    // Process each contig in parallel
    let mut results: Vec<(String, Vec<BedEntry>)> = contigs
        .into_par_iter()
        .filter_map(|chrom| {
            // Determine fetch end position from sizes or default to MAX
            let end_pos = sizes.get(&chrom).cloned().unwrap_or(std::u64::MAX);
            // Open fresh reader
            let merged_res: Result<Vec<BedEntry>, CliError> = Reader::from_path(&args.input)
                .map_err(CliError::Tbx)
                .and_then(|mut reader| {
                    // Resolve rid for this contig
                    let rid = reader
                        .tid(&chrom)
                        .map_err(|e| CliError::Tbx(e))?;
                    // Fetch with correct end
                    reader.fetch(rid, 0, end_pos)
                        .map_err(CliError::Tbx)?;
                    // Merge intervals
                    merge_intervals(&mut reader, &chrom)
                });

            match merged_res {
                Ok(merged) => {
                    info!("{}: merged {} intervals", chrom, merged.len());
                    Some((chrom, merged))
                }
                Err(e) => {
                    warn!("Skipping {}: {}", chrom, e);
                    None
                }
            }
        })
        .collect();

    // Sort for deterministic output
    results.sort_by(|(a, _), (b, _)| a.cmp(b));

    // Write out merged intervals
    let mut out = BufWriter::new(std::fs::File::create(&args.output)?);
    let mut total = 0;
    for (chrom, intervals) in results {
        for iv in intervals {
            let val = if (iv.value + 1.0).abs() < f64::EPSILON { 0.0 } else { iv.value };
            writeln!(out, "{}\t{}\t{}\t{}", chrom, iv.start, iv.end, val)?;
            total += 1;
        }
    }
    info!("Wrote {} merged entries to {}", total, args.output.display());
    Ok(())
}

/// Merge adjacent -1.0 intervals
fn merge_intervals(reader: &mut Reader, chrom: &str) -> Result<Vec<BedEntry>, CliError> {
    let mut merged = Vec::new();
    let mut current: Option<BedEntry> = None;
    for rec in reader.records() {
        let raw = rec?;
        let cols: Vec<&[u8]> = raw.split(|&b| b == b'\t').collect();
        if cols.len() < 6 { continue; }
        let start: u64 = std::str::from_utf8(cols[1]).unwrap().parse().unwrap();
        let endp:  u64 = std::str::from_utf8(cols[2]).unwrap().parse().unwrap();
        let value: f64 = std::str::from_utf8(cols[5]).unwrap().parse().unwrap();
        let entry = BedEntry::new(chrom.to_string(), start, endp, value);
        if let Some(mut cur) = current {
            if cur.end == entry.start && cur.value == -1.0 && entry.value == -1.0 {
                cur.end = entry.end;
                current = Some(cur);
            } else {
                merged.push(cur);
                current = Some(entry);
            }
        } else {
            current = Some(entry);
        }
    }
    if let Some(c) = current { merged.push(c); }
    Ok(merged)
}
