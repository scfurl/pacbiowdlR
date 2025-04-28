extern crate rayon;
extern crate flate2;
extern crate clap;
extern crate rust_htslib;
extern crate thiserror;
extern crate log;
extern crate env_logger;

// use std::fmt::Error;
use std::path::PathBuf;
use std::io::{self, Write, BufWriter};
use clap::Parser;
use rayon::prelude::*;
use rust_htslib::tbx::{Reader, Read};
use rust_htslib::errors::Error as TbxError;
use thiserror::Error;
use log::{info, warn};

/// CLI for merging -1.0 BED intervals from a Tabix-indexed BED file
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Tabix-indexed bgzipped BED file (.bed.gz)
    #[clap(value_parser)]
    input: PathBuf,

    /// Output file to write merged BED entries
    #[clap(value_parser)]
    output: PathBuf,

    /// Number of threads to use (default: all available)
    #[clap(short, long, value_parser)]
    threads: Option<usize>,
}

#[derive(Error, Debug)]
enum CliError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Tabix error: {0}")]
    Tbx(#[from] TbxError),
    #[error("Parse error: {0}")]
    Parse(String),
}

#[derive(Debug, Clone)]
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
    env_logger::init();
    let args = Args::parse();

    // Configure Rayon thread pool
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to build thread pool");
    }
    info!("Using {} threads", rayon::current_num_threads());

    // Validate input file
    if !args.input.exists() {
        return Err(CliError::Io(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Input file '{}' not found", args.input.display()),
        )));
    }
    if args.input.extension().and_then(|s| s.to_str()) != Some("gz") {
        warn!("Input file does not end with .gz");
    }

    // Extract sequence names
    let tbx0 = Reader::from_path(&args.input)?;
    let seqs: Vec<String> = tbx0.seqnames().into_iter().map(String::from).collect();
    info!("Found {} sequences", seqs.len());

    // Parallel fetch & merge per chromosome
    let mut results: Vec<(String, Vec<BedEntry>)> = seqs
        .into_par_iter()
        .map(|chrom| {
            let mut reader = Reader::from_path(&args.input)?;
            // Convert chrom name to tid index
            // let tid = reader.header()
            //     .name2rid(chrom.as_bytes())
            //     .ok_or_else(|| CliError::Parse(format!("Unknown contig: {}", chrom)))?;
            let tid = match reader.tid(&chrom) {
                Ok(tid) => tid,
                Err(_) => panic!("Could not resolve {} to contig ID", chrom),
            };
            let merged = stream_merge_chrom(&mut reader, tid, &chrom)?;
            info!("{}: merged {} entries", chrom, merged.len());
            Ok::<(std::string::String, Vec<BedEntry>), CliError>((chrom, merged))
        })
        .filter_map(Result::ok)
        .collect();

    // Sort for deterministic output
    results.sort_by(|(a, _), (b, _)| a.cmp(b));

    // Write merged output
    let output = std::fs::File::create(&args.output)?;
    let mut wtr = BufWriter::new(output);
    let mut count = 0;
    for (_chrom, entries) in results {
        for e in entries {
            let val = if (e.value - (-1.0)).abs() < f64::EPSILON { 0.0 } else { e.value };
            writeln!(wtr, "{}\t{}\t{}\t{}", e.chrom, e.start, e.end, val)?;
            count += 1;
        }
    }
    info!("Wrote {} merged entries to {}", count, args.output.display());

    Ok(())
}

/// Fetch & merge -1.0 intervals for one chromosome by tid index
fn stream_merge_chrom(
    reader: &mut Reader,
    tid: u64,
    chrom: &str,
) -> Result<Vec<BedEntry>, CliError> {
    reader.fetch(tid as u64, 0, std::u64::MAX)?;
    let mut merged = Vec::new();
    let mut current: Option<BedEntry> = None;

    for record in reader.records() {
        let raw = record?;
        // Split on tab byte
        let fields: Vec<&[u8]> = raw.split(|&b| b == b'\t').collect();
        if fields.len() < 6 { continue; }
        let chrom_str = chrom.to_string();
        let start: u64 = std::str::from_utf8(fields[1]).unwrap().parse()?;
        let end: u64 = std::str::from_utf8(fields[2]).unwrap().parse()?;
        let value: f64 = std::str::from_utf8(fields[5]).unwrap().parse()?;
        let entry = BedEntry::new(chrom_str, start, end, value);

        if let Some(mut cur) = current.take() {
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


    if let Some(cur) = current {
        merged.push(cur);
    }
    Ok(merged)
}
