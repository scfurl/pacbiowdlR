/** 

./motif_finder -i ~/refs/GRCh38/GRCh38.p13.genome.fa -n 50 -s GGAA -m 20
clear
**/

// use std::fs::File;
// use std::io::{self, BufReader};
use std::io;
use std::sync::mpsc;
use clap::Parser;
use bio::io::fasta::{IndexedReader, Index};
use rayon::prelude::*;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author, version, about = "Find fuzzy repeats in a FASTA file (indexed), parallelized by contig, customizable motifs")]
struct Args {
    /// Input FASTA file (must be indexed with samtools faidx)
    #[arg(short, long)]
    input: String,

    /// Min number of consecutive repeats
    #[arg(short = 'n', long, default_value_t = 3)]
    min_repeats: usize,

    /// Max mismatches allowed per motif
    #[arg(short = 'm', long, default_value_t = 1)]
    max_mismatches: usize,

    /// Motif sequences to search (forward strands), comma-delimited (e.g. GGAA,CATG)
    #[arg(short = 's', long = "motifs", value_delimiter = ',')]
    motifs: Vec<String>,
}

/// Compute Hamming distance between two equal-length byte slices
fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Compute the reverse complement of a DNA motif
fn revcomp(motif: &[u8]) -> Vec<u8> {
    motif.iter().rev().map(|&b| match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }).collect()
}

/// Find runs of a given motif allowing up to max_mismatches per motif
/// Returns tuples of (start, end, repeat_count)
fn find_runs(
    seq: &[u8],
    motif: &[u8],
    min_repeats: usize,
    max_mismatches: usize,
) -> Vec<(usize, usize, usize, usize)> {
    let mut runs = Vec::new();
    let mlen = motif.len();
    let seqlen = seq.len();
    let mut i = 0;

    while i + mlen * min_repeats <= seqlen {
        let mut total_mm = 0;
        let mut count = 0;
        let mut j = i;

        while j + mlen <= seqlen {
            let mismatches = hamming(&seq[j..j + mlen], motif);
            if total_mm + mismatches > max_mismatches {
                break;
            }
            total_mm += mismatches;
            count += 1;
            j += mlen;
        }

        if count >= min_repeats {
            runs.push((i, j, count, total_mm));
            i = j; // skip past this run
        } else {
            i += 1;
        }
    }
    runs
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let fai = args.input.clone() + ".fai";


    // Open Index
    let index = Index::from_file(&fai).unwrap();

    // Build motif list with forward and reverse strands
    let mut search_motifs = Vec::new();
    for motif_str in &args.motifs {
        let motif_bytes = motif_str.as_bytes().to_vec();
        search_motifs.push((motif_bytes.clone(), "+"));
        let rc = revcomp(&motif_bytes);
        search_motifs.push((rc, "-"));
    }

    // Print header
    println!("contig\tstrand\tstart\tend\trepeats\tsequence");

    // Channel for thread-safe output
    let (tx, rx) = mpsc::channel();

    // Parallelize scan over contigs from the index
    index.sequences().par_iter().for_each_with(tx.clone(), |tx, contig_rec| {
        let name = &contig_rec.name;
        let len = contig_rec.len;
        let mut reader = IndexedReader::from_file(&args.input).unwrap();
        // Fetch entire contig
        reader.fetch_all(&name).expect("Failed to fetch contig");
        let mut seq = Vec::with_capacity(len as usize);
        reader.read(&mut seq).expect("Failed to read contig");

        // Scan motifs
        for (motif, strand) in &search_motifs {
            let runs = find_runs(
                &seq,
                motif,
                args.min_repeats,
                args.max_mismatches,
            );
            for (start, end, count, total_mm) in runs {
                let matched_seq = &seq[start..end];
                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    name,
                    strand,
                    start + 1,
                    end,
                    count,
                    total_mm,
                    String::from_utf8_lossy(matched_seq)
                );
                tx.send(line).expect("Failed to send result");
            }
        }
    });

    // Close extra sender
    drop(tx);

    // Print results
    for line in rx {
        println!("{}", line);
    }

    Ok(())
}
