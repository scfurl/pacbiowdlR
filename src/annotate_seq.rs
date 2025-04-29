use std::collections::HashMap;
use std::error::Error;
use std::str;
use clap::Parser;
use rust_htslib::bam::{Reader, Read};
use rust_htslib::bam::record::Aux;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author, version, about = "Annotate m6A (mA) and m5C (mC) in reads using MM/ML tags")]
struct Args {
    /// Input BAM/CRAM file
    #[arg(short, long)]
    bam: String,

    /// Quality threshold (Phred-scale) for calling a modification
    #[arg(short, long, default_value_t = 18)]
    threshold: u32,

    /// Only process this specific read name (optional)
    #[arg(short = 'n', long)]
    read_name: Option<String>,
}

/// Parse one MM sub-tag like "A+a.,4,9,20" or "C+m,10,50"
/// into a vector of absolute read-positions (0-based).
fn decode_mm(mm_sub: &str) -> Vec<usize> {
    let parts = mm_sub.split(',');
    // header is e.g. "A+a." or "C+m,"
    let offsets: Vec<usize> = parts
        .skip(1)
        .filter_map(|s| s.parse().ok())
        .collect();

    if offsets.is_empty() {
        return vec![];
    }

    // convert run-lengths into absolute positions
    let mut pos = offsets[0];
    let mut result = vec![pos];
    for &off in &offsets[1..] {
        pos = pos + 1 + off;
        result.push(pos);
    }
    result
}

/// Annotate a BAM/CRAM with "mA" and "mC" for calls above `threshold`.
fn annotate_mods(path: &str, threshold: u32, filter_name: &Option<String>) -> Result<(), Box<dyn Error>> {
    let mut bam = Reader::from_path(path)?;

    for record in bam.records() {
        let rec = record?;
        let name = str::from_utf8(rec.qname())?;

        // If a specific read name is provided, skip others
        if let Some(ref target) = *filter_name {
            if name != target {
                continue;
            }
        }

        // raw sequence (A/C/G/T)
        let seq = rec.seq().as_bytes();

        // pull tags
        let mm_tag = match rec.aux(b"MM") {
            Ok(Aux::String(s)) => s,
            _ => continue,
        };
        let ml_tag = match rec.aux(b"ML") {
            Ok(Aux::ArrayU32(v)) => v,
            _ => continue,
        };

        // map read-position -> annotation
        let mut ann: HashMap<usize, &str> = HashMap::new();
        let mut ml_iter = ml_tag.iter();

        for sub in mm_tag.split(';') {
            let (annot, _base_code) = if sub.starts_with("A+a") {
                ("mA", b'A')
            } else if sub.starts_with("C+m") {
                ("mC", b'C')
            } else {
                continue;
            };

            let positions = decode_mm(sub);
            for pos in positions {
                if pos >= seq.len() {
                    break;
                }
                if let Some(q) = ml_iter.next() {
                    if q >= threshold {
                        ann.insert(pos, annot);
                    }
                }
            }
        }

        // build annotated sequence
        let mut out = String::with_capacity(seq.len() * 2);
        for (i, &b) in seq.iter().enumerate() {
            if let Some(&a) = ann.get(&i) {
                out.push_str(a);
            } else {
                out.push(b as char);
            }
        }

        println!("{}\t{}", name, out);
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    annotate_mods(&args.bam, args.threshold, &args.read_name)?;
    Ok(())
}
