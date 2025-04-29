/**
samtools view -s 0.00001 -b /Users/sfurlan/Desktop/CHLA10.GRCh38.haplotagged.bam | samtools view -b -h - -o /Users/sfurlan/Desktop/CHLA10.test.bam
./annotate_seq -b /Users/sfurlan/Desktop/CHLA10.GRCh38.haplotagged.bam -n m84265_250417_230407_s1/173412086/ccs
 **/
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
     threshold: u8,
 
     /// Only process this specific read name (optional)
     #[arg(short = 'n', long)]
     read_name: Option<String>,
 }
 
 /// Decode MM sub-tag into absolute read positions of the specified base
 fn decode_mm_offsets(mm_sub: &str, base: u8, seq: &[u8]) -> Vec<usize> {
     // Split off the header (e.g. "A+a.")
     let mut parts = mm_sub.split(',');
     let _header = parts.next();
     // Parse comma-separated offsets
     let offsets: Vec<usize> = parts.filter_map(|s| s.parse().ok()).collect();
 
     // Collect positions of the target base in the read
     let base_positions: Vec<usize> = seq.iter()
         .enumerate()
         .filter(|&(_, &b)| b == base)
         .map(|(i, _)| i)
         .collect();
 
     let mut result = Vec::new();
     // Track index into base_positions
     if offsets.is_empty() {
         return result;
     }
     // First modification occurs at the "offsets[0]"-th occurrence of `base`
     let mut idx = offsets[0];
     if let Some(&pos) = base_positions.get(idx) {
         result.push(pos);
     } else {
         return result;
     }
     // Subsequent offsets are counts of intervening base occurrences
     for &off in &offsets[1..] {
         idx = idx + off + 1;
         if let Some(&pos) = base_positions.get(idx) {
             result.push(pos);
         } else {
             break;
         }
     }
     result
 }
 
 /// Annotate a BAM/CRAM with "mA" and "mC" for calls above `threshold`.
 fn annotate_mods(path: &str, threshold: u8, filter_name: &Option<String>) -> Result<(), Box<dyn Error>> {
     let mut bam = Reader::from_path(path)?;
     for record in bam.records() {
         let rec = record?;
         let name = str::from_utf8(rec.qname())?;
 
         // If a specific read name is provided, skip others
         if let Some(target) = filter_name {
             if &name != target {
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
         // ML is stored as an array of u8
         let ml_tag = match rec.aux(b"ML") {
             Ok(Aux::ArrayU8(v)) => v,
             _ => continue,
         };
 
         // map read-position -> annotation
         let mut ann: HashMap<usize, &str> = HashMap::new();
         let mut ml_iter = ml_tag.iter();
 
         for sub in mm_tag.split(';') {
             let (annot, base): (&str, u8) = if sub.starts_with("A+a") {
                 ("mA", b'A')
             } else if sub.starts_with("C+m") {
                 ("mC", b'C')
             } else {
                 continue;
             };
 
             // calculate absolute positions of modifications
             let positions = decode_mm_offsets(sub, base, &seq);
             for pos in positions {
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
 