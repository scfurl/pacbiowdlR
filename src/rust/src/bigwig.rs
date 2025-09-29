#![allow(dead_code)]
// use ndarray::Array2;
// use rayon::prelude::*;
// use std::{collections::HashSet, path::Path};
use anyhow::*;
use noodles_bed as bed;
// use bigtools::BigWigRead;
use std::result::Result::Ok;
use bigtools::{BigWigRead, Value};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use noodles_bed::feature::record::Strand;
// use std::char::Char;
// 
#[derive(Clone, Debug)]
pub struct Region {
    pub chrom: String,   // "1" or "chr1"
    pub start: u32,      // 0-based, half-open
    pub end: u32,
    pub name: Option<String>,
    pub strand: Option<char>, // '+' or '-' or None
    pub coverage_bins: Option<Vec<f32>>, // per-bin aggregated values
}

#[derive(Clone, Copy, Debug)]
pub enum Reduce {
    Mean,
    Sum,
    Max,
    Min,
}

#[derive(Clone, Copy, Debug)]
pub enum Missing {
    NaN,
    Zero,
}

/// Parallel binner: fixed bin width in bp, with Reduce option.
/// Each thread opens its own BigWig reader for thread-safety.
pub fn read_bigwig(
    path: &str,
    mut regions: Vec<Region>,
    bin_width: u32,
    reduce: Reduce,
    num_threads: Option<usize>, // None => rayon default
    missing: Missing,           // <— NEW: NaN (deepTools default) or Zero (like --missingDataAsZero)
) -> anyhow::Result<Vec<Region>> {
    assert!(bin_width >= 1, "bin_width must be >= 1");

    // Inspect chrom style once (as you already do)
    let header_reader = BigWigRead::open_file(path)?;
    let chroms = header_reader.chroms();
    if chroms.is_empty() {
        return Err(anyhow!("BigWig file has no chromosomes"));
    }
    let bw_is_ucsc = chroms[1].name.starts_with("chr");

    // Build a LOCAL pool and run ALL parallel work inside install()
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or_else(|| rayon::current_num_threads()))
        .build()?;

    pool.install(|| {
        process_parallel(path, bw_is_ucsc, &mut regions, bin_width, reduce, missing)
    })?;

    Ok(regions)
}

fn process_parallel(
    path: &str,
    bw_is_ucsc: bool,
    regions: &mut [Region],
    bin_width: u32,
    reduce: Reduce,
    missing: Missing,            // <— NEW
) -> anyhow::Result<()> {
    let chunk_size = 512;

    regions.par_chunks_mut(chunk_size).try_for_each(|chunk| -> anyhow::Result<()> {
        let mut reader = BigWigRead::open_file(path)?;

        let map_chrom = |c: &str| -> String {
            let has_chr = c.starts_with("chr");
            match (has_chr, bw_is_ucsc) {
                (true,  true)  => c.to_string(),
                (false, false) => c.to_string(),
                (false, true)  => format!("chr{}", c),
                (true,  false) => c.trim_start_matches("chr").to_string(),
            }
        };

        for region in chunk.iter_mut() {
            if region.end <= region.start {
                region.coverage_bins = Some(vec![]);
                continue;
            }

            let chrom = map_chrom(&region.chrom);
            let len = region.end - region.start;
            let bins = ((len + bin_width - 1) / bin_width) as usize;
            
            // full span of each bin (for last partial bin)
            let mut span = vec![0.0f64; bins];
            for b in 0..bins {
                let b_start = region.start + (b as u32) * bin_width;
                let b_end   = (b_start + bin_width).min(region.end);
                span[b] = (b_end - b_start) as f64;
            }
            
            // Accumulators
            let mut num = match reduce {
                Reduce::Mean | Reduce::Sum => vec![0.0f64; bins],
                _ => Vec::new(),
            };
            let mut den = match reduce {
                Reduce::Mean => vec![0.0f64; bins],
                _ => Vec::new(),
            };
            let mut covered_sum = match reduce {
                Reduce::Sum => vec![false; bins], // mark if any bp in bin was covered
                _ => Vec::new(),
            };
            let mut agg = match reduce {
                Reduce::Max => vec![f32::NEG_INFINITY; bins],
                Reduce::Min => vec![f32::INFINITY; bins],
                _ => Vec::new(),
            };
            let mut seen = match reduce {
                Reduce::Max | Reduce::Min => vec![false; bins],
                _ => Vec::new(),
            };

            // Iterate sparse bigWig intervals overlapping this region
            for iv in reader.get_interval(&chrom, region.start, region.end)? {
                let Value { start, end, value } = iv?;
                let mut cur  = (start as u32).max(region.start);
                let stop     = (end   as u32).min(region.end);
                if stop <= cur { continue; }

                while cur < stop {
                    let rel    = cur - region.start;
                    let b      = (rel / bin_width) as usize;
                    let b_end  = (region.start + ((b as u32 + 1) * bin_width)).min(region.end);
                    let next   = b_end.min(stop);
                    if next > cur {
                        let w = (next - cur) as f64;
                        match reduce {
                            Reduce::Mean => { num[b] += w * (value as f64); den[b] += w; }
                            Reduce::Sum  => { num[b] += w * (value as f64); covered_sum[b] = true; }
                            Reduce::Max  => { agg[b] = agg[b].max(value);   seen[b] = true; }
                            Reduce::Min  => { agg[b] = agg[b].min(value);   seen[b] = true; }
                        }
                    }
                    cur = next;
                }
            }

            // Finalize with deepTools-like missing handling
            let default = match missing { Missing::NaN => f32::NAN, Missing::Zero => 0.0 };
            let mut out = vec![default; bins];

            match reduce {
                // Reduce::Mean => {
                //     for b in 0..bins {
                //         if den[b] > 0.0 {
                //             out[b] = (num[b] / den[b]) as f32;
                //         }
                //     }
                // }
                Reduce::Mean => {
                    for b in 0..bins {
                        out[b] = match missing {
                            Missing::NaN  => {
                                if den[b] > 0.0 { (num[b] / den[b]) as f32 } else { f32::NAN }
                            }
                            Missing::Zero => {
                                // divide by full bin length (deepTools behavior)
                                if span[b] > 0.0 { (num[b] / span[b]) as f32 } else { 0.0 }
                            }
                        };
                    }
                }
                Reduce::Sum => {
                    for b in 0..bins {
                        if covered_sum[b] {
                            out[b] = num[b] as f32; // sum(value * covered_bp)
                        }
                    }
                }
                Reduce::Max => {
                    for b in 0..bins {
                        if seen[b] {
                            out[b] = agg[b];
                        }
                    }
                }
                Reduce::Min => {
                    for b in 0..bins {
                        if seen[b] {
                            out[b] = agg[b];
                        }
                    }
                }
            }
            if matches!(region.strand, Some('-')) {
                out.reverse();
            }
            region.coverage_bins = Some(out);
        }

        Ok(())
    })?;

    Ok(())
}

// #![allow(dead_code)]
// // use ndarray::Array2;
// // use rayon::prelude::*;
// // use std::{collections::HashSet, path::Path};
// use anyhow::*;
// use noodles_bed as bed;
// // use bigtools::BigWigRead;
// use std::result::Result::Ok;
// use bigtools::{BigWigRead, Value};
// use rayon::ThreadPoolBuilder;
// use rayon::prelude::*;


pub fn read_bed_with_group(path: &std::path::Path) -> Result<Vec<Region>> {
    let _group = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("ranges")
        .to_string();

    let mut reader = bed::io::reader::Builder::<6>.build_from_path(path)?;
    let mut out = Vec::new();

    let mut record = bed::Record::default();
    let valid_chroms: std::collections::HashSet<&'static str> = [
        "1","2","3","4","5","6","7","8","9","10",
        "11","12","13","14","15","16","17","18","19",
        "20","21","22","X","Y","M",
        "chr1","chr2","chr3","chr4","chr5","chr6",
        "chr7","chr8","chr9","chr10","chr11","chr12",
        "chr13","chr14","chr15","chr16","chr17",
        "chr18","chr19","chr20","chr21","chr22",
        "chrX","chrY","chrM",
    ].into_iter().collect();


    while reader.read_record(&mut record)? != 0 {
        let chrom = record.reference_sequence_name().to_string();
        let start = record.feature_start().map(|p| u32::try_from(usize::from(p) - 1).unwrap()).unwrap_or(0);
        let end = record.feature_end().map(|p| u32::try_from(usize::from(p.unwrap())).unwrap()).unwrap_or(start + 1); // PROBLEMATIC unwrap
        let name = record.name().map(|s| s.to_string());
            let strand = record.strand().transpose().map(|s| {
        match s {
            Ok(Strand::Forward) => '+',
            Ok(Strand::Reverse) => '-',
            // Ok(Strand::Unknown) => '.',
            Err(_) => '.',
        }
    });
        if end > start && valid_chroms.contains(chrom.as_str()) {
            out.push(Region { chrom, start, end, name, strand, coverage_bins: None });
        }
    }
    // eprintln!("Read {} regions from {}", out.len(), path.display());
    Ok(out)
}



// #[derive(Clone, Debug)]
// pub struct Region {
//     pub chrom: String,   // "1" or "chr1"
//     pub start: u32,      // 0-based, half-open
//     pub end: u32,
//     pub name: Option<String>,
//     pub coverage_bins: Option<Vec<f32>>, // per-bin aggregated values
// }

// #[derive(Clone, Copy, Debug)]
// pub enum Reduce {
//     Mean,  // length-weighted
//     Sum,   // integral value * covered_bp
//     Max,   // max over covered sub-intervals
//     Min,   // min over covered sub-intervals
// }

// /// Parallel binner: fixed bin width in bp, with Reduce option.
// /// Each thread opens its own BigWig reader for thread-safety.
// pub fn read_bigwig(
//     path: &str,
//     mut regions: Vec<Region>,
//     bin_width: u32,
//     reduce: Reduce,
//     num_threads: Option<usize>, // None => rayon default
// ) -> anyhow::Result<Vec<Region>> {
//     assert!(bin_width >= 1, "bin_width must be >= 1");

//     // Open once just to detect chromosome style
//     let header_reader = BigWigRead::open_file(path)?;
//     let chroms = header_reader.chroms();
//     let mut bw_is_ucsc = true; // assume UCSC-style for now
//     if chroms.is_empty() {
//         return Err(anyhow!("BigWig file has no chromosomes"));
//     }
//     if chroms[1].name.starts_with("chr") {
//         eprintln!("BigWig has UCSC-style");
//         bw_is_ucsc = true;
//     } else {
//         eprintln!("BigWig has non-UCSC style");
//         bw_is_ucsc = false;
//     }

//     // let bw_is_ucsc = header_reader.chroms().keys().any(|c| c.starts_with("chr"));
//     // drop(header_reader);



//     let _map_chrom = move |c: &str| -> String {
//         let has_chr = c.starts_with("chr");
//         match (has_chr, true) {
//             (true,  true)  => c.to_string(),
//             (false, false) => c.to_string(),
//             (false, true)  => format!("chr{}", c),
//             (true,  false) => c.trim_start_matches("chr").to_string(),
//         }
//     };

//     // Build a LOCAL pool and run ALL parallel work inside install()
//     let pool = ThreadPoolBuilder::new().num_threads(num_threads.unwrap()).build()?;
//     pool.install(|| {
//         process_parallel(path, bw_is_ucsc, &mut regions, bin_width, reduce)
//     })?;

//     // process_parallel(path, bw_is_ucsc, &mut regions, bin_width, reduce)?;

//     Ok(regions)
// }


// fn process_parallel(
//     path: &str,
//     bw_is_ucsc: bool,
//     regions: &mut [Region],
//     bin_width: u32,
//     reduce: Reduce,
// ) -> anyhow::Result<()> {
//     let chunk_size = 512;

//     regions.par_chunks_mut(chunk_size).try_for_each(|chunk| -> anyhow::Result<()> {
//         let mut reader = BigWigRead::open_file(path)?;

//         // map chrom helper using the bool captured above
//         let map_chrom = |c: &str| -> String {
//             let has_chr = c.starts_with("chr");
//             match (has_chr, bw_is_ucsc) {
//                 (true,  true)  => c.to_string(),
//                 (false, false) => c.to_string(),
//                 (false, true)  => format!("chr{}", c),
//                 (true,  false) => c.trim_start_matches("chr").to_string(),
//             }
//         };

//         for (_i, region) in chunk.iter_mut().enumerate() {
//             if region.end <= region.start {
//                 region.coverage_bins = Some(vec![0.0]);
//                 continue;
//             }

//             let chrom = map_chrom(&region.chrom);
//             let len = region.end - region.start;
//             let bins = ((len + bin_width - 1) / bin_width) as usize;

//             let mut num = match reduce {
//                 Reduce::Mean | Reduce::Sum => vec![0.0f64; bins],
//                 Reduce::Max | Reduce::Min => Vec::new(),
//             };
//             let mut den = match reduce {
//                 Reduce::Mean => vec![0.0f64; bins],
//                 _ => Vec::new(),
//             };
//             let mut agg = match reduce {
//                 Reduce::Max => vec![f32::NEG_INFINITY; bins],
//                 Reduce::Min => vec![f32::INFINITY; bins],
//                 _ => Vec::new(),
//             };
//             let mut seen = match reduce {
//                 Reduce::Max | Reduce::Min => vec![false; bins],
//                 _ => Vec::new(),
//             };

//             for iv in reader.get_interval(&chrom, region.start as u32, region.end as u32)? {
//                 let Value { start, end, value } = iv?;
//                 let mut cur = (start as u32).max(region.start);
//                 let stop = (end as u32).min(region.end);
//                 if stop <= cur { continue; }

//                 while cur < stop {
//                     let rel = cur - region.start;
//                     let b = (rel / bin_width) as usize;
//                     let bin_end = (region.start + ((b as u32 + 1) * bin_width)).min(region.end);
//                     let next = bin_end.min(stop);
//                     if next > cur {
//                         let w = (next - cur) as f64;
//                         match reduce {
//                             Reduce::Mean => { num[b] += w * (value as f64); den[b] += w; }
//                             Reduce::Sum  => { num[b] += w * (value as f64); }
//                             Reduce::Max  => { agg[b] = agg[b].max(value); seen[b] = true; }
//                             Reduce::Min  => { agg[b] = agg[b].min(value); seen[b] = true; }
//                         }
//                     }
//                     cur = next;
//                 }
//             }

//             let mut out = vec![0.0f32; bins];
//             match reduce {
//                 Reduce::Mean => for b in 0..bins { out[b] = if den[b] > 0.0 { (num[b]/den[b]) as f32 } else { 0.0 }; },
//                 Reduce::Sum  => for b in 0..bins { out[b] = num[b] as f32; },
//                 Reduce::Max  => for b in 0..bins { out[b] = if seen[b] { agg[b] } else { 0.0 }; },
//                 Reduce::Min  => for b in 0..bins { out[b] = if seen[b] { agg[b] } else { 0.0 }; },
//             }

//             region.coverage_bins = Some(out);
//         }
//         Ok(())
//     })?;

//     Ok(())
// }



