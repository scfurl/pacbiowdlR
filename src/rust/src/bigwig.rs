use itertools::Itertools;
use ndarray::Array2;
use rayon::prelude::*;
use std::{collections::HashSet, path::Path};
use anyhow::*;
use noodles_bed as bed;

pub fn read_bed_with_group(path: &std::path::Path) -> Result<Vec<Region>> {
    let group = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("ranges")
        .to_string();

    let mut reader = bed::io::Reader::from_path(path)?;
    let mut out = Vec::new();

    for result in reader.records() {
        let rec = result?;
        let chrom = rec.reference_sequence_name().to_string();
        let start = rec.start_position().map(|p| u32::try_from(usize::from(p) - 1).unwrap()).unwrap_or(0);
        let end   = rec.end_position().map(|p| u32::try_from(usize::from(p)).unwrap()).unwrap_or(start+1);
        let name  = rec.name().map(|s| s.to_string());

        if end > start {
            out.push(Region { chrom, start, end, name, group: group.clone() });
        }
    }
    Ok(out)
}

/// A simple genomic interval; strand is optional for this use-case.
#[derive(Clone, Debug, serde::Serialize)]
pub struct Region {
    pub chrom: String,
    pub start: u32, // 0-based, half-open
    pub end:   u32,
    pub name:  Option<String>,
    pub group: String, // derived from the BED filename (basename)
}

/// Parameters mirroring the R function (bigWig-only subset)
pub struct Params {
    pub style: Style,               // "percentOfRegion" | "point"
    pub n_windows: usize,           // windows across region (or across flank window for point)
    pub bin_size: usize,            // only used for style=point post-binning
    pub distance_around: Option<f64>, // percent of region length, only used for percentOfRegion
    pub distance_up: Option<u32>,   // used for point
    pub distance_down: Option<u32>, // used for point
}

pub enum Style {
    PercentOfRegion,
    Point,
}

mod bigwig {
    use anyhow::*;

    // Adapt these to your actual bigtools API:
    pub struct Bw {
        inner: bigtools::bigwig::BigWigRead, // adjust if different
    }

    pub fn open(path: &std::path::Path) -> Result<Bw> {
        let inner = bigtools::bigwig::BigWigRead::open(path)
            .with_context(|| format!("open bigWig: {}", path.display()))?;
        Ok(Bw { inner })
    }

    /// Return set of chrom names contained in the bigWig.
    pub fn chrom_names(bw: &Bw) -> std::collections::HashSet<String> {
        bw.inner.chroms().keys().cloned().collect()
    }

    /// Iterate (start,end,value) triples over [start,end) on chrom.
    /// If your API has a “summary over N bins”, prefer that! You can then drop window-level iteration below.
    pub fn iter_values<'a>(
        bw: &'a Bw,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Box<dyn Iterator<Item = (u32,u32,f32)> + 'a>> {
        let it = bw.inner.interval_iter(chrom, start as u64, end as u64)
            .with_context(|| format!("interval_iter {}:{}-{}", chrom, start, end))?
            .map(|r| {
                let r = r.unwrap(); // or handle properly
                (r.start as u32, r.end as u32, r.value as f32)
            });
        Ok(Box::new(it))
    }
}


pub struct ChipProfile {
    /// Per bigWig: matrix rows = regions (concatenated in provided order), cols = windows
    pub per_bw_matrix: Vec<(String, Array2<f32>)>,
    /// Row metadata in same row order (for all groups concatenated):
    pub rows: Vec<Region>,
    /// Column edges for reference (same for all bigWigs if style is consistent)
    pub bin_edges: Vec<u32>, // last seen (varies per row if variable-length windows—see note below)
}

/// NOTE on bin edges:
/// - For `percentOfRegion`, edges differ per row (if region lengths differ). If you need per-row edges,
///   store Vec<Vec<u32>> instead. For speed/memory, many pipelines just store the matrix.

pub fn chipprofile_from_bigwigs(
    bigwigs: &[&Path],
    bed_files: &[&Path],
    params: &Params,
) -> anyhow::Result<ChipProfile> {
    // 1) Load bigWigs
    let mut bws = Vec::new();
    let mut bw_chrom_sets = Vec::new();
    for p in bigwigs {
        let bw = bigwig::open(p)?;
        let chroms = bigwig::chrom_names(&bw);
        bw_chrom_sets.push(chroms);
        bws.push((p.file_name().unwrap().to_string_lossy().to_string(), bw));
    }

    // 2) Find common chroms across bigWigs
    let common_chroms: HashSet<String> = bw_chrom_sets
        .iter()
        .skip(1)
        .fold(bw_chrom_sets[0].clone(), |acc, set| {
            acc.intersection(set).cloned().collect()
        });

    // Decide “style” for each bw (UCSC vs NCBI) once
    let bw_wants_ucsc: Vec<bool> = bws.iter()
        .map(|(_,bw)| is_ucsc_style(&bigwig::chrom_names(bw)))
        .collect();

    // 3) Read regions and filter to common chroms; assign group labels by bed filename
    let mut all_regions: Vec<Region> = Vec::new();
    for bed in bed_files {
        let mut rs = read_bed_with_group(bed)?;
        rs.retain(|r| common_chroms.contains(&r.chrom) || common_chroms.contains(&format!("chr{}", r.chrom)));
        all_regions.extend(rs);
    }

    // 4) For each bigWig, compute matrix
    let per_bw_matrix: Vec<(String, Array2<f32>)> = bws.into_par_iter().enumerate().map(|(i,(bw_name,bw))| {
        let want_ucsc = bw_wants_ucsc[i];

        // Build matrix [n_regions x n_windows]
        let n_regions = all_regions.len();
        let n_windows = params.n_windows;
        let mut mat = Array2::<f32>::zeros((n_regions, n_windows));

        for (row, r) in all_regions.iter().enumerate() {
            let chrom = to_style(&r.chrom, want_ucsc);

            let (win_start, win_end, edges) = match params.style {
                Style::PercentOfRegion => {
                    let (s,e) = match params.distance_around {
                        Some(p) => region_with_percent_flanks(r, p, None),
                        None    => (r.start, r.end),
                    };
                    let edges = linspace_u32(s, e, n_windows);
                    (s, e, edges)
                }
                Style::Point => {
                    // Use the center of the region as the “point”
                    let center = r.start / 2 + r.end / 2;
                    let up  = params.distance_up.unwrap_or(1000);
                    let dn  = params.distance_down.unwrap_or(1000);
                    let s = center.saturating_sub(up);
                    let e = center.saturating_add(dn);
                    let edges = linspace_u32(s, e, n_windows);
                    (s, e, edges)
                }
            };

            // Fill columns by coverage-weighted mean per bin
            for col in 0..n_windows {
                let s = edges[col];
                let e = edges[col+1].max(s+1);
                let m = mean_signal_in_window(&bw, &chrom, s, e)
                    .unwrap_or(f32::NAN);
                mat[(row, col)] = m;
            }
        }

        // Optional: bin columns for style=point
        let mat = if matches!(params.style, Style::Point) && params.bin_size > 1 {
            let bins = n_windows / params.bin_size;
            let mut binned = Array2::<f32>::zeros((n_regions, bins));
            for row in 0..n_regions {
                for b in 0..bins {
                    let c0 = b * params.bin_size;
                    let c1 = c0 + params.bin_size;
                    let mut acc = 0.0f32;
                    for c in c0..c1 {
                        acc += mat[(row, c)];
                    }
                    binned[(row, b)] = acc / (params.bin_size as f32);
                }
            }
            binned
        } else {
            mat
        };

        Ok::<_, anyhow::Error>((bw_name, mat))
    }).collect::<Result<Vec<_>>>()?;

    Ok(ChipProfile {
        per_bw_matrix,
        rows: all_regions,
        bin_edges: vec![], // omit to save memory; or store per-row edges if you need them
    })
}


fn mean_signal_in_window(
    bw: &bigwig::Bw,
    chrom: &str,
    w_start: u32,
    w_end: u32,
) -> anyhow::Result<f32> {
    if w_end <= w_start { return Ok(f32::NAN); }
    let mut num = 0.0f64;
    let mut den = 0.0f64;

    for (s,e,v) in bigwig::iter_values(bw, chrom, w_start, w_end)? {
        let ss = s.max(w_start);
        let ee = e.min(w_end);
        if ee > ss {
            let w = (ee - ss) as f64;
            num += w * (v as f64);
            den += w;
        }
    }
    if den > 0.0 { Ok((num/den) as f32) } else { Ok(0.0) }
}

#[inline]
fn linspace_u32(start: u32, end: u32, n: usize) -> Vec<u32> {
    // n edges (not bins) OR n+1 edges? We want n windows → n+1 edges.
    let n_edges = n + 1;
    (0..n_edges)
        .map(|i| {
            let t = i as f64 / (n as f64);
            let x = (start as f64) * (1.0 - t) + (end as f64) * t;
            x.floor().max(0.0) as u32
        })
        .collect()
}

fn region_with_percent_flanks(r: &Region, pct: f64, chrom_len: Option<u32>) -> (u32,u32) {
    let len = (r.end - r.start) as f64;
    let flank = (len * (pct / 100.0)).round() as i64;
    let s = (r.start as i64 - flank).max(0) as u32;
    let e = match chrom_len {
        Some(cl) => (r.end as i64 + flank).min(cl as i64) as u32,
        None => (r.end as i64 + flank).max(r.end as i64) as u32,
    };
    (s, e)
}


fn is_ucsc_style(chroms: &std::collections::HashSet<String>) -> bool {
    chroms.iter().any(|c| c.starts_with("chr"))
}

fn to_style(chrom: &str, want_ucsc: bool) -> String {
    let has_chr = chrom.starts_with("chr");
    match (has_chr, want_ucsc) {
        (true,  true)  => chrom.to_string(),
        (false, false) => chrom.to_string(),
        (false, true)  => {
            // don't add "chr" to odd contigs that already start with GL/KN/etc.
            if chrom.starts_with("GL") || chrom.starts_with("KI") { chrom.to_string() }
            else { format!("chr{}", chrom) }
        }
        (true, false)  => chrom.trim_start_matches("chr").to_string(),
    }
}
