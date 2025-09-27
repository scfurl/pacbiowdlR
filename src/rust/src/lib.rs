#![allow(non_snake_case)]
use extendr_api::prelude::*;
//  use rayon::result;

mod bigwig;
use bigwig::Region;
// use std::result::Result::Ok;
// use itertools::Itertools;
// use ndarray::Array2;
// use rayon::prelude::*;
// use std::{collections::HashSet, path::Path};
// use anyhow::*;
use std::result::Result::Ok;
// Read BED file with group handling
/// @export
/// @keywords internal
#[extendr]
pub fn read_bed(bed_path: Robj, bw_path: Robj, bin_size: Robj, method: Robj, nthreads: Robj) -> List {
    // Check inputs
    let bed_path_str = bed_path.as_str().unwrap();
    let bed_path = std::path::Path::new(bed_path_str);
    // check if file exists
    if !bed_path.exists() {
        eprintln!("File {} does not exist", bed_path.display());
        return list!(0);
    }
    eprint!("Reading BED file from {}...\n", bed_path.display());

    let bw_path_str = bw_path.as_str().unwrap();
    let bw_path = std::path::Path::new(bw_path_str);
    // check if file exists
    if !bw_path.exists() {
        eprintln!("File {} does not exist", bw_path.display());
        return list!(0);
    }
    eprint!("Reading bigWig file from {}...\n", bw_path.display());
    let regions = bigwig::read_bed_with_group(bed_path).unwrap();
    let n_regions = regions.len();
    println!("Read {} regions from BED file", n_regions);

    let bin_size = bin_size.as_real().unwrap() as u32;
    let method_str = method.as_str().unwrap();
    let method = match method_str {
        "mean" => bigwig::Reduce::Mean,
        "sum"  => bigwig::Reduce::Sum,
        "max"  => bigwig::Reduce::Max,
        "min"  => bigwig::Reduce::Min,
        _      => {
            eprintln!("Unknown method: {}. Using 'mean'", method_str);
            bigwig::Reduce::Mean
        }
    };
    let nthreads = nthreads.as_real().unwrap() as usize;
    if nthreads > 1 {
        eprint!("Using {} threads\n", nthreads);
        // rayon::ThreadPoolBuilder::new().num_threads(nthreads).build_global().unwrap();
    } else {
        eprint!("Using single thread\n");
    }
    let iresults = bigwig::read_bigwig(bw_path_str, regions, bin_size, method, Some(nthreads)).unwrap();

    // eprintln!("{:?}", results);
    println!("Done reading bigWig file");
    let results = convert_region_to_robj(&iresults);


    return results;
}



// fn granges_to_regions(gr: &Robj) -> anyhow::Result<Vec<Region>> {
//     // seqnames(gr) -> character()
//     let seq = call!("as.character", call!("seqnames", gr)?)?
//         .as_str_vector()
//         .ok_or_else(|| anyhow!("seqnames(): expected character vector"))?;

//     // start(gr), end(gr) -> integer()
//     let start = call!("start", gr)?
//         .as_integer_vector()
//         .ok_or_else(|| anyhow!("start(): expected integer vector"))?;
//     let end = call!("end", gr)?
//         .as_integer_vector()
//         .ok_or_else(|| anyhow!("end(): expected integer vector"))?;

//     // Optional names(gr)
//     let names_vec = call!("names", gr)?
//         .as_str_vector()
//         .unwrap_or_else(|| Vec::new());

//     if !(seq.len() == start.len() && start.len() == end.len()) {
//         bail!("GRanges vectors have inconsistent lengths");
//     }

//     // Valid chroms filter (same as your BED helper)
//     let valid: std::collections::HashSet<&'static str> = [
//         "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19",
//         "20","21","22","X","Y",
//         "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
//         "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"
//     ].into_iter().collect();

//     let mut out = Vec::with_capacity(seq.len());
//     for i in 0..seq.len() {
//         // R GRanges is 1-based inclusive: convert to 0-based half-open [start0, end0)
//         let s_i = start[i];
//         let e_i = end[i];
//         if s_i.is_na() || e_i.is_na() { continue; }

//         let s0 = (i32::from(s_i) - 1).max(0) as u32;
//         let e0 = (i32::from(e_i)).max(0) as u32;
//         if e0 <= s0 { continue; }

//         let chrom = seq[i].to_string();
//         if !valid.contains(chrom.as_str()) { continue; }

//         let name = if !names_vec.is_empty() && !names_vec[i].is_empty() {
//             Some(names_vec[i].to_string())
//         } else {
//             None
//         };

//         out.push(Region {
//             chrom,
//             start: s0,
//             end: e0,
//             name,
//             coverage_bins: None,
//         });
//     }
//     Ok(out)
// }


pub fn convert_region_to_robj(regions: &Vec<Region>) -> List {
    use extendr_api::prelude::*;

    // Pre-allocate; 1 R list element per region
    let mut elems: Vec<Robj> = Vec::with_capacity(regions.len());

    for r in regions {
        // Build one R list for this region
        let coverage = r.coverage_bins.clone().unwrap_or_default();

        // If you want named fields in each element:
        let region_list = list!(
            chrom    = r.chrom.clone(),
            start    = r.start as i64,   // R likes integers as i32/i64
            end      = r.end   as i64,
            name     = r.name.clone().unwrap_or_default(),
            coverage = coverage          // Vec<f32> converts to R numeric
        );

        elems.push(region_list.into_robj());
    }

    // One final allocation to make an R list from all region lists
    List::from_values(elems)
}




// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod pacbiowdlR;
    fn read_bed;
}
