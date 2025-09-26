#![allow(non_snake_case)]
use extendr_api::prelude::*;
//  use rayon::result;

mod bigwig;
use bigwig::Region;
// use itertools::Itertools;
// use ndarray::Array2;
// use rayon::prelude::*;
// use std::{collections::HashSet, path::Path};
// use anyhow::*;

// Read BED file with group handling
// @export
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
