#![allow(non_snake_case)]
use extendr_api::prelude::*;
//  use rayon::result;

mod bigwig;
use bigwig::Region;
use std::result::Result::Ok;


// Extract data from Bigwig using regions from a BED file
/// @export
/// @keywords internal
#[extendr]
pub fn read_bigwig_using_bed(bed_path: Robj, bw_path: Robj, bin_size: Robj, method: Robj, nthreads: Robj) -> List {
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
    let nthreads = nthreads.as_integer().unwrap() as usize;
    if nthreads > 1 {
        eprint!("Using {} threads\n", nthreads);
        // rayon::ThreadPoolBuilder::new().num_threads(nthreads).build_global().unwrap();
    } else {
        eprint!("Using single thread\n");
    }
    let iresults = bigwig::read_bigwig(bw_path_str, regions, bin_size, method, Some(nthreads), bigwig::Missing::NaN).unwrap();

    // eprintln!("{:?}", results);
    println!("Done reading bigWig file");
    let results = convert_region_to_robj(&iresults);


    return results;
}


// Extract data from Bigwig using granges from a GRanges object
/// @export
/// @keywords internal
#[extendr]
pub fn read_bigwig_using_gr(features: Robj, bw_path: Robj, bin_size: Robj, method: Robj, nthreads: Robj) -> List {
    // Check inputs
    let bw_path_str = bw_path.as_str().unwrap();
    let bw_path = std::path::Path::new(bw_path_str);
    // check if file exists
    if !bw_path.exists() {
        eprintln!("File {} does not exist", bw_path.display());
        return list!(0);
    }
    eprint!("Reading bigWig file from {}...\n", bw_path.display());

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
    let nthreads = nthreads.as_integer().unwrap() as usize;
    if nthreads > 1 {
        eprint!("Using {} threads\n", nthreads);
        // rayon::ThreadPoolBuilder::new().num_threads(nthreads).build_global().unwrap();
    } else {
        eprint!("Using single thread\n");
    }

    // ---- Parse features list ----
    let flist = match features.as_list() {
        Some(l) => l,
        None => {
            eprintln!("`features` must be a list: list(seqnames, start, start_plus_width, strand)");
            return list!(0);
        }
    };
    if flist.len() < 4 {
        eprintln!("`features` must have 4 elements: seqnames, start, start_plus_width, strand");
        return list!(0);
    }

    let mut it = flist.iter();
    let binding = it.next().unwrap();
    let seqs    = binding.1.as_str_vector().unwrap_or_default();
    let starts  = it.next().unwrap().1.as_integer_vector().unwrap_or_default(); // Vec<Rint>
    let ends_p1 = it.next().unwrap().1.as_integer_vector().unwrap_or_default(); // start+width
    let binding = it.next().unwrap();
    let strand  = binding.1.as_str_vector().unwrap_or_default();
    let binding = it.next().unwrap();
    let names = binding.1.as_str_vector().unwrap_or_default();

    // eprintln!("Seqs len: {}, starts len: {}, ends_p1 len: {}, strand len: {}", seqs.len(), starts.len(), ends_p1.len(), strand.len());

    if !(seqs.len() == starts.len() && starts.len() == ends_p1.len()) {
        eprintln!("lengths of seqnames, start, and start_plus_width differ");
        return list!(0);
    }

    // ---- Build Regions (0-based, half-open) ----
    let mut regions: Vec<bigwig::Region> = Vec::with_capacity(seqs.len());
    for i in 0..seqs.len() {
        let s_ri = starts[i];
        let e1_ri = ends_p1[i];
        if s_ri.is_na() || e1_ri.is_na() { continue; }

        let s1 = i32::from(s_ri);
        let e1 = i32::from(e1_ri);  

        if s1 <= 0 || e1 <= 0 { continue; }

        let start0 = (s1 - 1).max(0) as u32;
        let end0   = e1.max(s1) as u32;

        if end0 <= start0 { continue; }

        let chrom = seqs[i].to_string();
        regions.push(bigwig::Region {
            chrom,
            start: start0,
            end: end0,
            strand: strand.get(i).and_then(|s| {
                if s.len() == 1 {
                    let c = s.chars().next().unwrap();
                    if c == '+' || c == '-' || c == '.' {
                        Some(c)
                    } else {
                        None
                    }
                } else {
                    None
                }
            }),
            name: Some(names[i].to_string()),
            coverage_bins: None,
        });
    }

    
    let iresults = bigwig::read_bigwig(bw_path_str, regions, bin_size, method, Some(nthreads), bigwig::Missing::NaN).unwrap();

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

        let strand = r.strand.unwrap_or('.'); // Use '.' if strand is None
        let strand_char = strand.to_string();

        // If you want named fields in each element:
        let region_list = list!(
            chrom    = r.chrom.clone(),
            start    = r.start as i64,   // R likes integers as i32/i64
            end      = r.end   as i64,
            name     = r.name.clone().unwrap_or_default(),
            strand   = strand_char, // Convert char to String for R
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
    fn read_bigwig_using_bed;
    fn read_bigwig_using_gr;
}