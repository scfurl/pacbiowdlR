use extendr_api::prelude::*;

mod bigwig;
use itertools::Itertools;
use ndarray::Array2;
use rayon::prelude::*;
use std::{collections::HashSet, path::Path};
use anyhow::*;
// use noodles_bed as bed;


/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod pacbiowdlR;
    fn hello_world;
}


/// Process bigWig files over regions defined in a BED file.
/// @export
/// @param bigwig_paths Character vector of paths to bigWig files.
/// @param bed_path Path to a BED file defining regions.
/// @param params A list of parameters (see `Params` struct in Rust code).
/// @return A list of matrices (as data frames) with processed values for each bigWig file.
#[extendr]
fn process_bigwigs_over_bed(
    bigwig_paths: Vec<String>,
    bed_path: String,
    params: List,
) -> Result<List> {
    // Convert params List to Params struct
    let style_str: String = params.elt("style").unwrap().as_str().unwrap().to_string();
    let style = match style_str.as_str() {
        "percentOfRegion" => Style::PercentOfRegion,
        "point" => Style::Point,
        _ => return Err(anyhow::anyhow!("Invalid style parameter")),
    };     
    let n_windows: usize = params.elt("n_windows").unwrap().as_real().unwrap() as usize;
    let bin_size: usize = params.elt("bin_size").unwrap().as_real().unwrap() as usize;
    let distance_around: Option<f64> = params.elt("distance_around").and_then(|v| v.as_real());
    let distance_up: Option<u32> = params.elt("distance_up").and_then(|v| v.as_real().map(|x| x as u32));
    let distance_down: Option<u32> = params.elt("distance_down").and_then(|v| v.as_real().map(|x| x as u32));
    let params = Params {
        style,
        n_windows,
        bin_size,
        distance_around,
        distance_up,
        distance_down,
    };
    // Read BED file
    let bed_path = std::path::Path::new(&bed_path);
    let regions = read_bed_with_group(bed_path).map_err(|e| anyhow::anyhow!("Failed to read BED file: {}", e))?;
    if regions.is_empty() {
        return Err(anyhow::anyhow!("No valid regions found in BED file"));
    }
    // Process each bigWig file in parallel
    let results: Vec<(String, Array2<f64>)> = bigwig_paths.par_iter()
        .map(|bw_path| {
            let bw = Bw::open(Path::new(bw_path)).map_err(|e| anyhow::anyhow!("Failed to open bigWig file {}: {}", bw_path, e))?;
            let matrix = process_bigwig_over_regions(&bw, &regions, &params)
                .map_err(|e| anyhow::anyhow!("Error processing bigWig file {}: {}", bw_path, e))?;
            Ok((bw_path.clone(), matrix))
        })
        .collect::<Result<Vec<_>>>()?;
    // Convert results to List for R
    let mut output = List::new();
    for (bw_path, matrix) in results {
        let df = matrix_to_dataframe(&matrix);
        output.push(df, &bw_path);
    }
    Ok(output)
}       