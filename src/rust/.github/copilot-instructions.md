# PacbioWdlR Rust Extension Instructions

## Architecture Overview

This is a Rust crate that provides R bindings for processing BigWig genomic data files over BED-defined regions. The project uses `extendr-api` for R integration and `bigtools` for BigWig file handling.

### Core Components

- **`src/lib.rs`**: Main R interface with `process_bigwigs_over_bed()` function exported via `#[extendr]`
- **`src/bigwig.rs`**: Core genomic data processing logic with BigWig reading, BED parsing, and matrix computation
- **`Cargo.toml`**: Defines static library (`crate-type = ["staticlib"]`) for R FFI integration

### Data Flow Architecture

1. **Input**: R calls `process_bigwigs_over_bed(bigwig_paths, bed_path, params)`
2. **BED Processing**: `read_bed_with_group()` parses genomic regions with group labels from filename
3. **BigWig Processing**: Parallel processing via Rayon across multiple BigWig files
4. **Matrix Generation**: Creates `Array2<f64>` matrices (regions × windows) for each BigWig
5. **R Integration**: Converts matrices to R data frames via `matrix_to_dataframe()`

## Development Patterns

### Error Handling Strategy
- Use `anyhow::Result` for internal processing, convert to `extendr_api::Result` at R boundaries
- Pattern: `.map_err(|e| anyhow::anyhow!("Context: {}", e))?` for error chaining
- R function signatures must return `Result<List>` with extendr error types

### Genomic Data Conventions
- **Coordinates**: 0-based half-open intervals (start included, end excluded)
- **Chromosome Naming**: Auto-detects UCSC (`chr1`) vs NCBI (`1`) style, converts as needed
- **Processing Styles**: 
  - `PercentOfRegion`: Windows across region ± percentage flanks
  - `Point`: Windows around region center with fixed distances

### Memory & Performance
- Use `rayon::prelude::*` for parallel BigWig processing across files
- `ndarray::Array2` for efficient matrix operations
- Coverage-weighted means for signal aggregation within windows

## Current Build Issues (Fix These First)

1. **Missing Functions**: Implement `process_bigwig_over_regions()` and `matrix_to_dataframe()`
2. **Module Visibility**: Make `bigwig` module items public or re-export in `lib.rs`
3. **BigTools API**: Update `bigtools::bigwig::BigWigRead` to match actual crate API
4. **ExtendrAPI**: Fix `List::new()` usage and parameter parsing from R

## Key Dependencies

- `extendr-api`: R integration (function exports, type conversion)
- `bigtools = "0.5.6"`: BigWig file reading
- `noodles-bed`: BED file parsing  
- `ndarray`: Matrix operations
- `rayon`: Parallel processing

## Testing Workflow

```bash
cargo check                    # Fast syntax/type checking
cargo build                    # Full compilation
cargo test                     # Run unit tests
```

## R Integration Notes

- Parameters from R come as `List` objects, extract with `.elt("param_name")`
- Convert numeric types carefully: `.as_real().unwrap() as usize`
- Matrix output: Convert `Array2` to R-compatible data frame format
- Use `extendr_module!` macro to register exported functions

When modifying this code, preserve the parallel processing architecture and ensure all genomic coordinate handling maintains 0-based conventions.