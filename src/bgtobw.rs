use clap::Parser;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom, Write, BufWriter};
use std::path::PathBuf;
use std::collections::HashMap;

/// Simple BigWig writer (uncompressed, no index) scaffold
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Input bedGraph (chrom start end value)
    #[clap(value_parser)]
    input: PathBuf,
    /// Chrom sizes file (chrom\tlength)
    #[clap(value_parser)]
    chrom_sizes: PathBuf,
    /// Output bigWig path
    #[clap(value_parser)]
    output: PathBuf,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    // Load chrom sizes into Vec and map name -> id
    let mut contigs = Vec::new();
    let mut sizes = HashMap::new();
    {
        let reader = BufReader::new(File::open(&args.chrom_sizes)?);
        for line in reader.lines() {
            let line = line?;
            let mut parts = line.split_whitespace();
            if let (Some(name), Some(len)) = (parts.next(), parts.next()) {
                let id = contigs.len() as u32;
                contigs.push(name.to_string());
                let length: u64 = len.parse().unwrap_or(0);
                sizes.insert(name.to_string(), (id, length));
            }
        }
    }

    // Open output, reserve header space
    let mut f = BufWriter::new(File::create(&args.output)?);
    // Write empty 64-byte header
    f.write_all(&vec![0u8; 64])?;
    let data_start = f.seek(SeekFrom::Current(0))?;

    // Stream bedGraph, write each record as (u32 chromId, u32 start, u32 end, f32 value)
    let infile = BufReader::new(File::open(&args.input)?);
    for line in infile.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let mut cols = line.split_whitespace();
        if let (Some(chr), Some(s), Some(e), Some(v)) = (cols.next(), cols.next(), cols.next(), cols.next()) {
            if let Some(&(id, _len)) = sizes.get(chr) {
                let start: u32 = s.parse().unwrap_or(0);
                let end:   u32 = e.parse().unwrap_or(0);
                let val:  f32 = v.parse().unwrap_or(0.0);
                f.write_all(&id.to_le_bytes())?;
                f.write_all(&start.to_le_bytes())?;
                f.write_all(&end.to_le_bytes())?;
                f.write_all(&val.to_le_bytes())?;
            }
        }
    }
    // Flush and record end of data
    f.flush()?;
    let data_end = f.seek(SeekFrom::Current(0))?;

    // Now write minimal header with offsets
    let mut h = f.into_inner()?;
    // magic, version, summaryCount
    h.seek(SeekFrom::Start(0))?;
    let sig:    u32 = 0x888FFC26;
    let version:u16 = 4;
    let summary:u16 = 0;
    h.write_all(&sig.to_le_bytes())?;
    h.write_all(&version.to_le_bytes())?;
    h.write_all(&summary.to_le_bytes())?;
    // chromTreeOffset (0), dataOffset, indexOffset == data_end
    let zero64: u64 = 0;
    h.write_all(&zero64.to_le_bytes())?;
    h.write_all(&data_start.to_le_bytes())?;
    h.write_all(&data_end.to_le_bytes())?;
    // fieldCount, definedFieldCount
    let fc: u16 = 4;
    h.write_all(&fc.to_le_bytes())?;
    h.write_all(&fc.to_le_bytes())?;
    // autoSqlOffset, totalSummaryOffset
    h.write_all(&zero64.to_le_bytes())?;
    h.write_all(&zero64.to_le_bytes())?;
    // uncompressBufSize (0), nameIndexOffset (0)
    let zero32: u32 = 0;
    h.write_all(&zero32.to_le_bytes())?;
    h.write_all(&zero64.to_le_bytes())?;
    Ok(())
}
