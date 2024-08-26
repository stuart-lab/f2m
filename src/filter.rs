use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashSet;
use std::io::Write;
use clap::ArgMatches;

pub fn run(matches: &ArgMatches) -> std::io::Result<()> {
    // Get file paths from command-line arguments
    let cells_file = matches.get_one::<String>("cells").unwrap();
    let fragments_file = matches.get_one::<String>("fragments").unwrap();

    // Load the cell barcodes into a FxHashSet for fast lookups
    let cell_barcodes = load_cells(cells_file)?;

    // Filter the fragment file based on the cell barcodes
    filter_fragments(fragments_file, &cell_barcodes)?;

    Ok(())
}

fn load_cells<P: AsRef<Path>>(path: P) -> std::io::Result<FxHashSet<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut cell_barcodes = FxHashSet::default(); // Change this line

    for line in reader.lines() {
        let line = line?;
        cell_barcodes.insert(line);
    }

    Ok(cell_barcodes)
}

fn filter_fragments<P: AsRef<Path>>(
    fragments_path: P,
    cell_barcodes: &FxHashSet<String>,
) -> std::io::Result<()> {
    // Open the fragment file (gzipped)
    let fragments_file = File::open(fragments_path)?;
    let fragments_reader = BufReader::new(MultiGzDecoder::new(fragments_file));

    // Use stdout for writing uncompressed data
    let stdout = std::io::stdout();
    let mut output_writer = stdout.lock();

    for line in fragments_reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() > 3 && cell_barcodes.contains(fields[3]) {
            writeln!(output_writer, "{}", line)?;
        }
    }

    Ok(())
}