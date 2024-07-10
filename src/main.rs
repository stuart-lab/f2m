use std::{
    io,
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
    io::BufRead,
    io::Write,
    collections::HashMap,
};
use clap::{
    Arg,
    Command
};
use rust_lapper::{Interval, Lapper};
use flate2::read::MultiGzDecoder;
use log::error;
use log::info;
use pretty_env_logger;

type Iv = Interval<u32, usize>;

fn main() -> Result<(), Box<dyn Error>> {
    pretty_env_logger::init();

    let matches = Command::new("fcount")
        .version("0.1")
        .author("Tim Stuart")
        .about("Create a cell x feature matrix")
        .arg(
            Arg::new("fragments")
                .short('f')
                .long("fragments")
                .help("Path to the fragment file")
                .required(true),
        )
        .arg(
            Arg::new("bed")
                .short('b')
                .long("bed")
                .help("BED file containing regions to quantify")
                .required(true),
        )
        .arg(
            Arg::new("cells")
                .short('c')
                .long("cells")
                .help("List of cells to include")
                .required(true),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .help("Output file name")
                .required(true),
        )
        .get_matches();

    let frag_file = Path::new(matches.get_one::<String>("fragments").unwrap())
        .canonicalize()
        .expect("Can't find path to input fragment file");
    info!("Received fragment file: {:?}", frag_file);

    let bed_file = Path::new(matches.get_one::<String>("bed").unwrap())
        .canonicalize()
        .expect("Can't find path to input BED file");
    info!("Received BED file: {:?}", bed_file);

    let cell_file = Path::new(matches.get_one::<String>("cells").unwrap())
        .canonicalize()
        .expect("Can't find path to input cell file");
    info!("Received cell file: {:?}", cell_file);

    let output_file = matches.get_one::<String>("outfile").unwrap();
    info!("Received output file path: {:?}", output_file);

    fcount(&frag_file, &bed_file, &cell_file, &output_file)?;

    Ok(())
}

fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: &Path,
    outfile: &String,
) -> io::Result<()> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );

    // create output file 
    let mut outf = File::create(outfile)?;
    
    // create BED intervals for overlaps with fragment coordinates
    let peaks = match peak_intervals(bed_file) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read BED file: {}", e);
            return Err(e);
        }
    };

    // create hashmap for cell barcodes
    let cellreader = File::open(cell_file)
        .map(BufReader::new)?;
    
    let mut cells: HashMap<String, usize> = HashMap::new();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        cells.insert(line.clone(), index);
    }

    // HashMap to store counts for peak-cell combinations
    let mut peak_cell_counts: HashMap<(usize, usize), usize> = HashMap::new();

    // iterate over gzip fragments file
    let reader = BufReader::new(MultiGzDecoder::new(File::open(frag_file)?));

    // Progress counter
    let mut line_count = 0;
    let update_interval = 1_000_000;

    for line in reader.lines() {
        let line = line?;

        // Skip header lines that start with #
        if line.starts_with('#') {
            continue;
        }

        line_count += 1;
        if line_count % update_interval == 0 {
            print!("\rProcessed {} M fragments", line_count / 1_000_000 );
            std::io::stdout().flush().expect("Can't flush output");
        }
        
        // parse bed entry
        let fields: Vec<&str> = line.split('\t').collect();

        // check if cell is to be included
        let cell_barcode: &str = fields[3];
        if let Some(&cell_index) = cells.get(cell_barcode) {

            // create intervals from fragment entry
            let seqname: &str = fields[0];

            // check if the chromosome exists in peaks
            if peaks.contains_key(seqname) {
                let startpos: u32 = fields[1].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let endpos: u32 = fields[2].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    
                // find overlapping peaks
                if let Some(olap_start) = find_overlaps(&peaks, seqname, startpos, startpos+1) {
                    for peak_index in olap_start {
                        *peak_cell_counts.entry((peak_index, cell_index)).or_insert(0) += 1;
                    }
                }
                if let Some(olap_end) = find_overlaps(&peaks, seqname, endpos, endpos+1) {
                    for peak_index in olap_end {
                        *peak_cell_counts.entry((peak_index, cell_index)).or_insert(0) += 1;
                    }
                }
            }
        }

        // TODO set endpos so we don't need to search all peaks again for next entry
        // might speed things up

    }

    // write count matrix
    let num_peaks = peaks.values().map(|lapper| lapper.intervals.len()).sum::<usize>();
    info!("Writing output file: {:?}", outf);
    write_matrix_market(&mut outf, &peak_cell_counts, num_peaks, cells.len())?; // peaks stored as rows

    Ok(())
}

fn write_matrix_market<W: Write>(
    writer: &mut W,
    peak_cell_counts: &HashMap<(usize, usize), usize>,
    nrow: usize,
    ncol: usize,
) -> io::Result<()> {
    // Create a string buffer to collect all lines
    let mut output = String::new();

    // Write the header for the Matrix Market format
    output.push_str("%%MatrixMarket matrix coordinate integer general\n");
    output.push_str("%%metadata json: {{\"software_version\": \"f2m-0.1.0\"}}\n");
    output.push_str(&format!("{} {} {}\n", nrow, ncol, peak_cell_counts.len()));

    // Collect each peak-cell-count entry into the string buffer
    for ((peak_index, cell_index), count) in peak_cell_counts {
        output.push_str(&format!("{} {} {}\n", peak_index + 1, cell_index + 1, count)); // +1 to convert 0-based to 1-based indices
    }

    // Write the entire string buffer to the writer in one go
    writer.write_all(output.as_bytes())?;

    Ok(())
}

fn find_overlaps(
    lapper_map: &HashMap<String, Lapper<u32, usize>>, 
    chromosome: &str, 
    start: u32, 
    end: u32
) -> Option<Vec<usize>> {
    lapper_map.get(chromosome).map(|lapper| {
        lapper.find(start, end).map(|interval| interval.val).collect()
    })
}


fn peak_intervals(
    bed_file: &Path,
) -> io::Result<HashMap<String, Lapper<u32, usize>>> {
    
    // bed file reader
    let file = File::open(bed_file)?;
    let reader = BufReader::new(file);

    // chromosome trees
    let mut chromosome_trees: HashMap<String, Vec<Iv>> = HashMap::new();

    for (index, line) in reader.lines().enumerate() {
        match line {
            Ok(line) => {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() >= 3 {
                    let chromosome = fields[0].to_string();
                    let start: u32 = match fields[1].parse() {
                        Ok(num) => num,
                        Err(_) => {
                            error!("Line {}: Failed to parse start position", index + 1);
                            continue;
                        }
                    };
                    let end: u32 = match fields[2].parse() {
                        Ok(num) => num,
                        Err(_) => {
                            error!("Line {}: Failed to parse end position", index + 1);
                            continue;
                        }
                    };

                    let intervals = chromosome_trees.entry(chromosome).or_insert_with(Vec::new);
                    intervals.push(Interval { start, stop: end, val: index });
                } else {
                    error!("Line {}: Less than three fields", index + 1);
                }
            },
            Err(e) => {
                error!("Error reading line {}: {}", index + 1, e);
                break;
            }
        }
    }

    let lapper_map = chromosome_trees.into_iter()
        .map(|(chr, intervals)| (chr, Lapper::new(intervals)))
        .collect();

    Ok(lapper_map)
}