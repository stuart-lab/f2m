#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;
use std::{
    io,
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
    io::BufRead,
    io::Write,
};
use clap::{
    Arg,
    Command
};
use rust_lapper::{Interval, Lapper};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use log::error;
use log::info;
use rustc_hash::FxHashMap;


fn main() -> Result<(), Box<dyn Error>> {
    pretty_env_logger::init();

    let matches = Command::new("f2m")
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
        .arg(
            Arg::new("group")
                .long("group")
                .help("Group peaks by variable in fourth BED column")
                .action(clap::ArgAction::SetTrue),
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

    let group = matches.get_flag("group");
    info!("Grouping peaks: {:?}", group);

    fcount(&frag_file, &bed_file, &cell_file, output_file, group)?;

    Ok(())
}

fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: &Path,
    outfile: &String,
    group: bool,
) -> io::Result<()> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );

    // create BED intervals for overlaps with fragment coordinates
    let (total_peaks, peaks) = match peak_intervals(bed_file, group) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read BED file: {}", e);
            return Err(e);
        }
    };

    // create hashmap for cell barcodes
    let cellreader = File::open(cell_file)
        .map(BufReader::new)?;
    
    let mut cells: FxHashMap<String, usize> = FxHashMap::default();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        cells.insert(line.clone(), index);
    }

    // vector of features
    // each element is hashmap of cell: count
    let mut peak_cell_counts: Vec<FxHashMap<usize, u32>> = vec![FxHashMap::<usize, u32>::default(); total_peaks];
    
    // iterate over gzip fragments file
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(frag_file)?));

    // Progress counter
    let mut line_count = 0;
    let update_interval = 1_000_000;

    let mut line_str = String::new();
    let mut startpos: u32;
    let mut endpos: u32;
    let mut check_end: bool;

    loop {
        match reader.read_line(&mut line_str) {
            Ok(0) => break,
            Ok(_) => {},
            Err(e) => {
                error!("Error reading fragment file: {}", e);
                return Err(e);
            }
        }
        let line = &line_str[..line_str.len() - 1];

        // Skip header lines that start with #
        if line.starts_with('#') {
            line_str.clear();
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
            check_end = true;

            // create intervals from fragment entry
            let seqname: &str = fields[0];

            // check if the chromosome exists in peaks
            if peaks.contains_key(seqname) {
                startpos = fields[1].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                endpos = fields[2].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    
                if let Some(olap_start) = find_overlaps(&peaks, seqname, startpos, startpos+1) {

                    for peak_index in olap_start {
                        *peak_cell_counts[peak_index.0].entry(cell_index).or_insert(0) += 1;

                        // check if fragment end is behind peak end (if so, it overlaps and we don't need a full search)
                        if endpos < peak_index.1 {
                            check_end = false;
                            *peak_cell_counts[peak_index.0].entry(cell_index).or_insert(0) += 1;
                        }
                    }
                }
                if check_end {
                    if let Some(olap_end) = find_overlaps(&peaks, seqname, endpos, endpos+1) {
                        for peak_index in olap_end {
                            *peak_cell_counts[peak_index.0].entry(cell_index).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
        line_str.clear();
    }
    
    // write count matrix
    info!("Writing output file: {:?}", &outfile);
    write_matrix_market(&outfile, &peak_cell_counts, total_peaks, cells.len())?; // peaks stored as rows
    Ok(())
}

fn write_matrix_market(
    outfile: &String,
    peak_cell_counts: &[FxHashMap<usize, u32>],
    nrow: usize,
    ncol: usize,
) -> io::Result<()> {

    // get nonzero value count
    let nonzero: usize = peak_cell_counts.iter().map(|map| map.len()).sum();

    // create output file 
    let writer = File::create(outfile)?;
    let mut encoder = GzEncoder::new(writer, Compression::default());

    // Create a string buffer to collect all lines
    let mut output = String::new();

    // Write the header for the Matrix Market format
    output.push_str("%%MatrixMarket matrix coordinate integer general\n");
    output.push_str("%%metadata json: {{\"software_version\": \"f2m-0.1.0\"}}\n");
    output.push_str(&format!("{} {} {}\n", nrow, ncol, nonzero));
    encoder.write_all(output.as_bytes())?;
    output.clear();

    // Collect each peak-cell-count entry into the string buffer
    for (index, hashmap) in peak_cell_counts.iter().enumerate() {
        for (key, value) in hashmap.iter() {
            output.push_str(&format!("{} {} {}\n", index + 1, key + 1, value)); // +1 to convert 0-based to 1-based indices
        }
        // write chunk, clear string
        if index % 5000 == 0 {
            encoder.write_all(output.as_bytes())?;
            output.clear();
        }
    }

    // Write the remaining string buffer
    encoder.write_all(output.as_bytes())?;
    encoder.finish()?;

    Ok(())
}

fn find_overlaps(
    lapper_map: &FxHashMap<String, Lapper<u32, usize>>, 
    chromosome: &str, 
    start: u32, 
    end: u32
) -> Option<Vec<(usize, u32)>> {
    lapper_map.get(chromosome).map(|lapper| {
        lapper.find(start, end).map(|interval| (interval.val, interval.stop)).collect()
    })
}

fn peak_intervals(
    bed_file: &Path,
    group: bool,
) -> io::Result<(usize, FxHashMap<String, Lapper<u32, usize>>)> {
    
    // bed file reader
    let file = File::open(bed_file)?;
    let reader = BufReader::new(file);
    
    // hashmap of peak intervals for each chromosome
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();
    
    // Store peak group name and corresponding index
    let mut peak_group_index: FxHashMap<String, usize> = FxHashMap::default();
    
    // track total number of peaks
    let mut total_peaks: usize = 0;

    // index for peak groups
    let mut current_index: usize = 0;

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

                    if group && (fields.len() >= 4) {
                        let peakgroup: String = match fields[3].parse() {
                            Ok(num) => num,
                            Err(_) => {
                                error!("Line {}: Failed to parse group information", index + 1);
                                continue;
                            }
                        };

                        let group_index = peak_group_index.entry(peakgroup).or_insert_with(|| {
                            let idx: usize = current_index;
                            current_index += 1;
                            idx
                        });

                        intervals.push(Interval { start, stop: end, val: *group_index });
                    } else {
                        intervals.push(Interval { start, stop: end, val: index });
                    }
                    total_peaks += 1;
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

    if group {
        total_peaks = current_index;
    }

    Ok((total_peaks, lapper_map))
}