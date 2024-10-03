use std::{
    io,
    fs,
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
    io::BufRead,
    io::Write,
};
use rust_lapper::{Interval, Lapper};
use flate2::read::MultiGzDecoder;
use flate2::Compression;
use log::error;
use log::info;
use log::warn;
use rustc_hash::FxHashMap;
use gzp::{
    deflate::Gzip,
    ZWriter,
    par::compress::{ParCompress, ParCompressBuilder},
};

pub fn f2m(matches: &clap::ArgMatches) -> Result<(), Box<dyn Error>> {

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

    let output_directory = matches.get_one::<String>("outdir").unwrap();
    info!("Received output directory: {:?}", output_directory);

    let group =matches.get_one::<usize>("group").copied();
    info!("Grouping peaks: {:?}", group);

    let weight = matches.get_one::<usize>("weight").copied();
    info!("Peak weight column: {:?}", weight);
  
    let output_path = Path::new(output_directory);

    let num_threads = *matches.get_one::<usize>("threads").unwrap();

    // Create the directory if it does not exist
    if !output_path.exists() {
        if let Err(e) = fs::create_dir_all(output_path) {
            eprintln!("Failed to create output directory: {}", e);
            std::process::exit(1);
        }
    }

    // make sure output is a directory
    match fs::metadata(output_path) {
        Ok(metadata) => {
            if metadata.is_dir() {
                info!("{:?} is a directory.", output_path);
            } else {
                eprintln!("Provided output is not a directory: {}", output_path.display());
                std::process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("Failed to get metadata for {:?}: {}", output_path, e);
            std::process::exit(1);
        }
    }

    fcount(&frag_file, &bed_file, &cell_file, output_path, group, weight, num_threads)?;
    
    Ok(())
}

fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: &Path,
    output: &Path,
    group: Option<usize>,
    weight: Option<usize>, 
    num_threads: usize,
    
) -> io::Result<()> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );

    // create BED intervals for overlaps with fragment coordinates
    // returns hashmap with each key being chromosome name
    // each value is intervals for that chromosome
    // interval value gives the index of the feature
    // also writes features to output directory to avoid second iteration of file
    // write features
    let feature_path = output.join("features.tsv.gz");
    info!("Writing output feature file: {:?}", &feature_path);
   
    let (total_peaks, mut peaks, peak_vec) = match peak_intervals(bed_file, group, weight, &feature_path, num_threads) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read BED file: {}", e);
            return Err(e);
        }
    };
    // create hashmap for cell barcodes
    let cellreader = File::open(cell_file)
        .map(BufReader::new)?;
    // cell barcode, cell index 
    let mut cells: FxHashMap<String, u32> = FxHashMap::default();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        let index_u32 = index as u32;
        cells.insert(line, index_u32);
    }

    // vector of features
    // each element is hashmap of cell: count (is float)
    let mut peak_cell_counts: Vec<FxHashMap<u32, f32>> = vec![FxHashMap::<u32, f32>::default(); total_peaks];
    // frag file reading
    let frag_file = File::open(frag_file)?;
    let mut reader = BufReader::with_capacity(1024 * 1024, MultiGzDecoder::new(frag_file));

    let mut line_count: u64 = 0;
    let update_interval = 1_000_000;
    let mut line_str = String::new();
    let mut startpos: u32;
    let mut endpos: u32;

    let mut current_chrom = String::new();
    let mut current_lapper: Option<&mut Lapper<u32, usize>> = None;
    let mut cursor: usize = 0;
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

        // # Skip header lines that start with #
        if line.starts_with('#') {
            line_str.clear();
            continue;
        }

        line_count += 1;
        if line_count % update_interval == 0 {
            print!("\rProcessed {} M fragments", line_count / 1_000_000);
            std::io::stdout().flush().expect("Can't flush output");
        }

        // Parse BED entry
        let fields: Vec<&str> = line.split('\t').collect();

        // Check if cell is to be included
        let cell_barcode: &str = fields[3];
        if let Some(&cell_index) = cells.get(cell_barcode) {
            check_end = true;

            // create intervals from fragment entry
            let seqname: &str = fields[0];

            if seqname != current_chrom {
                current_chrom = seqname.to_string();
                current_lapper = peaks.get_mut(&current_chrom);
                cursor = 0;
            }

            // try to parse the coordinates, skip the line if parsing fails
            startpos = match fields[1].trim().parse() {
                Ok(num) => num,
                Err(e) => {
                    warn!("Failed to parse start position: {:?}. Error: {}", line_count, e);
                    line_str.clear();
                    continue;
                }
            };
            
            endpos = match fields[2].trim().parse() {
                Ok(num) => num,
                Err(e) => {
                    warn!("Failed to parse end position: {:?}. Error: {}", line_count, e);
                    line_str.clear();
                    continue;
                }
            };

            if let Some(lapper) = &mut current_lapper {
                // seems to be a problem with seek if lapper has one element
                // set cursor to 0
                if lapper.intervals.len() == 1 {
                    cursor = 0;
                }
                for interval in lapper.seek(startpos, startpos + 1, &mut cursor) {
                    
                    let (peak_index, peak_weight) = peak_vec[interval.val];
                    let peak_end = interval.stop;
                    *peak_cell_counts[peak_index].entry(cell_index).or_insert(0.0) += peak_weight;

                    // Check if fragment end is behind peak end (if so, it overlaps and we don't need a full search)
                    if endpos < peak_end {
                        check_end = false;
                        *peak_cell_counts[peak_index].entry(cell_index).or_insert(0.0) += peak_weight;
                    }
                }
                if check_end {
                    for interval in lapper.seek(endpos, endpos + 1, &mut cursor) {
                        let (peak_index, peak_weight) = peak_vec[interval.val];
                        *peak_cell_counts[peak_index].entry(cell_index).or_insert(0.0) += peak_weight;
                    }
                }
            }
        }
        line_str.clear();
    }
    eprintln!();
    
    // write count matrix
    let counts_path = output.join("matrix.mtx.gz");
    info!("Writing output counts file: {:?}", &counts_path);
    write_matrix_market(&counts_path, &peak_cell_counts, total_peaks, cells.len(), num_threads)
        .expect("Failed to write matrix"); // features stored as rows

    // write cells
    let cell_path = output.join("barcodes.tsv");
    info!("Writing output cells file: {:?}", &cell_path);
    write_cells(&cell_path, cell_file)
        .expect("Failed to write cells");

    Ok(())
}

fn write_cells(
    outfile: &Path,
    cells: &Path,
) -> io::Result<()> {
    // Copy cell barcodes to output directory
    match fs::copy(&cells, &outfile) {
        Ok(bytes_copied) => info!("Successfully copied {} bytes.", bytes_copied),
        Err(e) => eprintln!("Failed to copy file: {}", e),
    }
    Ok(())
}

fn write_matrix_market(
    outfile: &Path,
    peak_cell_counts: &[FxHashMap<u32, f32>],
    nrow: usize,
    ncol: usize,
    num_threads: usize,
) -> io::Result<()> {

    // get nonzero value count
    let nonzero: usize = peak_cell_counts.iter().map(|map| map.len()).sum();

    // create output file
    let writer = File::create(outfile)?;
    let mut encoder: ParCompress<Gzip> = ParCompressBuilder::new()
        .compression_level(Compression::default())  // Set compression level
        .num_threads(num_threads)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))? 
        .from_writer(writer);

    // Create a string buffer to collect all lines
    let mut output = String::new();

    // Write the header for the Matrix Market format
    output.push_str("%%MatrixMarket matrix coordinate integer general\n");
    output.push_str("%%metadata json: {{\"software_version\": \"fragtk-1.1.0\"}}\n");
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
    if !output.is_empty() {
        encoder.write_all(output.as_bytes())?;
    }

    encoder.finish().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    Ok(())
}

fn peak_intervals(
    bed_file: &Path,
    group: Option<usize>,
    weight:Option<usize>, 
    outfile: &Path,
    num_threads: usize,
) -> io::Result<(usize, FxHashMap<String, Lapper<u32, usize>>, Vec<(usize, f32)>)>{ 

    // feature file
    let writer = File::create(outfile)?;
    let mut writer: ParCompress<Gzip> = ParCompressBuilder::new()
        .compression_level(Compression::default())
        .num_threads(num_threads)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
        .from_writer(writer);
    
    // bed file reader
    let file = File::open(bed_file)?;
    let reader = BufReader::new(file);
    
    // hashmap of peak intervals for each chromosome
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();
    
    // Store peak group name and corresponding index
    let mut peak_group_index: FxHashMap<String, usize> = FxHashMap::default();

    // store index and (peak index, weight)
    let mut ind_peak: Vec<(usize, f32)> = Vec::new();
    // let mut ind_peak: FxHashMap<usize, (usize, f32)> = FxHashMap::default();

    // track total number of peaks
    let mut total_peaks: usize = 0;

    // index for peak groups
    let mut current_index: usize = 0;

    // skipped lines
    let mut skipped_lines: usize = 0;

    for (index, line) in reader.lines().enumerate() {

        match line {
            Ok(line) => {
                if line.starts_with('#') {
                    skipped_lines += 1;
                    continue;
                }
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

                    let intervals = chromosome_trees.entry(chromosome.clone()).or_insert_with(Vec::new);
                    // changes: add in weight_value 
                    
                    let mut peak_weight = 1.0; // Default weight 
                    // Update peak weight if weigh_column is specified
                    if let Some(column_number) = weight {
                        if fields.len() > column_number{
                            match fields[column_number].parse() { //change column number here 
                                Ok(num) => peak_weight = num, 
                                Err(_) => {
                                    error!("Line {}: Failed to parse weight value", index + 1);
                                    continue;  
                                }
                            }
                        }
                    }
                    
                    if let Some(group_col) = group{
                        if fields.len() > group_col {
                            let peakgroup: String = match fields[group_col].parse() {
                                Ok(num) => num,
                                Err(_) => {
                                    error!("Line {}: Failed to parse group information", index + 1);
                                    continue;
                                }
                            };

                            let group_index = peak_group_index.entry(peakgroup.clone()).or_insert_with(|| {
                                writeln!(writer, "{}", peakgroup).expect("Failed to write");
                                let idx: usize = current_index - skipped_lines;
                                current_index += 1;
                                idx
                            });

                            intervals.push(Interval {start, stop: end, val: total_peaks});
                            ind_peak.push((*group_index, peak_weight)); 
                                
                    } else { 
                            error!("Line {}: Group Column not found", index +1 );
                        } 
                    }else {
                        intervals.push(Interval { start, stop: end, val: total_peaks});
                        ind_peak.push((index-skipped_lines, peak_weight));
                        writeln!(writer, "{}-{}-{}", chromosome, start, end)?;
                    }
                    total_peaks += 1; // corresponds to if fields.len() >= 3 {
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

    let lapper_map: FxHashMap<String, Lapper<u32, usize>> = chromosome_trees.into_iter()
        .map(|(chr, intervals)| (chr, Lapper::new(intervals)))
        .collect();

    info!("\rTotal No. Peaks {} before change: ", total_peaks);
    assert_eq!(ind_peak.len(), total_peaks, "Length of ind_peak ({}) does not match total_peaks ({})", ind_peak.len(), total_peaks);
    let total_intervals: usize = lapper_map.values().map(|lapper| lapper.len()).sum();
    assert_eq!(total_intervals, total_peaks, "Total number of intervals in lapper_map ({}) does not match total_peaks ({})", total_intervals, total_peaks);

    if group.is_some() { 
        total_peaks = current_index;
    }

    // Finalize the compression, converting GzpError to io::Error
    writer.finish().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    Ok((total_peaks, lapper_map, ind_peak))
}
