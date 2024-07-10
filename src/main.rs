use std::{
    io,
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
    io::BufRead,
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
    let outf = File::create(outfile)?;
    
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

    // nonzero entry counter
    let mut nz = 0;

    // iterate over gzip fragments file
    let reader = BufReader::new(MultiGzDecoder::new(File::open(frag_file)?));

    for line in reader.lines() {
        let line = line?;
        
        // parse bed entry
        let fields: Vec<&str> = line.split('\t').collect();

        // check if cell is to be included
        let cb: &str = fields[3];
        if cells.contains_key(cb) {

            // create intervals from fragment entry
            let seqname: &str = fields[0];
            let startpos: u32 = fields[1].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let endpos: u32 = fields[2].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            // find overlapping peaks
            if let Some(olap_start) = find_overlaps(&peaks, seqname, startpos, startpos+1) {
                for interval in olap_start {
                    // record cell barcode and add count for start position fragment insertion

                    nz += 1;
                    println!("{}", interval);
                }
            };

            if let Some(olap_end) = find_overlaps(&peaks, seqname, endpos, endpos+1) {
                for interval in olap_end {
                    // record cell barcode and add count for end position fragment insertion
                    // println!("{}", interval);
                }
            };
        }

        nz += 1; // todo: update this 

        // TODO set endpos so we don't need to search all peaks again for next entry
        // might speed things up

    }
    println!("{}", nz);

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

// // steps:
// 1. read line
// 2. overlap with interal tree
// 3. if match, add to hashmap with peak, cell, count
// 4. read next line, if start position beyond end of last peak check for overlaps with peak set again
// 5. if start position of next line is less than end position of the previous peak, add to quantification for same peak
// 6. continue until no more lines
// 7. write hashmap containing peak, cell, count values to file as matrix market format

// time to iterate over all lines of 2.2 Gb file: 1 min (227M fragments)
//                                   200 mB file: 6 s   (20M fragments)
