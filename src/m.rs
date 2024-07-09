use std::{
    path::Path,
    error::Error,
    fs::File,
    io::Write,
    io::BufReader,
    io::BufRead,
    io::copy,
    io::Seek,
    io::SeekFrom,
    collections::HashMap,
};
use clap::{
    Arg,
    Command
};
use noodles::{
    tabix,
    bed,
    core::Region,
    };
use log::info;
use pretty_env_logger;
use tempfile::NamedTempFile;

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
) -> Result<(), Box<dyn Error>> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );

    // create temp file
    let mut temp_file = NamedTempFile::new()?;

    // create output file 
    let mut outf = File::create(outfile)?;

    // tabix file reader
    let mut tbxreader = tabix::io::indexed_reader::Builder::default().build_from_path(frag_file)?;
    // let index = tabix::read("./test/atac_pbmc_1k_nextgem_fragments.tsv.gz.tbi")?;
    // println!("{:?}", index);
    
    // bed file reader
    let mut bedreader = File::open(bed_file)
        .map(BufReader::new)
        .map(bed::io::Reader::new)?;

    // // read bed coordinates into vector
    // let mut bed_coords: Vec<Region> = Vec::new();
    // for result in bedreader.records::<3>() {   
    //     let record = result?;
    //     let region = Region::new(record.reference_sequence_name(), record.start_position()..=record.end_position());
    //     bed_coords.push(region);
    // }
    // println!("{:?}", bed_coords);

    //  create hashmap for cell barcodes
    let cellreader = File::open(cell_file)
        .map(BufReader::new)?;
    
    let mut cells: HashMap<String, usize> = HashMap::new();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        cells.insert(line.clone(), index);
    }

    // create region index variable for matrix row index
    let mut bed_idx = 1;

    // nonzero entry counter
    let mut nz = 0;

    let mut cb_vec = Vec::new();
    let mut bed_vec = Vec::new();
    let mut value_vec = Vec::new();

    for result in bedreader.records::<3>() {
        
        let record = result?;
        let region = Region::new(record.reference_sequence_name(), record.start_position()..=record.end_position());
        let query = tbxreader.query(&region)?;
        let mut cb = Vec::new();

        for entry in query {
            let frag_entry = entry?;
            let lines: Vec<&str> = frag_entry.as_ref().split('\t').collect();
            cb.push(lines[3].to_string());
        }
        
        // create frequency table and cell index hashmap
        let mut cb_freq = HashMap::new();
        for cell in cb {
            if let Some(&idx) = cells.get(&cell) {
                let counter = cb_freq.entry(idx).or_insert(0);
                *counter += 1;
            }
        }

        for (index, count) in cb_freq {
            nz = nz + 1;
            cb_vec.push(index+1);
            value_vec.push(count);
            bed_vec.push(bed_idx);
        }

        if cb_vec.len() > 1000000 {
            let mut combined_string = String::new();
            for i in 0..cb_vec.len() {
                combined_string.push_str(&format!("{} {} {}\n", bed_vec[i], cb_vec[i], value_vec[i]));
            }

            // write to file
            temp_file.write_all(combined_string.as_bytes())?;

            // clear vectors
            cb_vec.clear();
            bed_vec.clear();
            value_vec.clear();
        }

        // increment region index
        bed_idx = bed_idx + 1;
    }

    // write remaining values to file
    let mut combined_string = String::new();
    for i in 0..cb_vec.len() {
        combined_string.push_str(&format!("{} {} {}\n", bed_vec[i], cb_vec[i], value_vec[i]));
    }
    // write to file
    temp_file.write_all(combined_string.as_bytes())?;

    let _ = outf.write_all("%%MatrixMarket matrix coordinate integer general\n".as_bytes());
    let _ = outf.write_all("%%metadata json: {{\"software_version\": \"f2m-0.1.0\"}}\n".as_bytes());
    let formatted_string = format!("{} {} {}\n", bed_idx - 1, cells.len(), nz);
    outf.write_all(formatted_string.as_bytes())?;

    // cat header and output file
    temp_file.as_file_mut().seek(SeekFrom::Start(0))?;
    copy(&mut temp_file, &mut outf)?;

    Ok(())
}


// works but quite slow, 20 min for one file
// most of the time is taken by reading from dist with position-based query
// streaming through the whole file and finding fragments that overlap peak coordinates would be faster
// can stream through 2 Gb file in about 10 s
// reading all regions takes >20 min

// TODO
// change to stream through fragments in the file
// for each fragment, check overlap with peaks
// assume peaks are sorted so once we find a fragment with overlap we can remove all the preceeding peaks from further checks
// if the fragment overlaps a peak, record peak index, cell index
// how to get the total count of fragments for each peak...

// if streaming over gzip file we don't need the noodles library
