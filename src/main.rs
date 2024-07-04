use std::{
    path::Path,
    error::Error,
    fs::File,
    io::Write,
    io::BufReader,
    io::BufRead,
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

    fcount(&frag_file, &bed_file, &cell_file)?;

    Ok(())
}

// TODO add stdout for mtx row,column,value entries
fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: &Path,
) -> Result<(), Box<dyn Error>> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );

    // printing to stdout. TODO: change to file output
    let mut stdout = std::io::stdout().lock();
    writeln!(stdout, "%%MatrixMarket matrix coordinate integer general").unwrap();

    // tabix file reader
    let mut tbxreader = tabix::io::indexed_reader::Builder::default().build_from_path(frag_file)?;

    // bed file reader
    let mut bedreader = File::open(bed_file)
        .map(BufReader::new)
        .map(bed::io::Reader::new)?;

    //  create hashmap for cell barcodes
    let cellreader = File::open(cell_file)
        .map(BufReader::new)?;
    let mut cells = HashMap::new();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        cells.insert(line.clone(), index);
        // writeln!(stdout, "{}", line).unwrap();
    }

    // create region index variable for matrix row index
    let mut bed_idx = 1;

    // TODO
    // might be better to output to a directory rather than stdout as we need the features and cells txt files as well, and to add the matrix dimension to the header
    
    for result in bedreader.records::<3>() {

        let record = result?;
        let region = Region::new(record.reference_sequence_name(), record.start_position()..=record.end_position());

        let query = tbxreader.query(&region)?;
        for entry in query {
            let frag_entry = entry?;

            // TODO match cell barcodes in hashmap, create vector

            writeln!(stdout, "{}", frag_entry.as_ref()).unwrap();  // placeholder
        }

        // - match cell barcodes in entry to cell barcodes in list
        // - convert cell barcode to list index
        // - count fragments per cell barcode
        // - append row,column,value to mtx file
        // - increment region index
        bed_idx = bed_idx + 1;

        // need to track number of nonzero entries (total lines written) and number of bed entries, add this to mtx header at the end?

    }
    Ok(())
}
