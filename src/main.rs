use std::{
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
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
                .help("List of cells to include"),
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

    let cell_file = matches.get_one::<String>("cells").map(|cell_path| {
        Path::new(cell_path)
            .canonicalize()
            .expect("Can't find path to input cell file")
    });
    if let Some(ref path) = cell_file {
        info!("Received cell file: {:?}", path);
    }

    fcount(&frag_file, &bed_file, cell_file.as_deref())?;

    Ok(())
}

fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: Option<&Path>,
) -> Result<(), Box<dyn Error>> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );

    let mut tbxreader = tabix::io::indexed_reader::Builder::default().build_from_path(frag_file)?;

    let mut bedreader = File::open(bed_file)
        .map(BufReader::new)
        .map(bed::io::Reader::new)?;
    
    for result in bedreader.records::<3>() {

        let record = result?;
        let region = Region::new(record.reference_sequence_name(), record.start_position()..=record.end_position());

        let query = tbxreader.query(&region)?;
        for entry in query {
            let frag_entry = entry?;
            println!("{}", frag_entry.as_ref());
        }
    }
    Ok(())
}
