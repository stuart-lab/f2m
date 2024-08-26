use std::io;
use std::error::Error;
use std::path::Path;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;
use std::thread;
use std::sync::mpsc;
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use log::info;

pub fn cellselect(matches: &clap::ArgMatches) -> Result<(), Box<dyn Error>> {

    let frag_file = Path::new(matches.get_one::<String>("fragments").unwrap())
        .canonicalize()
        .expect("Can't find path to input fragment file");
    info!("Received fragment file: {:?}", frag_file);

    let output_file = matches.get_one::<String>("outfile").unwrap();
    info!("Output file: {:?}", output_file);

    let threshold: usize = matches
        .get_one::<String>("threshold")
        .expect("defaulted")
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Failed to parse threshold as usize");
            std::process::exit(1);
        });
    info!("Cell count cutoff: {:?}", threshold);

    let bc_count = count_barcodes(&frag_file)?;
    let selected = select_barcodes(&bc_count, &threshold)?;

    // Output results to the specified file
    let mut writer = File::create(output_file)?;
    let mut output = String::new();

    for (barcode, count) in &bc_count {
        output.push_str(&format!("{}\t{}\n", barcode, count));
    }
    writer.write_all(output.as_bytes())?;

    // print selected cells to stdout
    for cell in selected {
        println!("{}", cell);
    }

    Ok(())
}

fn select_barcodes(
    barcodes: &FxHashMap<String, usize>,
    count_cutoff: &usize,
) -> io::Result<Vec<String>> {

    // iterate over key, value entries, adding cells if count is greater than threshold
    let mut filtered_cells = Vec::new();
    for (cell_barcode, &count) in barcodes.iter() {
        if count > *count_cutoff {
            filtered_cells.push(cell_barcode.clone());
        }
    }

    Ok(filtered_cells)
}

fn count_barcodes(frag_file: &Path,) -> io::Result<FxHashMap<String, usize>> {

    // hashmap for cell barcode counts
    let mut cells: FxHashMap<String, usize> = FxHashMap::default();

    // Create a channel for communication between the decompression and processing threads
    let (tx, rx) = mpsc::sync_channel(500);

    // Spawn the decompression thread
    let frag_file = frag_file.to_path_buf();
    let decompress_handle = thread::spawn(move || {
        let reader = BufReader::new(MultiGzDecoder::new(File::open(frag_file).expect("Failed to open fragment file")));
        for line in reader.lines() {
            let line = line.expect("Failed to read line");
            if tx.send(line).is_err() {
                break;
            }
        }
    });

    // Progress counter
    let mut line_count: u64 = 0;
    let update_interval = 1_000_000;

    for line in rx {

        // Skip header lines that start with #
        if line.starts_with('#') {
            continue;
        }

        line_count += 1;
        if line_count % update_interval == 0 {
            eprint!("\rProcessed {} M fragments", line_count / 1_000_000 );
            std::io::stdout().flush().expect("Can't flush output");
        }

        // parse bed entry
        let fields: Vec<&str> = line.split('\t').collect();

        // update count for cell barcode
        if let Some(cell_barcode) = fields.get(3) {
            let cell_barcode = cell_barcode.to_string();
            *cells.entry(cell_barcode).or_insert(0) += 1;
        }
    }
    eprintln!();

    // Join thread to ensure it completes
    decompress_handle.join().expect("Failed to join decompression thread");

    Ok(cells)
}