#[cfg(not(target_os = "windows"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: std::alloc::System = std::alloc::System;

use clap::{Command, Arg, ArgAction};
use std::error::Error;

mod f2m;
mod cellselect;


fn main() -> Result<(), Box<dyn Error>> {

    let matches = Command::new("fragtk")
        .version("1.0")
        .author("Tim Stuart")
        .about("Fragment file processing tools")
        .arg_required_else_help(true)
        .subcommand(
            Command::new("matrix")
                .about("Create a feature x cell matrix from a fragment file")
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
                        .help("BED file containing non-overlapping genomic regions to quantify")
                        .required(true),
                )
                .arg(
                    Arg::new("cells")
                        .short('c')
                        .long("cells")
                        .help("File containing cell barcodes to include")
                        .required(true),
                )
                .arg(
                    Arg::new("outdir")
                        .short('o')
                        .long("outdir")
                        .help("Output directory name. Directory will be created if it does not exist.
        The output directory will contain matrix.mtx.gz, features.tsv, barcodes.tsv")
                        .required(true),
                )
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .help("Number of compression threads to use")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("4")
                        .required(false),
                )
                .arg(
                    Arg::new("group")
                        .long("group")
                        .help("Group peaks by variable in fourth BED column")
                        .action(ArgAction::SetTrue),
                )
        )
        .subcommand(
            Command::new("count")
            .about("Count number of fragments per cell barcode")
            .arg(
                Arg::new("fragments")
                    .short('f')
                    .long("fragments")
                    .value_name("FILE")
                    .help("Path to the fragment file")
                    .required(true),
            )
            .arg(
                Arg::new("outfile")
                    .short('o')
                    .long("outfile")
                    .value_name("FILE")
                    .help("Name of output file. Will contain each cell barcodes and its total count")
                    .required(true),
            )
            .arg(
                Arg::new("threshold")
                    .short('t')
                    .long("threshold")
                    .value_name("NUMBER")
                    .help("Sets the threshold value for filtering cell barcodes")
                    .default_value("200"),
            )
        )
        .get_matches();

    pretty_env_logger::init_timed();

    match matches.subcommand() {
        Some(("matrix", sub_matches)) => f2m::f2m(sub_matches)?,
        Some(("count", sub_matches)) => cellselect::cellselect(sub_matches)?,
        _ => {

        }
    }

    Ok(())
}