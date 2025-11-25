mod alignment;
mod finder;
mod io;

use clap::Parser;
use finder::TrfParams;
use std::time::Instant;

fn main() {
    let start_time = Instant::now();
    let args = Args::parse();

    // read in sequence
    let sequence = match io::read_fna(args.input.clone()) {
        Ok(sequence) => sequence,
        Err(error) => {
            eprintln!("Error reading file: {}", error);
            return;
        }
    };

    println!("Loaded sequence ({} bases)", sequence.len());
    println!("Running TRF algorithm...\n");

    let params = TrfParams {
        match_weight: 2,
        mismatch_penalty: 7,
        indel_penalty: 7,
        min_score: args.threshold,
        max_period: args.max_period,
    };

    // run trf and print results
    finder::print_trf_output(&sequence, params);
    let duration = start_time.elapsed();
    let secs = duration.as_secs();
    println!("Execution Time: {} seconds.", secs);
}

#[derive(Parser, Debug)]
#[command(
    name = "trf-rs",
    version = "1.0.0",
    about = "Tandem Repeats Finder - Rust implementation",
    long_about = "Find tandem repeats in DNA sequences using the TRF algorithm"
)]
struct Args {
    // input gene file
    #[arg(short, long)]
    input: String,

    // maximum length of repeats
    #[arg(short = 'l', long, default_value = "500")]
    max_period: usize,

    // threshold (for alignment scoring)
    #[arg(short, long, default_value = "50")]
    threshold: i32,
}
