mod alignment;
mod finder;
mod io;

use clap::Parser;
use finder::TrfParams;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    name = "trf_disease_repeats",
    version = "1.0.0",
    about = "Tandem Repeats Finder - Rust implementation (normalized by gene length)"
)]
struct Args {
    #[arg(short, long)]
    input: String,

    #[arg(short = 'l', long, default_value = "500")]
    max_period: usize,

    #[arg(short, long, default_value = "50")]
    threshold: i32,
}

fn main() {
    let start_time = Instant::now();
    let args = Args::parse();

    // read list of files to process
    let file_list = match io::read_file_list(&args.input) {
        Ok(list) => list,
        Err(error) => {
            eprintln!("Error reading file list: {}", error);
            return;
        }
    };

    let mut summaries = Vec::new();

    for file_path in file_list {
        println!("Processing file: {}", file_path);

        // load sequence from fna
        let sequence = match io::read_fna(file_path.clone()) {
            Ok(seq) => seq,
            Err(error) => {
                eprintln!("Error reading {}: {}", file_path, error);
                continue;
            }
        };

        println!("Loaded sequence ({} bases)", sequence.len());
        println!("Running TRF-like algorithm (normalized by gene length)...\n");

        // TRF parameters
        let params = TrfParams {
            match_weight: 2,
            mismatch_penalty: 7,
            indel_penalty: 7,
            min_score: args.threshold,
            max_period: args.max_period,
            ..Default::default()
        };

        // run TRF and collect summary of results
        let summary = finder::collect_trf_summary(&sequence, params);
        summaries.push((file_path, summary));
    }

    // output results
    println!("\n====== FINAL SUMMARY ======");
    for (file_name, summary) in &summaries {
        println!("File: {}", file_name);
        println!("# Repeats: {}", summary.num_repeats);
        println!("% Repeats of sequence: {:.4}%", summary.percent_repeats);
        println!("------------------------");
    }

    let duration = start_time.elapsed();
    println!(
        "Total Execution Time: {:.2} seconds",
        duration.as_secs_f64()
    );
}
