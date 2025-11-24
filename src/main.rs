mod alignment;
mod finder;
mod io;
use clap::Parser;

fn main() {
    let args = Args::parse();
    let sequence = match io::read_fna(args.input) {
        Ok(sequence) => sequence,
        Err(error) => {
            eprintln!("Error: {}", error);
            return;
        }
    };
    println!("Loaded sequence ({} bases)", sequence.len());
    println!("Running TRF algorithm...");

    // Use minimal summary output
    finder::print_repeat_summary("GeneName", &sequence, args.length, args.threshold);
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    input: String,
    #[arg(short, long)]
    length: usize,
    #[arg(short, long)]
    threshold: i32,
}
