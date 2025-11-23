pub mod alignment;
pub mod finder;
pub mod io;
use clap::Parser;

fn main() {
    let args = Args::parse();
    println!("Hello, world!");
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

    #[arg(short, long)]
    output: String,
}
