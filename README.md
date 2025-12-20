# trf_disease_repeats
Analysis of tandem repeats in disease-associated genes using a Rust implementation of a tandem repeat detection algorithm. This project examines baseline repeat quantity and length, to understand how repeat structures influence susceptibility to repeat-expansion disorders.

# Overview
This project implements a two-phase detection algorithm:
* Phase 1: Fast candidate detection using lightweight scoring.
* Phase 2: Refinement using a banded Smith-Waterman local alignment algorithm.
The algorithm considers start and end positions, alignment score, period length, copy number, and base composition in each tandem repeat. The output consists of summary statistics including quantity of repeats found and percent coverage of the original sequence.

# Quickstart
There are two ways to run the algorithm on your machine. You can either download the pre-built release, or clone the repository and build it locally. Note that the second method will require you to have Rust and the necessary crates installed.

## Pre-Built Release
Navigate to the "Releases" section on the right-hand side and download v1.0.0. Keep note of the file path when you download it. Copy the filepath and paste it into the terminal, and then provide an input txt file, max period length, and threshold of your choice. Please note that the algorithm requires the input txt file to have only one sequence .fna file per line. If you are on MacOS, you may need to fix security permissions to run the release. To do so, run "chmod +x trf_disease_repeats" in the terminal, in the same directory as your release download. You will likely also need to approve it from within the System Settings (System Settings -> Privacy & Security -> Security -> Allow trf_disease_repeats from unknown developer).

### Command Format
```bash
./<path-to-release> --input <path-to-txt-input-file-list> --max-period <max-period-length> --threshold <threshold-number>
```

### Example Command
```bash
./trf_disease_repeats --input files.txt --max-period 3 --threshold 20
```

### Sequences
Please see any .fna file for the expected sequence file format. The implementation expects a single line header, and a single line sequence on the next line. To test the release, you can download htt.fna from the Github repository, in the root directory. For a sample files.txt format, you can put the path to your htt.fna download on the first line and save it. To do multiple sequences at once for the same parameters, put the file path on a separate line. See files.txt for an example.

## Cloning & Running locally
Clone the repository:
```bash
git clone https://github.com/kyesteele/trf_disease_repeats.git
cd trf_disease_repeats
```

Ensure Rust is installed:
```bash
rustc --version
cargo --version
```
If Rust is not installed, follow instructions at [https://rust-lang.org/tools/install/](https://rust-lang.org/tools/install/).

Build the project using Cargo:
```bash
cargo build --release
```

The compiled binary will be located at target/release/trf_disease_repeats.

Run the program:
```bash
./target/release/trf_disease_repeats --input <path-to-txt-input-file-list> --max-period <max-period-length> --threshold <threshold-number>
```
