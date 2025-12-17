use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Read a text file containing one file path per line
/// Returns Vec<String> of file paths
pub fn read_file_list(path: &str) -> io::Result<Vec<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut files = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            files.push(trimmed.to_string());
        }
    }

    Ok(files)
}

/// Read a FASTA/FNA file and return the sequence as a single String
/// - Skips header lines starting with '>'
/// - Concatenates multiline sequences
/// - Uppercases bases
pub fn read_fna(path: String) -> io::Result<String> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Skip FASTA header
            continue;
        }
        sequence.push_str(line.trim());
    }

    Ok(sequence.to_ascii_uppercase())
}
