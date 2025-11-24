use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

// Read a .fna file with a single line header
pub fn read_fna(filename: String) -> Result<String, std::io::Error> {
    let f = File::open(filename)?;
    let mut reader = BufReader::new(f);
    let mut buffer: String = String::new();
    let mut sequence: String = String::new();
    let mut header: String = String::new();
    reader.read_line(&mut header)?;
    loop {
        let bytes_read = reader.read_line(&mut buffer)?;
        if bytes_read <= 0 {
            break;
        } else {
            sequence += &buffer;
        }
    }

    Ok(sequence)
}
