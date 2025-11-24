use crate::alignment;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::cmp;
use std::sync::{Arc, Mutex};

#[derive(Eq, PartialEq, Hash, Clone, Debug)]
pub struct Repeat {
    pub motif: String,
    pub position: usize,
    pub length: usize,
    pub copies: usize,
    pub score: i32,
}

/// Find the smallest repeating unit of a motif
fn smallest_motif(motif: &str) -> &str {
    for len in 1..=motif.len() {
        if motif.len() % len != 0 {
            continue;
        }
        let candidate = &motif[0..len];
        if candidate.repeat(motif.len() / len) == motif {
            return candidate;
        }
    }
    motif
}

/// Filter repeats: remove invalid ones and contained repeats
fn filter_repeats(repeats: Vec<Repeat>) -> Vec<Repeat> {
    let mut repeats_vec = repeats;
    repeats_vec.sort_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then(b.score.cmp(&a.score))
            .then(b.length.cmp(&a.length))
    });

    let mut filtered: Vec<Repeat> = Vec::new();
    for repeat in repeats_vec {
        if repeat.copies >= 2
            && repeat.length >= 2 * repeat.motif.len()
            && repeat.length % repeat.motif.len() == 0
        {
            let minimal_motif = smallest_motif(&repeat.motif).to_owned();
            let length = repeat.length - (repeat.length % minimal_motif.len());
            filtered.push(Repeat {
                motif: minimal_motif,
                position: repeat.position,
                length,
                copies: length / repeat.motif.len(),
                score: repeat.score,
            });
        }
    }

    // Sweep-line: remove contained repeats
    let mut final_repeats: Vec<Repeat> = Vec::new();
    let mut last_end = 0;
    for candidate in filtered {
        let candidate_end = candidate.position + candidate.length;
        if candidate.position >= last_end {
            last_end = candidate_end;
            final_repeats.push(candidate);
        }
    }
    final_repeats
}

/// Scan a sequence for repeats in parallel with progress bar
pub fn scan_sequence(
    sequence: &str,
    max_length: usize,
    threshold: i32,
    pb: &ProgressBar,
) -> Vec<Repeat> {
    let seq_bytes = sequence.as_bytes();
    let chunk_size = 10_000;
    let global_repeats = Arc::new(Mutex::new(Vec::new()));
    let pb = Arc::new(pb.clone());

    let positions: Vec<(usize, usize)> = (0..seq_bytes.len())
        .step_by(chunk_size)
        .map(|start| (start, cmp::min(start + chunk_size, seq_bytes.len())))
        .collect();

    positions.into_par_iter().for_each(|(start, end)| {
        let mut chunk_repeats = Vec::new();

        for i in start..end {
            for len in 1..=max_length {
                if i + len > seq_bytes.len() {
                    break;
                }
                let motif = &seq_bytes[i..i + len];
                let mut j = i;
                let mut copies = 0;
                let mut max_score = 0;

                while j + len <= seq_bytes.len() {
                    let window = &seq_bytes[j..j + len];
                    let score = alignment::align(motif, window);
                    if score >= threshold {
                        copies += 1;
                        max_score = cmp::max(score, max_score);
                        j += len;
                    } else {
                        break;
                    }
                }

                if copies >= 2 {
                    chunk_repeats.push(Repeat {
                        motif: String::from_utf8_lossy(motif).to_string(),
                        position: i,
                        length: copies * len,
                        copies,
                        score: max_score,
                    });
                }
            }
        }

        let filtered_chunk = filter_repeats(chunk_repeats);
        let mut global = global_repeats.lock().unwrap();
        global.extend(filtered_chunk);
        pb.inc((end - start) as u64);
    });

    pb.finish_with_message("Chunks scanned. Performing final sweep-line filtering...");

    // Final sweep-line filter
    let mut repeats = {
        let mut locked = global_repeats.lock().unwrap();
        locked.drain(..).collect::<Vec<_>>()
    };

    repeats.sort_by(|a, b| a.position.cmp(&b.position).then(b.score.cmp(&a.score)));

    let mut final_repeats: Vec<Repeat> = Vec::new();
    let mut last_end = 0;
    for repeat in repeats {
        let candidate_end = repeat.position + repeat.length;
        if repeat.position >= last_end {
            last_end = candidate_end;
            final_repeats.push(repeat);
        }
    }

    final_repeats
}

/// Minimal summary: number of repeats, max repeat length, and repeat density
pub fn print_repeat_summary(gene_name: &str, sequence: &str, max_length: usize, threshold: i32) {
    let pb = ProgressBar::new(sequence.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    let repeats = scan_sequence(sequence, max_length, threshold, &pb);

    let num_repeats = repeats.len();
    let max_repeat_length = repeats.iter().map(|r| r.length).max().unwrap_or(0);
    let repeat_density = if sequence.len() > 0 {
        repeats.iter().map(|r| r.length).sum::<usize>() as f64 / sequence.len() as f64
    } else {
        0.0
    };

    println!("Gene: {}", gene_name);
    println!("Number of repeats: {}", num_repeats);
    println!("Maximum repeat length: {}", max_repeat_length);
    println!("Repeat density: {:.4}", repeat_density);
}
