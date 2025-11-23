use crate::alignment;
use std::cmp;
use std::collections::HashSet;

#[derive(Eq, PartialEq, Hash)]
struct Repeat {
    motif: String,
    position: usize,
    length: usize,
    copies: usize,
    score: i32,
}

// only keep true, non-overlapping, non-fully-contained, long enough candidate repeats
fn filter_repeats(repeats: Vec<Repeat>) -> Vec<Repeat> {
    // need an ordered collection, but hashset lets us remove duplicates first
    let repeats_hashset: HashSet<Repeat> = HashSet::from_iter(repeats);
    let mut filtered: Vec<Repeat> = Vec::new();
    let mut repeats_vec: Vec<Repeat> = repeats_hashset.into_iter().collect();
    // sort repeats by position, then score, then length
    repeats_vec.sort_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then(b.score.cmp(&a.score))
            .then(b.length.cmp(&a.length))
    });
    // remove trivial repeats
    for repeat in repeats_vec {
        if repeat.copies >= 2 && repeat.length >= 2 * repeat.motif.len() {
            filtered.push(repeat);
        }
    }
    let mut final_repeats: Vec<Repeat> = Vec::new();
    // filter out fully contained repeats
    for candidate in filtered {
        let candidate_start = candidate.position;
        let candidate_end = candidate.position + candidate.length;
        let mut is_contained = false;
        for kept in &final_repeats {
            let kept_start = kept.position;
            let kept_end = kept.position + kept.length;
            if candidate_start >= kept_start && candidate_end <= kept_end {
                if kept.score >= candidate.score {
                    is_contained = true;
                    break;
                }
            }
        }
        if !is_contained {
            final_repeats.push(candidate);
        }
    }
    final_repeats
}

// scan a sequence for repeats
fn scan_sequence(sequence: &str, max_length: usize, threshold: i32) -> Vec<Repeat> {
    let motifs: HashSet<String> = generate_motifs(max_length, sequence);
    let mut repeats: Vec<Repeat> = Vec::new();
    for i in 0..sequence.len() {
        for motif in &motifs {
            if i + motif.len() > sequence.len() {
                continue;
            } else {
                let window = &sequence[i..i + motif.len()];
                let score = alignment::align(motif, window);
                if score >= threshold {
                    let repeat = extend(motif, sequence, i, threshold);
                    repeats.push(repeat);
                }
            }
        }
    }
    let filtered_repeats = filter_repeats(repeats);
    filtered_repeats
}

// attempt to extend a repeat
fn extend(motif: &str, sequence: &str, position: usize, threshold: i32) -> Repeat {
    let mut i: usize = position;
    let mut copies: usize = 0;
    let mut max_score: i32 = 0;
    loop {
        if i + motif.len() > sequence.len() {
            break;
        }
        let window = &sequence[i..i + motif.len()];
        let score = alignment::align(motif, window);
        if score >= threshold {
            copies += 1;
            max_score = cmp::max(score, max_score);
            i += motif.len();
        } else {
            break;
        }
    }
    Repeat {
        motif: motif.to_owned(),
        position: position,
        length: copies * motif.len(),
        copies: copies,
        score: max_score,
    }
}

// generate all motifs for a given sequence and max length of motif
fn generate_motifs(max_length: usize, sequence: &str) -> HashSet<String> {
    let mut motifs: HashSet<String> = HashSet::new();
    for i in 0..sequence.len() {
        let max_j = std::cmp::min(max_length, sequence.len() - i);
        for j in 1..=max_j {
            motifs.insert(sequence[i..i + j].to_owned());
        }
    }
    motifs
}
