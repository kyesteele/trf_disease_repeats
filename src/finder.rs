use indicatif::{ProgressBar, ProgressStyle};

// structure of tandem repeat
#[derive(Clone, Debug)]
pub struct Repeat {
    pub start: usize,
    pub end: usize,
    pub period_size: usize,
    pub copy_number: f32,
    pub score: i32,
    pub a_percent: f32,
    pub c_percent: f32,
    pub g_percent: f32,
    pub t_percent: f32,
    pub sequence: String,
}

// for percentages of each base appearance
#[inline]
fn calc_composition(seq: &[u8]) -> (f32, f32, f32, f32) {
    let mut a = 0u32;
    let mut c = 0u32;
    let mut g = 0u32;
    let mut t = 0u32;

    for &base in seq {
        match base.to_ascii_uppercase() {
            b'A' => a += 1,
            b'C' => c += 1,
            b'G' => g += 1,
            b'T' => t += 1,
            _ => {}
        }
    }

    let total = (a + c + g + t) as f32;
    if total == 0.0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    (
        (a as f32 / total) * 100.0,
        (c as f32 / total) * 100.0,
        (g as f32 / total) * 100.0,
        (t as f32 / total) * 100.0,
    )
}

// check if repeat exists at given position
fn detect_repeat_at_position(
    seq: &[u8],
    start: usize,
    period: usize,
    match_weight: i32,
    mismatch_penalty: i32,
    min_score: i32,
) -> Option<Repeat> {
    if start + period * 2 > seq.len() {
        return None;
    }

    let pattern = &seq[start..start + period];
    let mut end = start + period;
    let mut num_copies = 1;

    let min_match_score = (period as i32 * match_weight * 2) / 3;

    while end + period <= seq.len() {
        let next_copy = &seq[end..end + period];

        // Quick match count
        let mut matches = 0;
        for i in 0..period {
            if pattern[i] == next_copy[i] {
                matches += 1;
            }
        }

        let score =
            (matches as i32 * match_weight) - ((period - matches) as i32 * mismatch_penalty);

        if score < min_match_score {
            break;
        }

        num_copies += 1;
        end += period;
    }

    if num_copies < 2 {
        return None;
    }

    let total_score = num_copies as i32 * period as i32 * match_weight / 2;

    if total_score < min_score {
        return None;
    }

    let repeat_seq = &seq[start..end];
    let (a, c, g, t) = calc_composition(repeat_seq);

    Some(Repeat {
        start: start + 1, // 1-indexed
        end,
        period_size: period,
        copy_number: num_copies as f32,
        score: total_score,
        a_percent: a,
        c_percent: c,
        g_percent: g,
        t_percent: t,
        sequence: String::from_utf8_lossy(repeat_seq).to_string(),
    })
}

// scan entire sequence for repeats
pub fn scan_sequence_trf(
    sequence: &str,
    match_weight: i32,
    mismatch_penalty: i32,
    _indel_penalty: i32,
    min_score: i32,
    max_period: usize,
    pb: &ProgressBar,
) -> Vec<Repeat> {
    let seq_bytes = sequence.as_bytes();
    let seq_len = seq_bytes.len();
    let mut repeats = Vec::new();

    // update in chunks
    let chunk_size = 10000.max(seq_len / 100);
    let mut i = 0;

    while i < seq_len {
        let chunk_end = (i + chunk_size).min(seq_len);

        while i < chunk_end {
            let mut best_repeat: Option<Repeat> = None;

            // try different period sizes, but prioritize smaller ones
            for period in 1..=max_period.min(seq_len - i).min(2000) {
                if i + period * 2 > seq_len {
                    break;
                }

                if let Some(repeat) = detect_repeat_at_position(
                    seq_bytes,
                    i,
                    period,
                    match_weight,
                    mismatch_penalty,
                    min_score,
                ) {
                    if best_repeat.is_none() || repeat.score > best_repeat.as_ref().unwrap().score {
                        best_repeat = Some(repeat);
                    }
                }
            }

            if let Some(repeat) = best_repeat {
                i = repeat.end;
                repeats.push(repeat);
            } else {
                i += 1;
            }
        }

        pb.set_position(i as u64);
    }

    pb.finish_with_message("Sequence scanned.");
    repeats
}

// print results in a table with summary
pub fn print_trf_output(sequence: &str, params: TrfParams) {
    let pb = ProgressBar::new(sequence.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    let repeats = scan_sequence_trf(
        sequence,
        params.match_weight,
        params.mismatch_penalty,
        params.indel_penalty,
        params.min_score,
        params.max_period,
        &pb,
    );

    println!(
        "Parameters: {} {} {} {} {} ",
        params.match_weight,
        params.mismatch_penalty,
        params.indel_penalty,
        params.min_score,
        params.max_period
    );
    println!();

    if repeats.is_empty() {
        println!("No tandem repeats found.");
        return;
    }

    println!("Indices\tPeriod\tCopies\tConsensus\tMatches\tIndels\tScore\tA%\tC%\tG%\tT%\tEntropy\tConsensus\tSequence");

    for r in &repeats {
        println!(
            "{}-{}\t{}\t{:.1}\t{}\t{:.0}\t{:.0}\t{:.0}\t{:.0}\t{}",
            r.start,
            r.end,
            r.period_size,
            r.copy_number,
            r.score,
            r.a_percent,
            r.c_percent,
            r.g_percent,
            r.t_percent,
            if r.sequence.len() > 50 {
                &r.sequence[..50]
            } else {
                &r.sequence
            }
        );
    }

    println!("\n== Summary ==");
    println!("Total repeats found: {}", repeats.len());

    let total_bases: usize = repeats.iter().map(|r| r.end - r.start + 1).sum();
    println!("Total bases in repeats: {}", total_bases);
    println!(
        "Percentage of sequence in repeats: {:.2}%",
        (total_bases as f32 / sequence.len() as f32) * 100.0
    );
}

pub struct TrfParams {
    pub match_weight: i32,
    pub mismatch_penalty: i32,
    pub indel_penalty: i32,
    pub min_score: i32,
    pub max_period: usize,
}

impl Default for TrfParams {
    fn default() -> Self {
        TrfParams {
            match_weight: 2,
            mismatch_penalty: 7,
            indel_penalty: 7,
            min_score: 50,
            max_period: 500,
        }
    }
}
