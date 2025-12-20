use crate::alignment::banded_smith_waterman;
use indicatif::{ProgressBar, ProgressStyle};
use std::cmp;

/// Structure representing a detected tandem repeat
#[derive(Clone, Debug)]
pub struct Repeat {
    pub start: usize,
    pub end: usize,
    pub period_size: usize,
    pub copy_number: f32,
    pub score: f32,
    pub a_percent: f32,
    pub c_percent: f32,
    pub g_percent: f32,
    pub t_percent: f32,
    pub sequence: String,
}

/// Parameters for TRF
pub struct TrfParams {
    pub match_weight: i32,
    pub mismatch_penalty: i32,
    pub indel_penalty: i32,
    pub min_score: i32,
    pub max_period: usize,
    pub prefilter_fraction: f32,
    pub max_copies: usize,
    pub refine_flank: usize,
    pub refine_band: usize,
}

impl Default for TrfParams {
    fn default() -> Self {
        TrfParams {
            match_weight: 2,
            mismatch_penalty: 7,
            indel_penalty: 7,
            min_score: 50,
            max_period: 500,
            prefilter_fraction: 0.75,
            max_copies: 1000,
            refine_flank: 100,
            refine_band: 8,
        }
    }
}

/// Calculate base composition (percentages)
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

/// Quick scoring function
fn quick_copy_score(motif: &[u8], copy: &[u8], match_weight: i32, mismatch_penalty: i32) -> i32 {
    let len = cmp::min(motif.len(), copy.len());
    let mut score = 0i32;
    for i in 0..len {
        if motif[i] == copy[i] {
            score += match_weight;
        } else {
            score -= mismatch_penalty;
        }
    }
    score
}

/// Phase 1: fast scan for candidate repeats
fn phase1_detect(
    seq: &[u8],
    start: usize,
    period: usize,
    params: &TrfParams,
) -> Option<(usize, usize, i32)> {
    let seq_len = seq.len();
    if start + period * 2 > seq_len {
        return None;
    }

    let motif = &seq[start..start + period];
    let mut end = start + period;
    let mut copies = 1usize;
    let mut agg_score = 0i32;

    let min_exact = ((period as f32) * params.prefilter_fraction).ceil() as usize;

    // extend motif while next copy is sufficiently similar
    while end + period <= seq_len && copies < params.max_copies {
        let next_copy = &seq[end..end + period];
        let mut exact = 0usize;
        for i in 0..period {
            if motif[i] == next_copy[i] {
                exact += 1;
            }
        }
        if exact < min_exact {
            break;
        }

        let sc = quick_copy_score(
            motif,
            next_copy,
            params.match_weight,
            params.mismatch_penalty,
        );
        if sc < -params.mismatch_penalty * (period as i32 / 4) {
            break;
        }

        agg_score = agg_score.saturating_add(sc);
        copies += 1;
        end += period;
    }

    if copies < 2 {
        return None;
    }

    // compute total score
    let total_score = (copies as i32) * (period as i32) * params.match_weight / 2 + agg_score / 2;
    if total_score < params.min_score {
        return None;
    }

    Some((end, copies, total_score))
}

/// Refine a rough repeat with banded alignment
fn refine_repeat(
    seq: &[u8],
    rough_start: usize,
    rough_end: usize,
    period: usize,
    copies: usize,
    rough_score: i32,
    params: &TrfParams,
) -> Repeat {
    let seq_len = seq.len();
    let flank = params.refine_flank;
    let win_start = rough_start.saturating_sub(flank);
    let win_end = cmp::min(seq_len, rough_end + flank);
    let window = &seq[win_start..win_end];

    // create pattern for alignment
    let motif = &seq[rough_start..rough_start + period];
    let repeat_times = cmp::min(copies, 50);
    let mut pattern_repeated = Vec::with_capacity(period * repeat_times);
    for _ in 0..repeat_times {
        pattern_repeated.extend_from_slice(motif);
    }

    // boundaries to use with banded SW
    let (best_score, rel_start, rel_end) = banded_smith_waterman(
        window,
        &pattern_repeated,
        params.match_weight,
        params.mismatch_penalty,
        params.indel_penalty,
        params.refine_band,
    );

    let abs_start = win_start.saturating_add(rel_start);
    let abs_end = win_start.saturating_add(rel_end).saturating_add(1);

    let final_start = cmp::max(rough_start, abs_start);
    let final_end = cmp::min(seq_len, cmp::max(rough_end, abs_end));

    let repeat_seq = &seq[final_start..final_end];
    let (a, c, g, t) = calc_composition(repeat_seq);

    let copy_number = ((final_end - final_start) as f32 / period as f32) / seq_len as f32;
    let combined_score = (rough_score.saturating_add(best_score / 2)) as f32 / seq_len as f32;

    Repeat {
        start: final_start + 1,
        end: final_end,
        period_size: period,
        copy_number,
        score: combined_score,
        a_percent: a,
        c_percent: c,
        g_percent: g,
        t_percent: t,
        sequence: String::from_utf8_lossy(repeat_seq).to_string(),
    }
}

/// Scan sequence for tandem repeats
pub fn scan_sequence_trf(sequence: &str, params: &TrfParams, pb: &ProgressBar) -> Vec<Repeat> {
    let seq_bytes = sequence.as_bytes();
    let seq_len = seq_bytes.len();
    let mut repeats = Vec::new();
    let mut last_repeat_end = 0usize;
    let chunk_size = 10000.max(seq_len / 100);
    let mut i = 0usize;

    while i < seq_len {
        let chunk_end = (i + chunk_size).min(seq_len);
        while i < chunk_end {
            if i < last_repeat_end {
                i += 1;
                continue;
            }

            // find the best candidate at this position
            let mut best_candidate: Option<(usize, usize, i32, usize)> = None;
            let max_period = cmp::min(params.max_period, seq_len.saturating_sub(i) / 2);
            for period in 1..=max_period {
                if i + period * 2 > seq_len {
                    break;
                }
                if let Some((rough_end, copies, rough_score)) =
                    phase1_detect(seq_bytes, i, period, params)
                {
                    if best_candidate.is_none() || rough_score > best_candidate.as_ref().unwrap().2
                    {
                        best_candidate = Some((rough_end, copies, rough_score, period));
                    }
                }
            }

            // refine and store repeat if pass
            if let Some((rough_end, copies, rough_score, period)) = best_candidate {
                let rep =
                    refine_repeat(seq_bytes, i, rough_end, period, copies, rough_score, params);
                if rep.start.saturating_sub(1) >= last_repeat_end {
                    last_repeat_end = rep.end;
                    repeats.push(rep);
                    i = last_repeat_end;
                    continue;
                } else {
                    i += 1;
                    continue;
                }
            } else {
                i += 1;
            }
        }

        pb.set_position(i as u64);
    }

    pb.finish_with_message("Sequence scanned.");
    repeats
}

/// Summary struct
pub struct TrfSummary {
    pub num_repeats: usize,
    pub percent_repeats: f64,
}

/// Collect summary normalized by sequence length
pub fn collect_trf_summary(sequence: &str, params: TrfParams) -> TrfSummary {
    let pb = ProgressBar::new(sequence.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    let repeats = scan_sequence_trf(sequence, &params, &pb);

    let total_bases: usize = repeats
        .iter()
        .map(|r| r.end.saturating_sub(r.start - 1))
        .sum();
    let percent_repeats = if sequence.len() > 0 {
        (total_bases as f64 / sequence.len() as f64) * 100.0
    } else {
        0.0
    };

    TrfSummary {
        num_repeats: repeats.len(),
        percent_repeats,
    }
}
