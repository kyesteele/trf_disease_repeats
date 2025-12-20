// src/alignment.rs
use std::cmp;

/// Banded Smith–Waterman local alignment
/// Returns (best_score, start_j, end_j) relative to window_seq
pub fn banded_smith_waterman(
    window_seq: &[u8],
    pattern_repeated: &[u8],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
    band: usize,
) -> (i32, usize, usize) {
    let n = window_seq.len();
    let m = pattern_repeated.len();
    if n == 0 || m == 0 {
        return (0, 0, 0);
    }

    let band_width = 2 * band + 1;
    let mut dp_prev = vec![0i32; band_width];
    let mut dp_cur = vec![0i32; band_width];

    let mut best_score = 0i32;
    let mut best_j = 0usize;

    // iterate over the pattern
    for i in 1..=m {
        let diag_pos = ((i * n) as f64 / m as f64).round() as isize;
        let j_min = cmp::max(1isize, diag_pos - band as isize) as usize;
        let j_max = cmp::min(n, (diag_pos + band as isize) as usize);

        dp_cur.fill(0);

        // iterate over positions within band
        for j in j_min..=j_max {
            let band_idx = j as isize - diag_pos + band as isize;
            if band_idx < 0 || band_idx >= band_width as isize {
                continue;
            }
            let b = band_idx as usize;

            // diagonal (match/mismatch) score
            let diag = dp_prev[b]
                + if pattern_repeated[i - 1] == window_seq[j - 1] {
                    match_score
                } else {
                    -mismatch_penalty
                };

            // insertion
            let left = if b > 0 {
                dp_cur[b - 1] - gap_penalty
            } else {
                i32::MIN
            };

            // deletion
            let up = dp_prev[b] - gap_penalty;

            let val = cmp::max(0, cmp::max(diag, cmp::max(left, up)));
            dp_cur[b] = val;

            if val > best_score {
                best_score = val;
                best_j = j;
            }
        }

        std::mem::swap(&mut dp_prev, &mut dp_cur);
    }

    // approximate start position
    let mut start_j = best_j;
    while start_j > 0 && best_j - start_j <= band + 10 {
        start_j -= 1;
    }

    (best_score, start_j, best_j)
}
