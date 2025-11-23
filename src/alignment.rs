use std::cmp;

pub fn align(motif: &str, window: &str) -> i32 {
    // match +2, mismatch -2, gap -7
    let mut score_left: i32 = -7;
    let mut score_center: i32 = 0;
    let mut score_right: i32 = -7;
    let motif_bytes = motif.as_bytes();
    let window_bytes = window.as_bytes();

    // check -1 position
    for i in 1..motif.len() {
        if window_bytes[i - 1] == motif_bytes[i] {
            score_left += 2;
        } else {
            score_left -= 2;
        }
    }

    // check 0 position (gapless)
    for i in 0..motif.len() {
        if window_bytes[i] == motif_bytes[i] {
            score_center += 2;
        } else {
            score_center -= 2;
        }
    }

    // check +1 position
    for i in 0..(motif.len() - 1) {
        if window_bytes[i + 1] == motif_bytes[i] {
            score_right += 2;
        } else {
            score_right -= 2;
        }
    }

    // return best score
    let mut best_score = cmp::max(score_left, score_center);
    best_score = cmp::max(best_score, score_right);
    best_score
}
