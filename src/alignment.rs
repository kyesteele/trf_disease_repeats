use std::cmp;

pub fn align(motif: &[u8], window: &[u8]) -> i32 {
    // match +2, mismatch -2, gap -7
    let mut score_left: i32 = -7;
    let mut score_center: i32 = 0;
    let mut score_right: i32 = -7;

    // check -1 position
    for i in 1..motif.len() {
        if window[i - 1] == motif[i] {
            score_left += 2;
        } else {
            score_left -= 2;
        }
    }

    // check 0 position (gapless)
    for i in 0..motif.len() {
        if window[i] == motif[i] {
            score_center += 2;
        } else {
            score_center -= 2;
        }
    }

    // check +1 position
    for i in 0..(motif.len() - 1) {
        if window[i + 1] == motif[i] {
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
