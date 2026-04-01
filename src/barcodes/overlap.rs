//! Overlap resolution for aligned barcoding tags.

use super::{AlignedTags, MAX_TAG_HITS};

/// Move tags not in `indices` from `aligned_tags` into `removed`.
///
/// `indices` must be sorted ascending. `removed` is cleared before filling.
/// Discared tags are not guaranteed to be sorted by alignement position.
fn pop_not_indexed<'a>(
    aligned_tags: &mut AlignedTags<'a>,
    indices: &[usize],
    removed: &mut AlignedTags<'a>,
) {
    removed.clear();
    let mut index_seq = indices.iter().peekable();
    let mut write_pos = 0;
    for read_pos in 0..aligned_tags.len() {
        if index_seq.peek() == Some(&&read_pos) {
            if write_pos != read_pos {
                aligned_tags.swap(write_pos, read_pos);
            }
            write_pos += 1;
            index_seq.next();
        }
    }
    removed.extend(aligned_tags.drain(write_pos..));
}

/// Ensure the barcode has no overlapping tags, using three phases:
///
/// **Phase 1 — `nxt[i]`**: for each tag i (start-sorted), binary-search for
/// the first index j such that `tags[j]` starts at or after `tags[i]` ends.
///
/// **Phase 2 — suffix DP**: `sdp[i]` = optimal `(span, count)` achievable
/// using only `tags[i..]`, processed in reverse start order. Optimises:
///   (1) maximum total span, then (2) maximum number of tags.
///
/// **Phase 3 — greedy forward pass**: iterate tags in start order; include
/// the earliest compatible tag that keeps the remaining optimum reachable.
/// This greedy choice yields the lexicographically smallest start-position
/// sequence among all solutions with the same `(span, count)` — i.e. at the
/// first position where two solutions differ, pick the tag with the smaller
/// start.
///
/// All buffers are stack-allocated (no heap allocation).
/// Removed tags are written into `removed` (cleared before use).
///
/// Requires `aligned_tags` sorted by start position on entry.
pub(super) fn remove_overlaps<'a>(
    aligned_tags: &mut AlignedTags<'a>,
    removed: &mut AlignedTags<'a>,
) {
    removed.clear();
    let num_tags = aligned_tags.len();
    if num_tags == 0 {
        return;
    }

    // Fast path: tags are start-sorted; no overlap iff every adjacent pair is gap-free.
    if aligned_tags
        .windows(2)
        .all(|w| w[0].0 + w[0].1.len <= w[1].0)
    {
        return;
    }

    // --- Phase 1: nxt[i] ---
    // nxt[i] = first index j > i such that tags[j].start >= tags[i].end,
    //          or num_tags if no such tag exists.
    let mut nxt_buf = [0usize; MAX_TAG_HITS];
    let nxt = &mut nxt_buf[..num_tags];
    for i in 0..num_tags {
        let end_i = aligned_tags[i].0 + aligned_tags[i].1.len;
        let mut lo = i + 1;
        let mut hi = num_tags;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if aligned_tags[mid].0 >= end_i {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        nxt[i] = lo;
    }

    // --- Phase 2: suffix DP ---
    // sdp[i] = (span, count): best achievable from tags[i..] in start order.
    // sdp[num_tags] = (0, 0) by initialisation.
    let mut sdp_buf = [(0usize, 0usize); MAX_TAG_HITS + 1];
    let sdp = &mut sdp_buf[..num_tags + 1];
    for i in (0..num_tags).rev() {
        let tag_len = aligned_tags[i].1.len;
        let skip = sdp[i + 1];
        let (base_span, base_count) = sdp[nxt[i]];
        let include = (base_span + tag_len, base_count + 1);
        sdp[i] = if include.0 > skip.0 || (include.0 == skip.0 && include.1 > skip.1) {
            include
        } else {
            skip
        };
    }

    // --- Phase 3: greedy forward pass ---
    // At each step, pick the earliest-starting compatible tag that keeps the
    // remaining (span, count) target reachable.
    let (mut rem_span, mut rem_count) = sdp[0];
    let mut pos = 0usize;
    let mut sel_buf = [0usize; MAX_TAG_HITS];
    let mut sel_len = 0usize;

    for i in 0..num_tags {
        if rem_count == 0 {
            break;
        }
        let (tag_start, tag) = &aligned_tags[i];
        if *tag_start < pos {
            continue;
        }
        let (after_span, after_count) = sdp[nxt[i]];
        if after_span + tag.len == rem_span && after_count + 1 == rem_count {
            sel_buf[sel_len] = i;
            sel_len += 1;
            rem_span -= tag.len;
            rem_count -= 1;
            pos = tag_start + tag.len;
        }
    }

    pop_not_indexed(aligned_tags, &sel_buf[..sel_len], removed);
}

#[cfg(test)]
mod tests {
    use super::super::{AlignedTags, MAX_TAG_HITS};
    use super::remove_overlaps;
    use crate::tags::Tag;

    // ── helpers ──────────────────────────────────────────────────────────────

    /// Build a Tag whose sequence is `seq_len` 'A's and whose info is `info`.
    /// max_mism is irrelevant for overlap-removal tests.
    fn make_tag(seq_len: usize, info: &str) -> Tag {
        Tag::new(&"A".repeat(seq_len), info, 0)
    }

    /// Thin wrapper around remove_overlaps for readability.
    fn resolve<'a>(hits: &mut AlignedTags<'a>) -> AlignedTags<'a> {
        let mut removed = Vec::with_capacity(MAX_TAG_HITS);
        remove_overlaps(hits, &mut removed);
        removed
    }

    /// Return the info strings of every tag in the slice, in order.
    fn infos<'a>(tags: &'a AlignedTags<'a>) -> Vec<&'a str> {
        tags.iter().map(|(_, t)| t.get_info()).collect()
    }

    // ── fast-path cases ───────────────────────────────────────────────────────

    #[test]
    fn no_tags() {
        let mut hits: AlignedTags = vec![];
        let removed = resolve(&mut hits);
        assert!(hits.is_empty());
        assert!(removed.is_empty());
    }

    #[test]
    fn single_tag() {
        let t = make_tag(8, "T1");
        let mut hits: AlignedTags = vec![(0, &t)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["T1"]);
        assert!(removed.is_empty());
    }

    #[test]
    fn adjacent_tags() {
        // [T1 @ 0..8]  [T2 @ 8..16] — adjacent, non-overlapping.
        let t1 = make_tag(8, "T1");
        let t2 = make_tag(8, "T2");
        let mut hits: AlignedTags = vec![(0, &t1), (8, &t2)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["T1", "T2"]);
        assert!(removed.is_empty());
    }

    #[test]
    fn no_overlap_three_tags() {
        // [T1 @ 0..5]  [T2 @ 5..10]  [T3 @ 10..15]
        let t1 = make_tag(5, "T1");
        let t2 = make_tag(5, "T2");
        let t3 = make_tag(5, "T3");
        let mut hits: AlignedTags = vec![(0, &t1), (5, &t2), (10, &t3)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["T1", "T2", "T3"]);
        assert!(removed.is_empty());
    }

    // ── tiebreak 1: maximum span ──────────────────────────────────────────────

    #[test]
    fn larger_tag_wins() {
        // [SHORT @ 0..3] overlaps [LONG @ 0..10].
        // LONG has greater span (10 > 3), so LONG is kept.
        let short = make_tag(3, "SHORT");
        let long = make_tag(10, "LONG");
        let mut hits: AlignedTags = vec![(0, &short), (0, &long)];
        hits.sort_by_key(|(pos, _)| *pos);
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["LONG"]);
        assert_eq!(infos(&removed), ["SHORT"]);
    }

    #[test]
    fn higher_span_wins() {
        // [A @ 0..10]  [B @ 0..6]  [C @ 9..15]
        // A overlaps B (both start at 0) and C (A ends at 10 > C's start 9).
        // B and C do not overlap (B ends at 6 ≤ C's start 9).
        // Solutions: [A] → span=10; [B,C] → span=12.
        // B+C has greater span, so it wins outright (no tiebreak needed).
        let a = make_tag(10, "A");
        let b = make_tag(6, "B");
        let c = make_tag(6, "C");
        let mut hits: AlignedTags = vec![(0, &a), (0, &b), (9, &c)];
        hits.sort_by_key(|(pos, _)| *pos);
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["B", "C"]);
        assert_eq!(infos(&removed), ["A"]);
    }

    #[test]
    fn higher_combined_span_wins() {
        // [K1 @ 0..15]  [K2 @ 15..30]  [C @ 0..10]  [D @ 10..20]  [E @ 20..29]
        //
        // Overlaps:
        //   K1 overlaps C (both at 0) and D (K1 ends 15 > D starts 10).
        //   K2 overlaps D (K2 starts 15 < D ends 20) and E (K2 starts 15 < E ends 29).
        //   C, D, E are mutually non-overlapping.
        //
        // Solutions: [K1,K2] → span=30; [C,D,E] → span=29.
        // [K1,K2] wins outright on span (30 > 29) despite having fewer tags.
        let k1 = make_tag(15, "K1");
        let k2 = make_tag(15, "K2");
        let c = make_tag(10, "C");
        let d = make_tag(10, "D");
        let e = make_tag(9, "E");
        let mut hits: AlignedTags = vec![(0, &k1), (0, &c), (10, &d), (15, &k2), (20, &e)];
        hits.sort_by_key(|(pos, _)| *pos);
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["K1", "K2"]);
        let mut r = infos(&removed);
        r.sort_unstable();
        assert_eq!(r, ["C", "D", "E"]);
    }

    // ── tiebreak 2: maximum number of tags ────────────────────────────────────

    #[test]
    fn more_tags_win() {
        // A @0..10, B @0..5, C @5..10 — A overlaps B and C; B and C don't overlap.
        // B+C: span=10, count=2.  A alone: span=10, count=1.
        // Tiebreak 2: more tags → B+C wins.
        let a = make_tag(10, "A");
        let b = make_tag(5, "B");
        let c = make_tag(5, "C");
        let mut hits: AlignedTags = vec![(0, &a), (0, &b), (5, &c)];
        hits.sort_by_key(|(pos, _)| *pos);
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["B", "C"]);
        assert_eq!(infos(&removed), ["A"]);
    }

    #[test]
    fn more_tags_win_three_over_two() {
        // [P @ 0..15]  [Q @ 15..30]  [X @ 0..10]  [Y @ 10..20]  [Z @ 20..30]
        //
        // Overlaps:
        //   P overlaps X (both at 0) and Y (P ends 15 > Y starts 10).
        //   Q overlaps Y (Q starts 15 < Y ends 20) and Z (Q starts 15 < Z ends 30).
        //   X, Y, Z are mutually non-overlapping.
        //
        // Solutions: [P,Q] → span=30, count=2; [X,Y,Z] → span=30, count=3.
        // Equal span → tiebreak 2: more tags wins → [X,Y,Z] wins.
        let p = make_tag(15, "P");
        let q = make_tag(15, "Q");
        let x = make_tag(10, "X");
        let y = make_tag(10, "Y");
        let z = make_tag(10, "Z");
        let mut hits: AlignedTags = vec![(0, &p), (0, &x), (10, &y), (15, &q), (20, &z)];
        hits.sort_by_key(|(pos, _)| *pos);
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["X", "Y", "Z"]);
        let mut r = infos(&removed);
        r.sort_unstable();
        assert_eq!(r, ["P", "Q"]);
    }

    // ── tiebreak 3: lexicographic start position ──────────────────────────────

    #[test]
    fn lex_overlap_is_second_tag() {
        // [T1 @ 0..5]  [T2 @ 5..10]  [T3 @ 8..13] — overlap between T2 and T3.
        // Solutions: [T1,T2] or [T1,T3], both span=10, count=2.
        // T2 starts at 5, T3 starts at 8 → T2 < T3 → keep [T1, T2].
        let t1 = make_tag(5, "T1");
        let t2 = make_tag(5, "T2");
        let t3 = make_tag(5, "T3");
        let mut hits: AlignedTags = vec![(0, &t1), (5, &t2), (8, &t3)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["T1", "T2"]);
        assert_eq!(infos(&removed), ["T3"]);
    }

    #[test]
    fn lex_overlap_is_first_tag() {
        // [T1 @ 0..5]  [T2 @ 3..8]  [T3 @ 8..13] — overlap between T1 and T2.
        // Solutions: [T1,T3] or [T2,T3], both span=10, count=2.
        // T1 starts at 0, T2 starts at 3 → T1 < T2 → keep [T1, T3].
        let t1 = make_tag(5, "T1");
        let t2 = make_tag(5, "T2");
        let t3 = make_tag(5, "T3");
        let mut hits: AlignedTags = vec![(0, &t1), (3, &t2), (8, &t3)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["T1", "T3"]);
        assert_eq!(infos(&removed), ["T2"]);
    }

    #[test]
    fn lex_overlap_is_middle_tag() {
        // [T1 @ 0..5]  [T2 @ 5..10]  [T3 @ 8..13]  [T4 @ 13..18]
        // Overlap between T2 and T3 only; T4 is compatible with both.
        // Solutions: [T1,T2,T4] or [T1,T3,T4], both span=15, count=3.
        // T2 starts at 5, T3 starts at 8 → T2 < T3 → keep [T1, T2, T4].
        let t1 = make_tag(5, "T1");
        let t2 = make_tag(5, "T2");
        let t3 = make_tag(5, "T3");
        let t4 = make_tag(5, "T4");
        let mut hits: AlignedTags = vec![(0, &t1), (5, &t2), (8, &t3), (13, &t4)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["T1", "T2", "T4"]);
        assert_eq!(infos(&removed), ["T3"]);
    }

    // ── complex cases ─────────────────────────────────────────────────────────

    #[test]
    fn multiple_overlap_groups() {
        // Two independent overlap pairs.
        // Group 1: [T1 @ 0..5] vs [T2 @ 3..8] — T1 starts earlier.
        // Group 2: [T3 @ 10..15] vs [T4 @ 13..18] — T3 starts earlier.
        // Optimal: [T1, T3], span=10, count=2.
        let t1 = make_tag(5, "T1");
        let t2 = make_tag(5, "T2");
        let t3 = make_tag(5, "T3");
        let t4 = make_tag(5, "T4");
        let mut hits: AlignedTags = vec![(0, &t1), (3, &t2), (10, &t3), (13, &t4)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["T1", "T3"]);
        assert_eq!(infos(&removed), ["T2", "T4"]);
    }

    #[test]
    fn chain_overlap() {
        // Tags overlap in a chain: A-B, B-C, C-D.
        // [A @ 0..10]  [B @ 8..18]  [C @ 16..26]  [D @ 24..34]
        // {A,C}: span=20, {A,D}: span=20, {B,D}: span=20 — all tied at span=20, count=2.
        // Lex: A(0) < B(8) → A goes first. After A (ends 10), C(16) comes before D(24).
        // → keep [A, C].
        let a = make_tag(10, "A");
        let b = make_tag(10, "B");
        let c = make_tag(10, "C");
        let d = make_tag(10, "D");
        let mut hits: AlignedTags = vec![(0, &a), (8, &b), (16, &c), (24, &d)];
        let removed = resolve(&mut hits);
        assert_eq!(infos(&hits), ["A", "C"]);
        assert_eq!(infos(&removed), ["B", "D"]);
    }
}
