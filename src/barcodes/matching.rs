//! Bitap fuzzy pattern matching for DNA/cDNA reads.

use super::MAX_TAG_MISM;
use crate::tags::Tag;

/// Look for a tag in a DNA/cDNA read using the bitap algorithm.
///
/// NOTE: Currently tag length must be 63 bp or lower.
pub(super) fn match_pattern(sequence: &[u8], tag: &Tag) -> Option<usize> {
    let mut state: [u64; MAX_TAG_MISM + 1] = [!1u64; MAX_TAG_MISM + 1];
    let match_mask: u64 = 1u64 << tag.len;

    for (i, base) in sequence.iter().enumerate() {
        let mut old_state = state[0];
        let base_patt = tag.patterns.get(&base);
        state[0] |= base_patt;
        state[0] <<= 1u8;
        // Loop unrolling cuts time in half at the cost of readability.
        match tag.max_mism {
            0 => {
                if (state[0] & match_mask) == 0 {
                    return Some(1 + i - tag.len);
                }
            }
            1 => {
                state[1] = (old_state & (state[1] | base_patt)) << 1;
                if (state[1] & match_mask) == 0 {
                    return Some(1 + i - tag.len);
                }
            }
            2 => {
                let tmp_state = state[1];
                state[1] = (old_state & (state[1] | base_patt)) << 1;
                state[2] = (tmp_state & (state[2] | base_patt)) << 1;
                if (state[2] & match_mask) == 0 {
                    return Some(1 + i - tag.len);
                }
            }
            _ => {
                for m in 1..=tag.max_mism {
                    let tmp_state = state[m];
                    state[m] = (old_state & (state[m] | base_patt)) << 1;
                    old_state = tmp_state;
                }
                if (state[tag.max_mism] & match_mask) == 0 {
                    return Some(1 + i - tag.len);
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::match_pattern;
    use crate::tags::Tag;

    fn make_tag(seq: &str, max_mism: usize) -> Tag {
        Tag::new(seq, "", max_mism)
    }

    // ── exact matching (max_mism = 0) ─────────────────────────────────────────

    #[test]
    fn exact_match_at_start() {
        let tag = make_tag("ACGT", 0);
        assert_eq!(match_pattern(b"ACGTTTTT", &tag), Some(0));
    }

    #[test]
    fn exact_match_in_middle() {
        let tag = make_tag("ACGT", 0);
        assert_eq!(match_pattern(b"TTTTACGTTTTT", &tag), Some(4));
    }

    #[test]
    fn exact_match_at_end() {
        let tag = make_tag("ACGT", 0);
        assert_eq!(match_pattern(b"TTTTACGT", &tag), Some(4));
    }

    #[test]
    fn no_match() {
        let tag = make_tag("ACGT", 0);
        assert_eq!(match_pattern(b"TTTTTTTT", &tag), None);
    }

    #[test]
    fn sequence_shorter_than_tag() {
        let tag = make_tag("ACGTACGT", 0);
        assert_eq!(match_pattern(b"ACG", &tag), None);
    }

    #[test]
    fn returns_first_match() {
        // tag appears twice; first occurrence is returned
        let tag = make_tag("ACGT", 0);
        assert_eq!(match_pattern(b"ACGTTTACGT", &tag), Some(0));
    }

    // ── 1-mismatch matching ───────────────────────────────────────────────────

    #[test]
    fn one_mismatch_found() {
        // ACGT vs ACTT — 1 mismatch at position 2
        let tag = make_tag("ACGT", 1);
        assert_eq!(match_pattern(b"ACTT", &tag), Some(0));
    }

    #[test]
    fn one_mismatch_not_allowed() {
        // max_mism=0, 1 mismatch → no match
        let tag = make_tag("ACGT", 0);
        assert_eq!(match_pattern(b"TTACTTTTTT", &tag), None);
    }

    // ── 2-mismatch matching ───────────────────────────────────────────────────

    #[test]
    fn two_mismatches_found() {
        // ACGT vs AATT — 2 mismatches at positions 1 and 2
        let tag = make_tag("ACGT", 2);
        assert_eq!(match_pattern(b"AATT", &tag), Some(0));
    }

    #[test]
    fn two_mismatches_not_allowed() {
        // max_mism=1, 2 mismatches → no match
        let tag = make_tag("ACGT", 1);
        assert_eq!(match_pattern(b"AATT", &tag), None);
    }

    // ── N in tag (always a mismatch) ──────────────────────────────────────────

    #[test]
    fn n_in_tag_counts_as_mismatch() {
        // ANCG: N always mismatches, so vs AACG there is exactly 1 mismatch
        let tag = make_tag("ANCG", 1);
        assert_eq!(match_pattern(b"AACG", &tag), Some(0));
    }

    #[test]
    fn n_in_tag_not_allowed() {
        // N always mismatches; max_mism=0 → no match even if rest is identical
        let tag = make_tag("ANCG", 0);
        assert_eq!(match_pattern(b"AACG", &tag), None);
    }

    // ── Truncated tag at read sides ───────────────────────────────────────────

    #[test]
    fn truncated_tag_start() {
        // Since the first base of the tag is truncated, not mismatched, the
        // algorithm should not find the tag even allowing for a mismatch.
        let tag = make_tag("AACG", 1);
        assert_eq!(match_pattern(b"ACGTTTT", &tag), None);
    }

    #[test]
    fn truncated_tag_end() {
        // Since the last base of the tag is truncated, not mismatched, the
        // algorithm should not find the tag even allowing for a mismatch.
        let tag = make_tag("AACG", 1);
        assert_eq!(match_pattern(b"TTTTAAC", &tag), None);
    }
}
