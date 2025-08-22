//! Single read barcoding functionalities (batching, bitap matching, overlap resolution)

use crate::tags::Tag;

const MAX_TAG_MISM: usize = 2; // Maximum number of mismatches per tag
const MAX_TAG_HITS: usize = 64; // Maximum number of tags aligned per read

type AlignedTags<'a> = Vec<(usize, &'a Tag)>;

/// Look for a tag in a DNA/cDNA read using the bitap algorithm
/// NOTE: Currently tag lenght must be 16 bp or lower
fn match_pattern(sequence: &[u8], tag: &Tag) -> Option<usize> {
    // Setup initial arrays (1 for exact match, 1 for each allowed mismatch).
    // NOTE: this allocation seems to be trivial speed-wise.
    let mut state: [u16; MAX_TAG_MISM + 1] = [!1u16; MAX_TAG_MISM + 1];
    let match_mask: u16 = 1u16 << tag.len;

    for (i, base) in sequence.iter().enumerate() {
        let mut old_state = state[0];
        let base_patt = tag.patterns.get(&base);
        state[0] |= base_patt;
        state[0] <<= 1u8;
        // Loop unrolling does not look pretty at all, but it cuts time
        // in half so it is, unfortunately, worth the loss in readability.
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
                // With max mismatch set to 2 this branch should never happen
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
    return None;
}

/// Return a bool depending of the presence of overlaps in a set of tags
///
/// For each tag, check whether its starting position is bigger than the
/// previous tag's end. Returns as soon as one overlap is found.
fn does_set_overlap(aligned_tags: &mut AlignedTags) -> bool {
    let mut prev_pos: usize = aligned_tags[0].0 + aligned_tags[0].1.len;
    for ind in 1..aligned_tags.len() {
        if aligned_tags[ind].0 < prev_pos {
            return true;
        }
        prev_pos = aligned_tags[ind].0 + aligned_tags[ind].1.len;
    }
    return false;
}

/// Try to solve overlaps by removing a single tag.
///
/// In turns, try to remove each tag from the set. Select as solutions the
/// maximum non overlapping subset. If no one is found, return false.
fn leave_one_out(aligned_tags: &mut AlignedTags) -> bool {
    let mut rm_index: i8 = -1;
    let mut max_span: usize = 0;

    // Test all possible -1 subsets
    for ind in 0..aligned_tags.len() {
        // NOTE: This collect should not affect performance significantly, but
        // simplifies overlap checking significantly. Replace if needed.
        let mut subset: AlignedTags = aligned_tags
            .iter()
            .enumerate()
            .filter(|(idx, _)| *idx != ind)
            .map(|(_, val)| *val)
            .collect();

        // There are still overlaps, so it is not a valid subset
        if does_set_overlap(&mut subset) {
            continue;
        }

        let span: usize = subset.iter().map(|val| val.1.len).sum();

        if span > max_span {
            max_span = span;
            rm_index = ind as i8;
        }
    }

    // No match without overlap was found
    if rm_index == -1 {
        return false;
    }

    // Remove chosen tag
    aligned_tags.remove(rm_index as usize);
    return true;
}

/// Ensure that the barcode is valid (no overlapping tags)
///
/// Tags are removed from the barcode until there are no more overlaps.
/// The problem is defined as finding the maximal non-overlapping set of spans.
/// Returned solution is guaranteed to be (one of) the optimal.
///
/// TODO: Need to add logging for debugging and understanding potentially
/// problematic pairs of tags.
fn sanitize_barcode(aligned_tags: &mut AlignedTags) {
    // Since the zero overlap and singe overlap cases will be by far the most
    // common, those cases are handled directly with more efficient functions.
    // For other cases, uses a more general but less efficient approach.

    // Not a single overlap, expected to be vast majority of the cases
    if (aligned_tags.len() == 0) || !does_set_overlap(aligned_tags) {
        return;
    }

    // Return here if the problem could be solved by removing a single tag.
    if leave_one_out(aligned_tags) {
        return;
    }

    // TODO: Add proper removal of multiple overlaps, currently just give up on the tag.
    aligned_tags.clear();
}

/// Reusable struct to find tags in a read sequence.
///
/// TODO: Add some stats for diagnostic purposes (num barcoded, overlaps...)
pub struct BarcodeBuilder<'a> {
    tags: &'a Vec<Tag>,
    hits: Vec<(usize, &'a Tag)>,
}

impl<'a> BarcodeBuilder<'a> {
    /// Initialize a new instance of the struct from a vector of tags.
    pub fn new(tags: &'a Vec<Tag>) -> Self {
        Self {
            tags,
            hits: Vec::with_capacity(MAX_TAG_HITS),
        }
    }

    /// Find all tags in a read and resolve potential overlaps.
    pub fn get_barcode(&mut self, seq: &[u8]) -> &Vec<(usize, &'a Tag)> {
        self.hits.clear();
        let seq_len = seq.len();

        for tag in self.tags.iter() {
            let mut start_pos: usize = 0;
            while start_pos + tag.len <= seq_len {
                if let Some(relative_pos) = match_pattern(&seq[start_pos..], tag) {
                    // Returned position from match is relative to str slice start
                    let absolute_pos = start_pos + relative_pos;
                    self.hits.push((absolute_pos, tag));
                    start_pos = absolute_pos + tag.len; // Skip over current match
                } else {
                    break;
                }
            }
        }
        self.hits.sort_by_key(|(pos, _)| *pos);

        // Remove overlapping tags, if any
        sanitize_barcode(&mut self.hits);

        return &self.hits;
    }
}
