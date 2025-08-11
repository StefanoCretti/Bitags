use crate::tags::Tag;

const MAX_MISM: usize = 2;

fn match_pattern(sequence: &[u8], tag: &Tag) -> Option<u8> {
    // Moving this allocation outside of the function does not
    // seem to do much, so left here for now.
    let mut state: [u16; MAX_MISM + 1] = [!1u16; MAX_MISM + 1];
    let match_mask = 1u16 << tag.len;
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
                    return Some((1 + i - tag.len) as u8);
                }
            }
            1 => {
                state[1] = (old_state & (state[1] | base_patt)) << 1;
                if (state[1] & match_mask) == 0 {
                    return Some((1 + i - tag.len) as u8);
                }
            }
            2 => {
                let tmp_state = state[1];
                state[1] = (old_state & (state[1] | base_patt)) << 1;
                state[2] = (tmp_state & (state[2] | base_patt)) << 1;
                if (state[2] & match_mask) == 0 {
                    return Some((1 + i - tag.len) as u8);
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
                    return Some((1 + i - tag.len) as u8);
                }
            }
        }
    }
    return None;
}

// Return a bool depending of the presence of overlaps in a set of tags
//
// For each tag, check whether its starting position is bigger than the
// previous tag's end. Returns as soon as one overlap is found.
fn does_set_overlap(aligned_tags: &mut Vec<(u8, &Tag)>) -> bool {
    let mut prev_pos: u8 = aligned_tags[0].0 + aligned_tags[0].1.len as u8;
    for ind in 1..aligned_tags.len() {
        if aligned_tags[ind].0 < prev_pos {
            return true;
        }
        prev_pos = aligned_tags[ind].0 + aligned_tags[ind].1.len as u8;
    }
    return false;
}

// Try to solve overlaps by removing a single tag.
//
// In turns, try to remove each tag from the set. Select as solutions the
// maximum non overlapping subset. If no one is found, return false.
fn leave_one_out(aligned_tags: &mut Vec<(u8, &Tag)>) -> bool {
    let mut rm_index: i8 = -1;
    let mut max_span: usize = 0;

    // Test all possible -1 subsets
    for ind in 0..aligned_tags.len() {
        // NOTE: This collect should not affect performance significantly, but
        // simplifies overlap checking significantly. Replace if needed.
        let mut subset: Vec<(u8, &Tag)> = aligned_tags
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

fn sanitize_barcode(aligned_tags: &mut Vec<(u8, &Tag)>) {
    // Not a single overlap, expected to be vast majority of the cases
    if (aligned_tags.len() == 0) || !does_set_overlap(aligned_tags) {
        return;
    }

    // Returns here if the problem could be solved by removing a single tag.
    if leave_one_out(aligned_tags) {
        return;
    }

    // TODO: Add proper removal of multiple overlaps
    // This function will be way more computationally expensive, but it
    // should be called so rarely that it should not matter.
    println!("Unable to solve multiple overlaps, removing all tags from the read");
    aligned_tags.clear();
}

// Reusable struct to find tags in a read sequence.
//
// TODO: Maybe add some stats for diagnostic purposes (num barcoded, overlaps...)
pub struct BarcodeBuilder<'a> {
    bitap_tags: Vec<(u8, &'a Tag)>,
    hit_buffer: Vec<(u8, &'a Tag)>,
    candidates: Vec<(u8, &'a Tag)>,
}

impl<'a> BarcodeBuilder<'a> {
    pub fn new(tags: &'a [Tag]) -> Self {
        let bitap_tags: Vec<(u8, &Tag)> = tags.into_iter().map(|tag| (0u8, tag)).collect();

        // Though technically reachable in extremely unlikely scenarios,
        // this is a VERY generous upper bound to the number of hits.
        let vec_length: usize = bitap_tags.len();

        Self {
            bitap_tags,
            hit_buffer: Vec::with_capacity(vec_length),
            candidates: Vec::with_capacity(vec_length),
        }
    }

    // Find all tags in a read and resolve potential overlaps.
    pub fn get_barcode(&mut self, seq: &[u8]) -> &Vec<(u8, &'a Tag)> {
        self.hit_buffer.clear();
        self.candidates.clear();
        self.candidates.extend_from_slice(&self.bitap_tags); // Still not the best

        let seq_len: u8 = seq.len() as u8;

        while !self.candidates.is_empty() {
            self.candidates.retain_mut(|(pos, tag)| {
                if let Some(new_pos) = match_pattern(&seq[*pos as usize..], tag) {
                    let match_pos = *pos + new_pos;
                    self.hit_buffer.push((match_pos, tag));

                    *pos = match_pos + tag.len as u8;
                    *pos < seq_len
                } else {
                    false
                }
            });
        }
        self.hit_buffer.sort_by_key(|(pos, _)| *pos);

        // Remove overlapping tags, if any
        sanitize_barcode(&mut self.hit_buffer);

        return &self.hit_buffer;
    }
}
