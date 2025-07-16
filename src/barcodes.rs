use crate::tags::Tag;

const MAX_MISM: usize = 2;

fn match_pattern(sequence: &str, tag: &Tag) -> Option<u8> {
    // Moving this allocation outside of the function does not
    // seem to do much, so left here for now.
    let mut state: [u16; MAX_MISM + 1] = [!1u16; MAX_MISM + 1];
    let match_mask = 1u16 << tag.len;

    for (i, base) in sequence.bytes().enumerate() {
        let mut old_state = state[0];
        let base_patt = tag.patterns.get(base);

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

pub struct BarcodeBuilder {
    bitap_tags: Vec<(u8, Tag)>,
    hit_buffer: Vec<(u8, Tag)>,
    candidates: Vec<(u8, Tag)>,
}

impl BarcodeBuilder {
    pub fn new(tags: &[Tag]) -> Self {
        let bitap_tags: Vec<(u8, Tag)> = tags.into_iter().map(|tag| (0u8, tag.clone())).collect();

        // Though technically reachable in extremely unlikely scenarios,
        // this is a VERY generous upper bound to the number of hits.
        let vec_length: usize = bitap_tags.len();

        Self {
            bitap_tags,
            hit_buffer: Vec::with_capacity(vec_length),
            candidates: Vec::with_capacity(vec_length),
        }
    }

    pub fn get_barcode(&mut self, seq: &str) -> &Vec<(u8, Tag)> {
        self.hit_buffer.clear();
        self.candidates.clear();
        self.candidates.extend_from_slice(&self.bitap_tags); // Still not the best

        let seq_len: u8 = seq.bytes().len() as u8;

        while !self.candidates.is_empty() {
            self.candidates.retain_mut(|(pos, tag)| {
                if let Some(new_pos) = match_pattern(&seq[*pos as usize..], tag) {
                    let match_pos = *pos + new_pos;
                    self.hit_buffer.push((match_pos, tag.clone()));

                    *pos = match_pos + tag.len as u8;
                    *pos < seq_len
                } else {
                    false
                }
            });
        }

        // TODO: Handle overlapping tags
        // TODO: Handle terminal tags false positives

        self.hit_buffer.sort_by_key(|(pos, _)| *pos);
        return &self.hit_buffer;
    }
}
