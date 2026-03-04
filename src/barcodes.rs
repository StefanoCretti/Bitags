//! Single read barcoding functionalities (batching, bitap matching, overlap resolution)

use crate::fastq_io::{Writeable, get_buff_writer};
use crate::tags::Tag;

use itertools::Itertools;
use std::io::Write;
use std::{io, path};

const MAX_TAG_MISM: usize = 2; // Maximum number of mismatches per tag
const MAX_TAG_HITS: usize = 64; // Maximum number of tags aligned per read

type AlignedTags<'a> = Vec<(usize, &'a Tag)>;

/// Look for a tag in a DNA/cDNA read using the bitap algorithm
///
/// NOTE: Currently tag lenght must be 32 bp or lower
fn match_pattern(sequence: &[u8], tag: &Tag) -> Option<usize> {
    // Setup initial arrays (1 for exact match, 1 for each allowed mismatch).
    // NOTE: this allocation seems to be trivial speed-wise.
    let mut state: [u32; MAX_TAG_MISM + 1] = [!1u32; MAX_TAG_MISM + 1];
    let match_mask: u32 = 1u32 << tag.len;

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

/// Return the size of the non-overlapping span of a set of tags.
///
/// If the tags do not overlap, return the span (sum of all tag lengths).
/// As soon as an overlap is found, return zero.
fn comb_span(aligned_tags: &AlignedTags, subset: &[usize]) -> usize {
    let mut prev_end = 0;
    let mut tot_span = 0;

    for &ind in subset {
        let (pos, tag) = &aligned_tags[ind];
        if *pos < prev_end {
            return 0; // There is an overlap
        }
        prev_end = pos + tag.len;
        tot_span += tag.len;
    }
    tot_span
}

/// Remove the tags which do not belong to the best subset and return them
///
/// The vector of indices contains the indices of the aligned tags to keep.
/// It is assumed that indices are sorted in ascending order.
/// Sorting is also kept for both input and output tags.
fn pop_not_indexed<'a>(aligned_tags: &mut AlignedTags<'a>, indices: Vec<usize>) -> AlignedTags<'a> {
    let mut popped_tags = Vec::with_capacity(aligned_tags.len() - indices.len());
    let mut index_seq = indices.into_iter().peekable();

    let mut read_pos = 0;
    let mut write_pos = 0;
    while read_pos < aligned_tags.len() {
        if index_seq.peek() == Some(&read_pos) {
            if write_pos != read_pos {
                aligned_tags.swap(write_pos, read_pos);
            }
            write_pos += 1;
            index_seq.next();
        }
        read_pos += 1;
    }

    // Remove excess elements and collect them
    popped_tags.extend(aligned_tags.drain(write_pos..));
    popped_tags
}

/// Ensure that the barcode does not contain overlaps.
///
/// The problem is defined as finding the maximal non-overlapping set of spans.
/// Returned solution is guaranteed to be the optimal. If multiple solutions
/// have the same span, return the one favoring tag positions close to the
/// read start.
///
/// NOTE: assumes that tags are sorted by alignment start position.
fn remove_overlaps<'a>(aligned_tags: &mut AlignedTags<'a>) -> AlignedTags<'a> {
    // NOTE: Dynamic programming approach is technically algorithmically better
    // but it saves a negligible amount of time, so sticking to this for clarity.

    let num_tags: usize = aligned_tags.len();
    if num_tags == 0 {
        return Vec::with_capacity(0);
    }

    let mut num_removed: usize = 0;
    let mut best_comb: Vec<usize> = Vec::with_capacity(num_tags);
    let mut best_span: usize = 0;

    while num_removed < num_tags {
        // Test ALL combinations for a set size to guarantee optimal solution.
        for comb in (0..num_tags).combinations(num_tags - num_removed) {
            let span = comb_span(aligned_tags, &comb);
            if span > best_span {
                best_span = span;
                best_comb = comb;
            }
        }

        if best_span > 0 {
            break; // At least one non overlapping set was found
        }

        num_removed += 1;
    }

    // Actually remove overlapping tags and return them for logging purposes
    pop_not_indexed(aligned_tags, best_comb)
}

/// Utility function to convert a collection of aligned tags to a loggable string.
fn tags_to_string(aligned_tags: &AlignedTags) -> String {
    aligned_tags
        .iter()
        .map(|(pos, tag)| format!("{}:{}:{}", tag.get_seq(), tag.get_info(), pos))
        .join("|")
}

/// Reusable struct to find tags in a read sequence.
pub struct BarcodeBuilder<'a> {
    tags: &'a Vec<Tag>,
    hits: Vec<(usize, &'a Tag)>,
    logs: io::BufWriter<Writeable>,
}

impl<'a> BarcodeBuilder<'a> {
    /// Initialize a new instance of the struct from a vector of tags.
    pub fn new(tags: &'a Vec<Tag>, logs_file: path::PathBuf) -> Self {
        Self {
            tags,
            hits: Vec::with_capacity(MAX_TAG_HITS),
            logs: get_buff_writer(logs_file, true).unwrap(),
        }
    }

    /// Find all tags in a read and resolve potential overlaps.
    pub fn get_barcode(&mut self, seq: &[u8]) -> &Vec<(usize, &'a Tag)> {
        self.hits.clear();
        let seq_len = seq.len();

        // Find all matching tags
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

        // Sanitize the barcode
        self.hits.sort_by_key(|(pos, _)| *pos);
        let removed_tags = remove_overlaps(&mut self.hits);
        if !removed_tags.is_empty() {
            writeln!(
                self.logs,
                "{}\t{}",
                tags_to_string(&self.hits),
                tags_to_string(&removed_tags)
            )
            .unwrap();
        }

        return &self.hits;
    }
}
