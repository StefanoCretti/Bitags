//! BarcodeBuilder: finds and resolves tags in a single read.

use crate::fastq_io::{Writeable, get_buff_writer};
use crate::tags::Tag;

use itertools::Itertools;
use std::io::Write;
use std::{io, path};

use super::matching::match_pattern;
use super::overlap::remove_overlaps;
use super::{AlignedTags, MAX_TAG_HITS};

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
    removed_hits: Vec<(usize, &'a Tag)>,
    logs: io::BufWriter<Writeable>,
}

impl<'a> BarcodeBuilder<'a> {
    /// Initialize a new instance of the struct from a vector of tags.
    pub fn new(tags: &'a Vec<Tag>, logs_file: path::PathBuf) -> Self {
        Self {
            tags,
            hits: Vec::with_capacity(MAX_TAG_HITS),
            removed_hits: Vec::with_capacity(MAX_TAG_HITS),
            logs: get_buff_writer(logs_file, true).unwrap(),
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
                    let absolute_pos = start_pos + relative_pos;
                    self.hits.push((absolute_pos, tag));
                    start_pos = absolute_pos + tag.len;
                } else {
                    break;
                }
            }
        }

        self.hits.sort_by_key(|(pos, _)| *pos);
        remove_overlaps(&mut self.hits, &mut self.removed_hits);
        if !self.removed_hits.is_empty() {
            writeln!(
                self.logs,
                "{}\t{}",
                tags_to_string(&self.hits),
                tags_to_string(&self.removed_hits)
            )
            .unwrap();
        }

        &self.hits
    }
}

#[cfg(test)]
mod tests {
    use super::BarcodeBuilder;
    use crate::tags::Tag;
    use std::path::PathBuf;

    fn sink() -> PathBuf {
        PathBuf::from("/dev/null")
    }

    fn make_tag(seq: &str, max_mism: usize) -> Tag {
        Tag::new(seq, "info", max_mism)
    }

    // ── basic finding ─────────────────────────────────────────────────────────

    #[test]
    fn no_tags_found() {
        let tags = vec![make_tag("ACGT", 0)];
        let mut b = BarcodeBuilder::new(&tags, sink());
        assert!(b.get_barcode(b"TTTTTTTT").is_empty());
    }

    #[test]
    fn single_tag_found_at_correct_position() {
        let tags = vec![make_tag("ACGT", 0)];
        let mut b = BarcodeBuilder::new(&tags, sink());
        let hits = b.get_barcode(b"TTTTACGTTT");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].0, 4);
        assert_eq!(hits[0].1.get_seq(), "ACGT");
    }

    #[test]
    fn same_tag_found_twice() {
        let tags = vec![make_tag("ACGT", 0)];
        let mut b = BarcodeBuilder::new(&tags, sink());
        let hits = b.get_barcode(b"ACGTTTTTACGT");
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].0, 0);
        assert_eq!(hits[1].0, 8);
    }

    #[test]
    fn two_different_tags_found() {
        let tags = vec![make_tag("AAAA", 0), make_tag("CCCC", 0)];
        let mut b = BarcodeBuilder::new(&tags, sink());
        let hits = b.get_barcode(b"AAAATTTTCCCC");
        assert_eq!(hits.len(), 2);
        // results are sorted by position
        assert_eq!(hits[0].1.get_seq(), "AAAA");
        assert_eq!(hits[1].1.get_seq(), "CCCC");
    }

    // ── mismatch ──────────────────────────────────────────────────────────────

    #[test]
    fn tag_found_with_one_mismatch() {
        // ACGT vs ACTT — 1 mismatch, max_mism=1 → found
        let tags = vec![make_tag("ACGT", 1)];
        let mut b = BarcodeBuilder::new(&tags, sink());
        let hits = b.get_barcode(b"ACTT");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].0, 0);
    }

    // ── overlap resolution ────────────────────────────────────────────────────

    #[test]
    fn overlapping_tags_resolved() {
        // AAAAA @0..5 and AAAC @3..7 overlap; AAAAA has larger span → kept
        let tags = vec![make_tag("AAAAA", 0), make_tag("AAAC", 0)];
        let mut b = BarcodeBuilder::new(&tags, sink());
        let hits = b.get_barcode(b"AAAAAC");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].1.get_seq(), "AAAAA");
    }

    // ── reusability ───────────────────────────────────────────────────────────

    #[test]
    fn builder_reused_across_reads() {
        let tags = vec![make_tag("ACGT", 0)];
        let mut b = BarcodeBuilder::new(&tags, sink());
        assert!(b.get_barcode(b"TTTTTTTT").is_empty());
        let hits = b.get_barcode(b"ACGT");
        assert_eq!(hits.len(), 1);
        assert!(b.get_barcode(b"TTTTTTTT").is_empty());
    }
}
