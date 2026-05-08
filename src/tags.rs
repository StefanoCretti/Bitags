//! Tag objects for bitap pattern matching — ported from bitags.

const MAX_TAG_LEN: usize = 64;
const MAX_INFO_LEN: usize = 32;

/// Fixed-size bitap pattern container for the four DNA bases.
///
/// Each field holds a u64 bitmask where bit i is 0 if the corresponding base
/// appears at position i in the tag sequence, and 1 otherwise.
/// Fixed size keeps the struct on the stack and avoids hash-map overhead.
#[derive(Clone)]
pub struct BitapPatterns {
    a: u64,
    c: u64,
    t: u64,
    g: u64,
}

impl BitapPatterns {
    pub fn new(sequence: &str) -> Self {
        let mut new_self = Self { a: !0u64, c: !0u64, t: !0u64, g: !0u64 };
        for (i, b) in sequence.bytes().enumerate() {
            let mask = !(1u64 << i);
            match b {
                b'A' => new_self.a &= mask,
                b'C' => new_self.c &= mask,
                b'T' => new_self.t &= mask,
                b'G' => new_self.g &= mask,
                b'N' => {}
                _ => panic!("Invalid base: {b}"),
            }
        }
        new_self
    }

    pub fn get(&self, base: &u8) -> &u64 {
        match base {
            b'A' => &self.a,
            b'C' => &self.c,
            b'T' => &self.t,
            b'G' => &self.g,
            b'N' => &!0u64,
            _ => panic!("Invalid base: {base}"),
        }
    }
}

/// Tag with precomputed bitap patterns for fuzzy matching.
#[derive(Clone)]
pub struct Tag {
    seq_arr: [u8; MAX_TAG_LEN],
    info_arr: [u8; MAX_INFO_LEN],
    info_len: usize,
    pub len: usize,
    pub max_mism: usize,
    pub patterns: BitapPatterns,
}

impl Tag {
    pub fn new(seq: &str, info: &str, max_mism: usize) -> Self {
        let seq_len = seq.len();
        let mut seq_arr = [0u8; MAX_TAG_LEN];
        seq_arr[..seq_len].copy_from_slice(seq.as_bytes());

        let info_len = info.len();
        let mut info_arr = [0u8; MAX_INFO_LEN];
        info_arr[..info_len].copy_from_slice(info.as_bytes());

        Tag {
            seq_arr,
            info_arr,
            info_len,
            len: seq_len,
            max_mism,
            patterns: BitapPatterns::new(seq),
        }
    }

    pub fn get_seq(&self) -> &str {
        std::str::from_utf8(&self.seq_arr[..self.len]).unwrap()
    }

    pub fn get_info(&self) -> &str {
        std::str::from_utf8(&self.info_arr[..self.info_len]).unwrap()
    }
}
