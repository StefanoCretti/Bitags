//! Functionalities to create tag objects for bitap pattern matching.

const MAX_TAG_LEN: usize = 32;
const MAX_INFO_LEN: usize = 32;

/// Fixed size container for the bitap patterns of DNA bases.
///
/// Creates and stores the bitap pattern for each DNA base (A, C, T, G).
/// The bitap pattern for a letter is a u32 with all 1s except in the
/// positions where the letter is found in the sequence.
///
/// As an example, the sequence "ACAG" will yield the following patterns:
/// - "A" -> 0b1111111111110101
/// - "C" -> 0b1111111111111011
/// - "G" -> 0b1111111111111110
///
/// This is definitely not the most elegant or flexible approach but:
/// - Hashmaps are too slow due to collision checking
/// - Positional arrays are not very readable plus the characters of the
///   bases are not contiguous in ASCII, which is annoying for indexing.
/// - This has fixed size so it is on the stack and thus can be copied/cloned.
/// - Allows to handle corner cases (e.g. N) in the same place as all
///   the other pattern logic.
///
/// Notes:
/// - Due to the usage of u32, max pattern length allowed is 32 characters.
#[derive(Clone)]
pub struct BitapPatterns {
    a: u32,
    c: u32,
    t: u32,
    g: u32,
}

impl BitapPatterns {
    pub fn new(sequence: &str) -> Self {
        let mut new_self = Self {
            a: !0u32,
            c: !0u32,
            t: !0u32,
            g: !0u32,
        };
        for (i, b) in sequence.bytes().enumerate() {
            let mask = !(1u32 << i);
            match b {
                b'A' => new_self.a &= mask,
                b'C' => new_self.c &= mask,
                b'T' => new_self.t &= mask,
                b'G' => new_self.g &= mask,
                b'N' => {} // Will always be a mismatch but allowed
                _ => panic!("Invalid base: {b}"),
            }
        }
        return new_self;
    }
    pub fn get(&self, base: &u8) -> &u32 {
        match base {
            b'A' => &self.a,
            b'C' => &self.c,
            b'T' => &self.t,
            b'G' => &self.g,
            b'N' => &!0u32, // Will always be a mismatch
            _ => panic!("Invalid base: {base}"),
        }
    }
}

/// Searcheable tag via bitap algorithm
///
/// Tag which was obtained from a tag db file. Stores the bitap patterns
/// and the maximum number of allowed mismatches for alignment purposes.
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
        let seq_len = seq.bytes().len();
        let mut seq_arr = [0u8; MAX_TAG_LEN];
        seq_arr[..seq_len].copy_from_slice(&seq.as_bytes()[..seq_len]);

        let info_len = info.bytes().len();
        let mut info_arr = [0u8; MAX_INFO_LEN];
        info_arr[..info_len].copy_from_slice(&info.as_bytes()[..info_len]);

        let patterns = BitapPatterns::new(&seq);

        return Tag {
            seq_arr,
            info_arr,
            info_len,
            len: seq_len,
            max_mism,
            patterns,
        };
    }

    /// Return tag sequence in string form.
    pub fn get_seq(&self) -> &str {
        std::str::from_utf8(&self.seq_arr[..self.len]).unwrap()
    }

    /// Return tag info in string form
    pub fn get_info(&self) -> &str {
        std::str::from_utf8(&self.info_arr[..self.info_len]).unwrap()
    }
}
