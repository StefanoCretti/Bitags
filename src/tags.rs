//! Functionalities to create tag objects for bitap pattern matching.

const MAX_TAG_NAME_LEN: usize = 64;

/// Fixed size container for the bitap patterns of DNA bases.
///
/// Creates and stores the bitap pattern for each DNA base (A, C, T, G).
/// The bitap pattern for a letter is a u16 with all 1s except in the
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
/// - Due to the usage of u16, max pattern length allowed is 16 characters.
#[derive(Clone)]
pub struct BitapPatterns {
    a: u16,
    c: u16,
    t: u16,
    g: u16,
}

impl BitapPatterns {
    pub fn new(sequence: &str) -> Self {
        let mut new_self = Self {
            a: !0u16,
            c: !0u16,
            t: !0u16,
            g: !0u16,
        };
        for (i, b) in sequence.bytes().enumerate() {
            let mask = !(1u16 << i);
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
    pub fn get(&self, base: &u8) -> &u16 {
        match base {
            b'A' => &self.a,
            b'C' => &self.c,
            b'T' => &self.t,
            b'G' => &self.g,
            b'N' => &!0u16, // Will always be a mismatch
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
    name_arr: [u8; MAX_TAG_NAME_LEN],
    name_len: usize,
    pub len: usize,
    pub max_mism: usize,
    pub patterns: BitapPatterns,
}

impl Tag {
    pub fn new(name: &str, seq: &str, max_mism: usize) -> Self {
        let mut name_arr = [0u8; MAX_TAG_NAME_LEN];
        let name_byt = name.as_bytes();
        let name_len = name_byt.len().min(MAX_TAG_NAME_LEN);
        let len = seq.bytes().len();
        let patterns = BitapPatterns::new(&seq);

        name_arr[..name_len].copy_from_slice(&name_byt[..name_len]);

        return Tag {
            name_arr,
            name_len,
            len,
            max_mism,
            patterns,
        };
    }

    /// Return tag name in string form.
    pub fn get_name(&self) -> &str {
        std::str::from_utf8(&self.name_arr[..self.name_len]).unwrap()
    }
}
