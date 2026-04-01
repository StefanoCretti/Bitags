//! Single read barcoding functionalities (batching, bitap matching, overlap resolution)

mod builder;
mod matching;
mod overlap;

pub use builder::BarcodeBuilder;

use crate::tags::Tag;

const MAX_TAG_HITS: usize = 64; // Maximum number of tags aligned per read
const MAX_TAG_MISM: usize = 2; // Maximum number of mismatches per tag

type AlignedTags<'a> = Vec<(usize, &'a Tag)>;
