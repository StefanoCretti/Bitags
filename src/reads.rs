//! Everything regarding reads, including I/O
//!
//! NOTE: At some point it will probably be replaced by oxbow

use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::tags::Tag;

pub struct Read {
    pub name: String,
    pub seq: String,
    pub quality: String,
}

impl Read {
    pub fn to_fq_entry(&self, hits: &Vec<(u8, Tag)>) -> String {
        let mut tags: Vec<String> = Vec::with_capacity(hits.len());
        for (pos, tag) in hits {
            tags.push(format!("[{}:{}]", pos, tag.get_name()));
        }

        format!(
            "{}::{}\n{}\n+\n{}",
            self.name,
            tags.join(""),
            self.seq,
            self.quality
        )
    }
}

pub struct ReadIterator {
    line_buffer: BufReader<File>,
}

impl ReadIterator {
    pub fn new(fastq_file: &std::path::PathBuf) -> Self {
        let file = File::open(fastq_file).expect("Could not open fastq file.");
        let line_buffer = BufReader::new(file);
        ReadIterator { line_buffer }
    }
}

impl Iterator for ReadIterator {
    type Item = Read;

    fn next(&mut self) -> Option<Self::Item> {
        let mut lines = [String::new(), String::new(), String::new(), String::new()];

        for i in 0..4 {
            lines[i].clear();
            match self.line_buffer.read_line(&mut lines[i]) {
                Ok(0) => {
                    // EOF
                    if i == 0 {
                        return None;
                    } else {
                        eprintln!("Warning: Incomplete FASTQ record at EOF");
                        return None;
                    }
                }
                Ok(_) => {
                    if lines[i].ends_with('\n') {
                        lines[i].pop();
                        if lines[i].ends_with('\r') {
                            lines[i].pop();
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Error reading line: {}", e);
                    return None;
                }
            }
        }

        let [name, seq, _, quality] = lines;
        return Some(Read { name, seq, quality });
    }
}
