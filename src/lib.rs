mod barcodes;
mod tags;

use barcodes::BarcodeBuilder;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use noodles::fastq;

use std::{fs::File, io::BufReader, io::BufWriter};

pub fn barcode_reads(
    raw_fastq: std::path::PathBuf,
    new_fastq: std::path::PathBuf,
    tags_file: std::path::PathBuf,
) {
    let mut reads_reader = File::open(raw_fastq)
        .map(MultiGzDecoder::new)
        .map(BufReader::new)
        .map(fastq::io::Reader::new)
        .unwrap();

    let mut reads_writer = File::create(new_fastq)
        .map(|f| GzEncoder::new(f, Compression::default()))
        .map(BufWriter::new)
        .map(fastq::io::Writer::new)
        .unwrap();

    let name_pos: usize = 1;
    let seq_pos: usize = 2;
    let mism_pos: usize = 3;

    let tags_vec = tags::get_bitap_tags(&tags_file, name_pos, seq_pos, mism_pos);
    let mut builder = BarcodeBuilder::new(&tags_vec);
    for (i, res) in reads_reader.records().enumerate() {
        if i % 1_000_000 == 0 {
            println!("Processing read {}", i);
        }

        let mut read = res.unwrap();
        let hits = builder.get_barcode(read.sequence());

        read.description_mut().extend_from_slice(b"::");
        read.description_mut().extend(
            hits.into_iter()
                .map(|(pos, tag)| format!("[{}:{}]", pos, tag.get_name()))
                .collect::<Vec<String>>()
                .join("")
                .into_bytes(),
        );

        reads_writer.write_record(&read).unwrap();
    }

    drop(reads_writer);
}
