mod barcodes;
mod reads;
mod tags;

use barcodes::BarcodeBuilder;

use reads::ReadIterator;
use std::fs;
use std::io::{Result, Write};

fn write_read_buffer(read_buffer: &mut Vec<String>, out_file: &std::path::PathBuf) -> Result<()> {
    let file = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(out_file);

    match file {
        Ok(mut f) => {
            writeln!(f, "{}", read_buffer.join("\n"))?;
        }
        Err(e) => {
            eprintln!("Could not write to output file: {}", e);
        }
    }

    read_buffer.clear();

    Ok(())
}

pub fn barcode_reads(
    raw_fastq: std::path::PathBuf,
    new_fastq: std::path::PathBuf,
    tags_file: std::path::PathBuf,
    buffer_size: usize,
) {
    let reads_iter = ReadIterator::new(&raw_fastq);
    let mut res_buffer: Vec<String> = Vec::with_capacity(buffer_size);

    let name_pos: usize = 1;
    let seq_pos: usize = 2;
    let mism_pos: usize = 3;

    let tags_vec = tags::get_bitap_tags(&tags_file, name_pos, seq_pos, mism_pos);
    let mut builder = BarcodeBuilder::new(&tags_vec);
    for read in reads_iter {
        let hits = builder.get_barcode(&read.seq);
        res_buffer.push(read.to_fq_entry(hits));

        if res_buffer.len() == buffer_size {
            let _ = write_read_buffer(&mut res_buffer, &new_fastq);
        }
    }

    if res_buffer.len() > 0 {
        let _ = write_read_buffer(&mut res_buffer, &new_fastq);
    }
}
