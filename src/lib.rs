mod barcodes;
mod file_ops;
mod tags;

use barcodes::BarcodeBuilder;
use file_ops::{Writeable, concat_files, fastq_writer, read_batches};
use noodles::fastq;
use std::{io, path};
use tempfile::{TempDir, tempdir};

fn barcode_batch(
    read_batch: Vec<fastq::Record>,
    mut builder: BarcodeBuilder,
    mut file_writer: fastq::io::Writer<io::BufWriter<Writeable>>,
) {
    for mut read in read_batch.into_iter() {
        let hits = builder.get_barcode(read.sequence());

        read.description_mut().extend_from_slice(b"::");
        read.description_mut().extend(
            hits.into_iter()
                .map(|(pos, tag)| format!("[{}:{}]", pos, tag.get_name()))
                .collect::<Vec<String>>()
                .join("")
                .into_bytes(),
        );

        file_writer.write_record(&read).unwrap();
    }
    drop(file_writer);
}

pub fn barcode_reads<P: AsRef<path::Path>>(
    raw_fastq: P,
    new_fastq: P,
    tags_file: P,
) -> io::Result<()> {
    let name_pos: usize = 1;
    let seq_pos: usize = 2;
    let mism_pos: usize = 3;
    let batch_size: usize = 500_000;

    let chunks_dir: TempDir = tempdir().unwrap();
    let tags_vec = tags::get_bitap_tags(&tags_file, name_pos, seq_pos, mism_pos);
    for (j, batch) in read_batches(&raw_fastq, batch_size).unwrap().enumerate() {
        println!("Starting to process chunk {:05}", j);

        let builder = BarcodeBuilder::new(&tags_vec);
        let chunk_name = format!("chunk_{:05}.fq.gz", j);
        let file_writer = fastq_writer(&chunks_dir.path().join(chunk_name))?;
        barcode_batch(batch.unwrap(), builder, file_writer);

        println!("Done processing chunk {:05}", j);
    }

    println!("Merging chunks...");
    concat_files(chunks_dir, new_fastq)?;
    println!("Done merging chunks!");

    Ok(())
}
