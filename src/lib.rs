mod barcodes;
mod file_ops;
mod tags;

use barcodes::BarcodeBuilder;
use crossbeam::channel::{Receiver, Sender, bounded, unbounded};
use file_ops::{Writeable, concat_files, fastq_writer, read_batches};
use noodles::fastq;
use std::{io, path, thread};
use tempfile::{TempDir, tempdir};

use crate::tags::Tag;

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

fn spawn_worker(
    chunk_rx: Receiver<(usize, Vec<fastq::Record>)>,
    result_tx: Sender<(usize, Result<(), io::Error>)>,
    tags_vec: Vec<Tag>,
    out_dir: path::PathBuf,
) -> thread::JoinHandle<()> {
    let tags_vec = tags_vec.clone();

    thread::spawn(move || {
        while let Ok((chunk_index, batch)) = chunk_rx.recv() {
            println!("Starting to process chunk {:05}", chunk_index);

            let builder = BarcodeBuilder::new(&tags_vec);
            let chunk_name = format!("chunk_{:05}.fq.gz", chunk_index);
            let file_writer = match fastq_writer(out_dir.join(chunk_name)) {
                Ok(writer) => writer,
                Err(e) => {
                    eprintln!("Failed to create writer for chunk {}: {}", chunk_index, e);
                    result_tx.send((chunk_index, Err(e))).unwrap();
                    continue;
                }
            };

            barcode_batch(batch, builder, file_writer);
            println!("Done processing chunk {:05}", chunk_index);

            result_tx.send((chunk_index, Ok(()))).unwrap();
        }
    })
}

fn spawn_reader(
    fastq_file: path::PathBuf,
    batch_size: usize,
    chunk_tx: Sender<(usize, Vec<fastq::Record>)>,
) -> thread::JoinHandle<Result<usize, std::io::Error>> {
    thread::spawn(move || -> io::Result<usize> {
        let mut total_chunks = 0;
        for (chunk_num, batch) in read_batches(fastq_file, batch_size).unwrap().enumerate() {
            match chunk_tx.send((chunk_num, batch.unwrap())) {
                Ok(_) => total_chunks += 1,
                Err(_) => break,
            }
        }
        Ok(total_chunks)
    })
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
    let num_workers: usize = 8;

    let chunks_dir: TempDir = tempdir().unwrap();
    let tags_vec = tags::get_bitap_tags(&tags_file, name_pos, seq_pos, mism_pos);

    // Create channels for thread communication
    let (data_tx, data_rx) = bounded::<(usize, Vec<fastq::Record>)>(num_workers);
    let (done_tx, done_rx) = unbounded::<(usize, Result<(), io::Error>)>();

    // Create worker threads
    let mut worker_handles = Vec::<thread::JoinHandle<()>>::with_capacity(num_workers);
    for _ in 0..num_workers {
        let worker = spawn_worker(
            data_rx.clone(),
            done_tx.clone(),
            tags_vec.clone(),
            chunks_dir.path().to_path_buf(),
        );
        worker_handles.push(worker);
    }
    drop(done_tx); // Need to remove all references to a channel to close it

    // Chunks are still read serially by the reader thread
    let reader_handle = spawn_reader(raw_fastq.as_ref().to_path_buf(), batch_size, data_tx);
    let total_chunks = reader_handle.join().unwrap()?;

    // Collect completion results
    let mut completed_chunks = 0;
    let mut errors = Vec::new();

    while completed_chunks < total_chunks {
        match done_rx.recv() {
            Ok((chunk_index, result)) => {
                completed_chunks += 1;
                if let Err(e) = result {
                    errors.push((chunk_index, e));
                }
            }
            Err(_) => break, // All workers finished and dropped their senders
        }
    }

    // Wait for all workers to finish
    for handle in worker_handles {
        handle.join().unwrap();
    }

    // Check for errors
    if !errors.is_empty() {
        eprintln!("Errors occurred in {} chunks:", errors.len());
        for (chunk_index, error) in errors {
            eprintln!("  Chunk {}: {}", chunk_index, error);
        }
        return Err(io::Error::new(
            io::ErrorKind::Other,
            "Some chunks failed to process",
        ));
    }

    println!("Merging chunks...");
    concat_files(chunks_dir, new_fastq)?;
    println!("Done merging chunks!");

    Ok(())
}
