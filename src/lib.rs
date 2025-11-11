//! Parallel barcoding of DNA/cDNA reads using specialized bitap fuzzy matching.

mod barcodes;
mod fastq_io;
mod tags;

use crossbeam::channel::{Receiver, Sender, bounded, unbounded};
use indicatif::{ProgressBar, ProgressStyle};
use noodles::fastq;

use itertools::Itertools;
use std::{fs, io, path, sync::Arc, thread};
use tempfile::{TempDir, tempdir};

use crate::barcodes::BarcodeBuilder;
use crate::fastq_io::{Writeable, concat_files, fastq_writer, read_batches};
use crate::tags::Tag;

const NUM_PAD_ZEROS: usize = 6;
// NOTE: Indirectly binds the number of chunks allowed, because if the number
// of chunks is 10^NUM_PAD_ZEROS or more, lexycographic order might not be
// preserved when joining chunks into a single output. Maybe add a check.

/// Load all tags present in a tags db file as bitap tags in a single vector.
pub fn get_bitap_tags<P: AsRef<path::Path>>(
    tags_file: P,
    info_pos: usize,
    seq_pos: usize,
    mism_pos: usize,
) -> Vec<Tag> {
    fs::read_to_string(tags_file)
        .unwrap()
        .lines()
        .filter(|x| !x.starts_with("#"))
        .map(|x| x.split("\t").collect::<Vec<&str>>())
        .map(|x| {
            Tag::new(
                x[seq_pos],
                x[info_pos],
                x[mism_pos].parse::<usize>().unwrap(),
            )
        })
        .collect::<Vec<Tag>>()
}

/// Process a batch of reads and save them to a dedicated output file.
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
                .map(|(pos, tag)| format!("{}:{}:{}", tag.get_seq(), tag.get_info(), pos))
                .join("|")
                .into_bytes(),
        );

        file_writer.write_record(&read).unwrap();
    }
    drop(file_writer);
}

/// Spawn a worker thread which barcodes batches of DNA/cDNA reads.
///
/// Worker threads are completely independent from each other (aside
/// from the progress bar) in order to facilitate parallelization.
/// Memory usage of the worker itself is negligible, hence the number
/// of workers is realistically bound only by the number of chunks which
/// can be loaded into memory at once.
fn spawn_worker(
    chunk_rx: Receiver<(usize, Vec<fastq::Record>)>,
    result_tx: Sender<(usize, Result<(), io::Error>)>,
    tags_vec: Vec<Tag>,
    out_dir: path::PathBuf,
    log_dir: path::PathBuf,
    progress_bar: Arc<ProgressBar>,
) -> thread::JoinHandle<()> {
    let tags_vec = tags_vec.clone(); // Almost zero cost and avoids lifetimes

    thread::spawn(move || {
        while let Ok((chunk_index, batch)) = chunk_rx.recv() {
            let chunk_name = format!("chunk_{:0NUM_PAD_ZEROS$}.fq.gz", chunk_index);
            let builder = BarcodeBuilder::new(&tags_vec, log_dir.join(&chunk_name));
            let file_writer = match fastq_writer(out_dir.join(&chunk_name)) {
                Ok(writer) => writer,
                Err(e) => {
                    eprintln!(
                        "Failed to create writer for chunk {:0NUM_PAD_ZEROS$}: {}",
                        chunk_index, e
                    );
                    result_tx.send((chunk_index, Err(e))).unwrap();
                    continue;
                }
            };

            barcode_batch(batch, builder, file_writer);

            progress_bar.set_message(format!("Done with chunk {}", chunk_index));
            progress_bar.inc(1);

            result_tx.send((chunk_index, Ok(()))).unwrap();
        }
    })
}

/// Spawn a reader thread that yields record batches from the input file
///
/// The reader creates batches of records from a fastq, but only up to the
/// number allowed by backpressure (i.e. the number of messages that can
/// be queued in the output channel). If the queue is full, the thread
/// is blocked and waits for new queue space to free up.
fn spawn_reader(
    fastq_file: path::PathBuf,
    batch_size: usize,
    chunk_tx: Sender<(usize, Vec<fastq::Record>)>,
    progress_bar: Arc<ProgressBar>,
) -> thread::JoinHandle<Result<usize, std::io::Error>> {
    // NOTE: Decompression is performed sequentially, without taking advantage
    // of the bgzip files. This is easier to implement and currently sufficient
    // to create batches faster than they are consumed. If IO ever becomes a
    // limiting factor, consider switching to async decompression.
    thread::spawn(move || -> io::Result<usize> {
        let mut total_chunks = 0;
        for (chunk_num, batch) in read_batches(fastq_file, batch_size).unwrap().enumerate() {
            match chunk_tx.send((chunk_num, batch.unwrap())) {
                Ok(_) => {
                    total_chunks += 1;
                    progress_bar.set_length(total_chunks as u64);
                }
                Err(_) => break,
            }
        }
        Ok(total_chunks)
    })
}

/// API to perform read barcoding using the bitap algorithm.
///
/// This function is the direct entry point for read barcoding.
/// If used directly you need to take care of input validation yourself.
/// There should be no reason to use this rather than the CLI tool.
pub fn barcode_reads<P: AsRef<path::Path>>(
    raw_fastq: P,
    new_fastq: P,
    bitap_tags: Vec<Tag>,
    num_workers: usize,
    batch_size: usize,
    overlap_log: P,
) -> io::Result<()> {
    // Tmp storage for the generated chunks
    let chunks_dir: TempDir = tempdir().unwrap();
    let logs_dir: TempDir = tempdir().unwrap();

    // Create channels for thread communication
    let (data_tx, data_rx) = bounded::<(usize, Vec<fastq::Record>)>(num_workers);
    let (done_tx, done_rx) = unbounded::<(usize, Result<(), io::Error>)>();

    // Create progress to pass to the
    let progress_bar = Arc::new(ProgressBar::new(0));
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{bar:40.cyan/blue}] Processed {pos}/{len} queued chunks ({msg})")
            .unwrap()
            .progress_chars("#>-"),
    );
    progress_bar.set_message("Starting to queue chunks.");

    // Create worker threads
    let mut worker_handles = Vec::<thread::JoinHandle<()>>::with_capacity(num_workers);
    for _ in 0..num_workers {
        let worker = spawn_worker(
            data_rx.clone(),
            done_tx.clone(),
            bitap_tags.clone(),
            chunks_dir.path().to_path_buf(),
            logs_dir.path().to_path_buf(),
            Arc::clone(&progress_bar),
        );
        worker_handles.push(worker);
    }
    drop(done_tx); // Need to remove all references to a channel to close it

    // Chunks are still read serially by the reader thread
    let reader_handle = spawn_reader(
        raw_fastq.as_ref().to_path_buf(),
        batch_size,
        data_tx,
        Arc::clone(&progress_bar),
    );
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

    // Finish progress bar
    progress_bar.finish_with_message("Processed all chunks.");

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
    concat_files(logs_dir, overlap_log)?;
    println!("Done merging chunks.");

    Ok(())
}
