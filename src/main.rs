//! CLI to simplify the usage of the tool, especially in a pipeline.

use bitags;
use clap::Parser;

#[derive(Parser)]
#[command(about = "A CLI app to find barcodes in DNA/cDNA reads using the bitap algorithm.")]
struct Cli {
    /// Path to the compressed fastq file to barcode.
    #[arg(short = 'i', long = "input")]
    raw_fastq: std::path::PathBuf,

    /// Path to the output compressed fastq file.
    #[arg(short = 'o', long = "output")]
    new_fastq: std::path::PathBuf,

    /// Path to the tags file (tsv with tags information).
    #[arg(short = 't', long = "tags")]
    tags_file: std::path::PathBuf,

    /// Index (0-based) of the column containing tag names.
    #[arg(short = 'n', long = "tags-name-pos", default_value = "1")]
    name_pos: usize,

    /// Index (0-based) of the column containing tag sequences.
    #[arg(short = 's', long = "tags-seq-pos", default_value = "2")]
    seq_pos: usize,

    /// Index (0-based) of the column containing the number of allowed mismatches per tag.
    #[arg(short = 'm', long = "tags-mism-pos", default_value = "3")]
    mism_pos: usize,

    /// Number of reads in each batch passed to a worker thread.
    #[arg(short = 'b', long, default_value = "500000")]
    batch_size: usize,

    /// Number of parallel worker threads.
    #[arg(short = 'w', long, default_value = "4")]
    num_workers: usize,

    /// Path to the tags overlap resolution log file.
    #[arg(short = 'l', long, default_value = "tag_overlap.log")]
    overlap_log: std::path::PathBuf,
}

fn main() -> std::io::Result<()> {
    let args = Cli::parse();
    let tags = bitags::get_bitap_tags(args.tags_file, args.name_pos, args.seq_pos, args.mism_pos);

    bitags::barcode_reads(
        args.raw_fastq,
        args.new_fastq,
        tags,
        args.num_workers,
        args.batch_size,
        args.overlap_log,
    )?;

    Ok(())
}
