use clap::Parser;

#[derive(Parser)]
struct Cli {
    raw_fastq: std::path::PathBuf,
    new_fastq: std::path::PathBuf,
    tags_file: std::path::PathBuf,
    buffer_size: usize,
}

fn main() {
    let args = Cli::parse();
    bitags::barcode_reads(
        args.raw_fastq,
        args.new_fastq,
        args.tags_file,
        args.buffer_size,
    );
}
