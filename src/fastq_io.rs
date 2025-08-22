//! Custom I/O functions for fastq files.

use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use noodles::fastq;

use std::io::{Read, Write}; // Inherited traits
use std::{fs, io, path};

const MB: usize = 1024 * 1024;
const BUFFER_CAPACITY: usize = 10 * MB;

//////////////////////////////////////////////////////////////////////////////
// Generic buffers for IO operations with and without compression ////////////
//////////////////////////////////////////////////////////////////////////////

/// Enumeration of readable sources to wrap in a buffer.
enum Readable {
    Plain(fs::File),
    Bgzip(MultiGzDecoder<fs::File>),
}

impl Read for Readable {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Readable::Plain(file) => file.read(buf),
            Readable::Bgzip(decoder) => decoder.read(buf),
        }
    }
}

/// Get a buffered reader from a file path (supports bgzip compression).
fn get_buff_reader<P: AsRef<path::Path>>(
    source: P,
    to_decompress: bool,
) -> io::Result<io::BufReader<Readable>> {
    let src = fs::File::open(source)?;
    let readable = match to_decompress {
        true => Readable::Bgzip(MultiGzDecoder::new(src)),
        false => Readable::Plain(src),
    };
    Ok(io::BufReader::with_capacity(BUFFER_CAPACITY, readable))
}

/// Enumeration of writeable sources to wrap in a buffer.
/// TODO: Remove the pub to make the API consistent.
pub enum Writeable {
    Plain(fs::File),
    Bgzip(GzEncoder<fs::File>),
}

impl Write for Writeable {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Writeable::Plain(file) => file.write(buf),
            Writeable::Bgzip(encoder) => encoder.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Writeable::Plain(file) => file.flush(),
            Writeable::Bgzip(encoder) => encoder.flush(),
        }
    }
}

/// Get a buffered writer to a file path (supports bgzip compression).
///
/// NOTE: Currently using default compression level. This is a bit
/// time consuming but it halves file size with respect to the fast
/// compression level.
fn get_buff_writer<P: AsRef<path::Path>>(
    sink: P,
    to_compress: bool,
) -> io::Result<io::BufWriter<Writeable>> {
    let out = fs::File::create(sink)?;
    let writeable = match to_compress {
        true => Writeable::Bgzip(GzEncoder::new(out, Compression::default())),
        false => Writeable::Plain(out),
    };
    Ok(io::BufWriter::with_capacity(BUFFER_CAPACITY, writeable))
}

//////////////////////////////////////////////////////////////////////////////
// Functions to read noodle records in batches ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// NOTE: Might replace with oxbow, but currently not enough documentation
// NOTE: Currently readers and writers are sequential, might need to use async
//       noodles if IO becomes a significant bottleneck of the pipeline.

/// Iterator of batches of noodle records for parallel processing dispatch.
pub struct ReadBatches {
    fastq_reader: fastq::io::Reader<io::BufReader<Readable>>,
    batch_size: usize,
}

impl ReadBatches {
    /// Get the next batch of records from the source.
    fn next_batch(&mut self) -> std::io::Result<Option<Vec<fastq::Record>>> {
        let mut batch = Vec::with_capacity(self.batch_size);

        for _ in 0..self.batch_size {
            let mut record = fastq::Record::default();
            match self.fastq_reader.read_record(&mut record) {
                Ok(0) => break, // EOF
                Ok(_) => batch.push(record),
                Err(e) => return Err(e),
            }
        }

        if batch.is_empty() {
            Ok(None)
        } else {
            Ok(Some(batch))
        }
    }
}

impl Iterator for ReadBatches {
    type Item = std::io::Result<Vec<fastq::Record>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_batch() {
            Ok(Some(batch)) => Some(Ok(batch)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Create and return an iterator of noodle record batches of the desired size.
pub fn read_batches<P: AsRef<path::Path>>(
    fastq_file: P,
    batch_size: usize,
) -> io::Result<ReadBatches> {
    let input_buffer = get_buff_reader(&fastq_file, true)?;
    let fastq_reader = fastq::io::Reader::new(input_buffer);
    Ok(ReadBatches {
        fastq_reader,
        batch_size,
    })
}

/// Return a handle to write processed reads to a gzipped fastq file
/// TODO: Change to hide writeable and make the API consistent.
pub fn fastq_writer<P: AsRef<path::Path>>(
    fastq_file: P,
) -> io::Result<fastq::io::Writer<io::BufWriter<Writeable>>> {
    let output_buffer = get_buff_writer(fastq_file, true)?;
    let fastq_writer = fastq::io::Writer::new(output_buffer);
    Ok(fastq_writer)
}

//////////////////////////////////////////////////////////////////////////////
// Utility functions /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// Concatenate the content of all files in a directory into a single one.
///
/// Source files are merged in lexicographic order.
/// Bytes are copied as is, meaning that compression does not affect the result.
///
/// NOTE: Currently this function is not very efficient in terms of disk usage.
/// If we denote with _n_ the size of the original unchunked input, currently
/// disk usage is _3n_ (original input, chunks from it, merged output). This
/// could be easily changed to _2n_ by removing the chunk after merging it into
/// the output, but manual file cleanup from a tmp dir feels very hacky. Since
/// the function will probably become _n_ at some point by implementing a
/// consumer model or something of the sort using async, left as _3n_ for now
/// as disk space is not likely to be a bottleneck anyway.
pub fn concat_files<P1: AsRef<path::Path>, P2: AsRef<path::Path>>(
    chunks_dir: P1,
    merged_file: P2,
) -> io::Result<()> {
    let mut out_writer = get_buff_writer(merged_file, false)?;

    let mut chunks: Vec<path::PathBuf> = fs::read_dir(&chunks_dir)?
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|path| path.is_file())
        .collect();
    chunks.sort();

    for chunk in chunks {
        let mut chunk_reader = get_buff_reader(chunk, false).unwrap();
        io::copy(&mut chunk_reader, &mut out_writer).unwrap();
    }

    out_writer.flush()?;

    Ok(())
}
