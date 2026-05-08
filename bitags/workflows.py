import json

import polars as pl

from bitags._typing import ReadType
from bitags.barcoding import barcode_reads
from bitags.classification import classify_reads
from bitags.manipulation import embed_barcode, to_unpaired
from bitags.read_io import scan_fastq, sink_fastq
from bitags.trimming import trim_reads
from bitags.visualization import render_reads

__all__ = [
    "barcode",
    "embed",
    "fastq_to_parquet",
    "parquet_to_fastq",
    "split",
    "trim",
    "visualize",
]


def fastq_to_parquet(r1: str, r2: str | None = None, *, out: str) -> None:
    """Convert single or paired fastq.gz files to a parquet file.

    Columns are suffixed with '_r1' (and '_r2' for paired reads).
    The operation is streamed without collecting in memory.
    """
    scan_fastq(r1, r2).sink_parquet(out)


def parquet_to_fastq(src: str, *, r1: str, r2: str | None = None) -> None:
    """Convert a parquet file back to fastq.gz files.

    Expects columns suffixed with '_r1' (and '_r2' for paired reads).
    Raises if the parquet contains _r2 columns but no r2 path is provided.
    """
    sink_fastq(pl.scan_parquet(src), r1=r1, r2=r2)


def barcode(src: str, out: str, *, tags: str, read: ReadType | None = None) -> None:
    """Apply bitap tag matching to reads in a parquet file and write the result.

    Loads tags from a TSV file (columns: seq, info, max_mism) and runs fuzzy
    matching against the sequence columns. If read is given, only that read is
    barcoded; otherwise both r1 and r2 are processed.
    """

    read_types = (read,) if isinstance(read, str) else ("r1", "r2")
    tags_df = pl.read_csv(tags, separator="\t", comment_prefix="#")

    lf = pl.scan_parquet(src)
    for rt in read_types:
        lf = barcode_reads(lf, tags_df, read=rt)

    lf.sink_parquet(out)


def split(src: str, out_dir: str, *, regex_json: str, read: ReadType = "r1") -> None:
    """Classify reads by R2 barcode and write per-category parquet partitions.

    Reads are classified as DNA, RNA, or Other based on the R2 tag_type column.
    Output: hive-partitioned parquet under out_dir (category=DNA/, category=RNA/, ...).
    """

    split_col = f"tag_type_{read}"
    with open(src, "r") as handle:
        regexes = json.load(handle)

    (
        pl.scan_parquet(src)
        .pipe(classify_reads, regexes, col=split_col, into="category")
        .sink_parquet(pl.PartitionBy(out_dir, key=split_col), mkdir=True)
    )


def trim(
    src: str,
    out: str,
    *,
    regex: str,
    read: ReadType = "r1",
    tags_only: bool = False,
) -> None:
    """Trim a read in a parquet file using a two-group regex and write the result.

    The regex is matched against tag_type; group 1 defines tags to remove from the
    5' end, group 2 from the 3' end. Sequence and quality are trimmed alongside
    the tag columns and tag positions are adjusted to the new sequence start.
    When tags_only=True, only the tag columns are trimmed; sequence and quality are untouched.
    """
    pl.scan_parquet(src).pipe(
        trim_reads, regex=regex, read=read, tags_only=tags_only
    ).sink_parquet(out)


def embed(
    src: str,
    out: str,
    *,
    source_read: ReadType,
    target_read: ReadType,
    separator: str = "::",
    remove_source: bool = True,
) -> None:
    """Embed the barcode from source_read into the read name of target_read.

    The barcode (tag_seq_{source_read}) is appended to name_{target_read} using
    separator. If remove_source=True, all source_read columns are dropped; if the
    source was r1, the remaining r2 columns are renamed to r1 so the output is
    treated as unpaired.
    """
    lf = embed_barcode(pl.scan_parquet(src), source_read, target_read, separator)
    if remove_source:
        lf = to_unpaired(lf, source_read)
    lf.sink_parquet(out)


def visualize(
    src: str,
    *,
    n: int = 10,
    color_map: str | None = None,
    out: str | None = None,
) -> None:
    """Display n reads from src with tag regions highlighted by type.

    If color_map is provided, it must be a path to a JSON file mapping tag type
    names to hex color strings (e.g. {"Adap_l": "#555555"}). If omitted, a
    divergent palette is generated automatically from the tag types in the data.
    If out is provided, saves to file instead of printing to terminal.
    Supported formats: .html, .svg (inferred from extension).
    """
    with open(src, "r") as handle:
        cmap = json.load(handle)
    pl.scan_parquet(src).pipe(render_reads, num_rows=n, color_map=cmap, out=out)
