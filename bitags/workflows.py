import io

import polars as pl
from rich.console import Console

from bitags._typing import MoleculeType, ReadType
from bitags.classification import classify_reads, split
from bitags.patterns import r1_trim_regex, r2_barcode_regex, r2_regex, r2_trim_regex
from bitags.read_io import scan_fastq, sink_fastq
from bitags.barcoding import barcode_reads
from bitags.trimming import extract_barcode, trim_reads
from bitags.viz import render_read_pair


def split_by_type(src: str, *, schema: str, out_dir: str) -> None:
    """Classify reads by R2 barcode and write per-category parquet partitions.

    Reads are classified as DNA, RNA, or Other based on the R2 tag_type column.
    Output: hive-partitioned parquet under out_dir (category=DNA/, category=RNA/, ...).
    """
    regexes = {
        "DNA": r2_regex(schema, "dna"),
        "RNA": r2_regex(schema, "rna"),
    }

    (
        pl.scan_parquet(src)
        .pipe(classify_reads, regexes, col="tag_type_r2", into="category")
        .pipe(split, on="category", out_dir=out_dir)
    )


def trim(
    src: str,
    *,
    schema: str,
    molecule: MoleculeType,
    out: str,
    reads: tuple[ReadType, ...] = ("r1", "r2"),
) -> None:
    """Trim barcodes from the specified reads of a single molecule type.

    Expects input pre-filtered to one molecule type (e.g. output of classify_and_split).
    """
    _regex = {"r1": r1_trim_regex, "r2": r2_trim_regex}
    lf = pl.scan_parquet(src)
    for read in reads:
        lf = lf.pipe(trim_reads, _regex[read](schema, molecule), read=read)
    lf.sink_parquet(out)


_BARCODE_REGEX = {"r2": r2_barcode_regex}


def transfer_barcode(
    src: str,
    src_read: ReadType,
    dst_read: ReadType,
    *,
    schema: str,
    molecule: MoleculeType,
    out: str,
) -> None:
    """Extract the barcode from src_read and store it in dst_read's tag_seq.

    Writes only dst_read columns. tag_type and tag_pos for dst_read are cleared
    since positional info no longer applies after transfer.
    """
    lf = pl.scan_parquet(src)
    dst_cols = [c for c in lf.collect_schema().names() if c.endswith(f"_{dst_read}")]

    (
        lf.pipe(
            extract_barcode,
            _BARCODE_REGEX[src_read](schema, molecule),
            read=src_read,
            into="barcode",
        )
        .with_columns(
            pl.col("barcode").alias(f"tag_seq_{dst_read}"),
            pl.lit("").alias(f"tag_type_{dst_read}"),
            pl.lit([], dtype=pl.List(pl.UInt32)).alias(f"tag_pos_{dst_read}"),
        )
        .drop("barcode")
        .select(dst_cols)
        .sink_parquet(out)
    )


def visualize(
    src: str,
    *,
    n: int = 10,
    color_map: dict[str, str] | None = None,
    out: str | None = None,
) -> None:
    """Display n reads from src with tag regions highlighted by type.

    If out is provided, saves to file instead of printing to terminal.
    Supported formats: .html, .svg (inferred from extension).
    """

    console = Console(record=True, file=io.StringIO()) if out else Console()

    rows = pl.scan_parquet(src).head(n).collect().to_dicts()
    for row in rows:
        render_read_pair(row, console=console, color_map=color_map)
        console.rule()

    if out:
        ext = out.rsplit(".", 1)[-1].lower()
        if ext == "html":
            content = console.export_html()
        elif ext == "svg":
            content = console.export_svg()
        else:
            raise ValueError(f"Unsupported format '.{ext}'. Use .html or .svg")
        with open(out, "w") as f:
            f.write(content)


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
