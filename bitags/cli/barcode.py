import click
import polars as pl

from bitags._barcoding import barcode_reads
from bitags._typing import ReadType

from . import cli


@cli.command()
@click.argument("src", type=click.Path(exists=True))
@click.argument("out", type=click.Path(exists=False))
@click.option(
    "-t",
    "--tags",
    required=True,
    type=click.Path(exists=True),
    help="""TSV file with tag definitions. Columns must be named `seq`, `info`
    and `max_mism`. `seq` is the tag nucleotide sequence. `info` is any
    metadata about the tag (name, type...). `max_mism` is the maximum number
    of allowed mismatches. Mismatches do NOT include indels.""",
)
@click.option(
    "-r",
    "--read",
    type=click.Choice(["r1", "r2"]),
    default=None,
    help="Read to barcode. If omitted, both r1 and r2 are processed.",
)
def barcode(src: str, out: str, tags: str, read: ReadType | None) -> None:
    """Use bitap algorithm to fuzzy match tags in the reads.

    SRC is the parquet file with the reads to barcode.
    OUT is the parquet file in which to save the results (must not exist).
    """
    read_types = (read,) if isinstance(read, str) else ("r1", "r2")
    tags_df = pl.read_csv(tags, separator="\t", comment_prefix="#")

    lf = pl.scan_parquet(src)
    for rt in read_types:
        lf = barcode_reads(lf, tags_df, read=rt)

    lf.sink_parquet(out)
