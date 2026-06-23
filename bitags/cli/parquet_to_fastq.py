import click
import polars as pl

from bitags._read_io import sink_fastq
from bitags._typing import SamTag

from . import cli


@cli.command("parquet-to-fastq")
@click.argument("src", type=click.Path(exists=True))
@click.option("--r1", required=True, help="Output path for the R1 fastq.gz file.")
@click.option(
    "--r2",
    default=None,
    help="Output path for the R2 fastq.gz file (paired reads only).",
)
@click.option(
    "-t",
    "--tag",
    "tags",
    nargs=3,
    multiple=True,
    metavar="NAME TYPE COLUMN",
    help="Append a column as a SAM tag to the R1 description. Repeatable.",
)
@click.option(
    "--keep-description",
    is_flag=True,
    default=False,
    help="Prepend the original description to the SAM tags (may interfere with alignment tools).",
)
def parquet_to_fastq(
    src: str, r1: str, r2: str | None, tags: tuple[SamTag, ...], keep_description: bool
) -> None:
    """Convert a parquet file back to fastq.gz files.

    SRC is the parquet file with the reads to convert.
    """
    sink_fastq(
        pl.scan_parquet(src),
        r1=r1,
        r2=r2,
        tags=tags or None,
        keep_description=keep_description,
    )
