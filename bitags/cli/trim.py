import click
import polars as pl

from bitags._trimming import trim_reads
from bitags._typing import ReadType

from . import cli


@cli.command()
@click.argument("src", type=click.Path(exists=True))
@click.argument("out", type=click.Path(exists=False))
@click.option(
    "-r",
    "--regex",
    required=True,
    help="Two-group regex: group 1 removes 5' tags, group 2 removes 3' tags.",
)
@click.option(
    "--read",
    type=click.Choice(["r1", "r2"]),
    default="r1",
    show_default=True,
    help="Read to trim.",
)
@click.option(
    "--tags-only",
    is_flag=True,
    default=False,
    help="Trim only tag columns; leave sequence and quality untouched.",
)
def trim(src: str, out: str, regex: str, read: ReadType, tags_only: bool) -> None:
    """Trim a read using a two-group regex.

    SRC is the parquet file with the barcoded reads.
    OUT is the parquet file in which to save the results (must not exist).
    """
    pl.scan_parquet(src).pipe(
        trim_reads, regex=regex, read=read, tags_only=tags_only
    ).sink_parquet(out)
