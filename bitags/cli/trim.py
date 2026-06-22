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
@click.option(
    "--fill-empty",
    type=bool,
    default=True,
    show_default=True,
    help="Replace empty trimmed sequences with a single N/I placeholder.",
)
@click.option(
    "-m",
    "--min-length",
    type=int,
    default=None,
    help="Discard reads whose trimmed sequence is shorter than this length.",
)
@click.option(
    "-e",
    "--excluded",
    type=click.Path(exists=False),
    default=None,
    help="Optional path to write reads that do not pass the --min-length cutoff.",
)
def trim(
    src: str,
    out: str,
    regex: str,
    read: ReadType,
    tags_only: bool,
    fill_empty: bool,
    min_length: int | None,
    excluded: str | None,
) -> None:
    """Trim a read using a two-group regex.

    SRC is the parquet file with the barcoded reads.
    OUT is the parquet file in which to save the results (must not exist).
    """
    if excluded is not None and min_length is None:
        raise click.UsageError("-e/--excluded requires -m/--min-length to be set.")
    lf = pl.scan_parquet(src).pipe(
        trim_reads, regex=regex, read=read, tags_only=tags_only, fill_empty=fill_empty
    )
    if min_length is not None:
        length_expr = pl.col(f"sequence_{read}").str.len_chars() >= min_length
        lf.filter(length_expr).sink_parquet(out)
        if excluded is not None:
            lf.filter(~length_expr).sink_parquet(excluded)
    else:
        lf.sink_parquet(out)
