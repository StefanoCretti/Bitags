import click
import polars as pl

from . import cli


@cli.command()
@click.argument("src", type=click.Path(exists=True))
@click.argument("out", type=click.Path(exists=False))
@click.option(
    "-c",
    "--column",
    required=True,
    help="Column to match against.",
)
@click.option(
    "-p",
    "--pattern",
    required=True,
    help="Value or pattern to match against the column.",
)
@click.option(
    "--regex",
    is_flag=True,
    default=False,
    help="Treat PATTERN as a regex instead of an exact match.",
)
@click.option(
    "-e",
    "--excluded",
    type=click.Path(exists=False),
    default=None,
    help="Optional path to write reads that do not pass the filter.",
)
def filter(
    src: str, out: str, column: str, pattern: str, regex: bool, excluded: str | None
) -> None:
    """Keep reads whose COLUMN matches PATTERN and write them to OUT.

    SRC is the input parquet file.
    OUT is the output parquet file path.
    """
    col = pl.col(column)
    expr = col.cast(pl.String).str.contains(pattern) if regex else col == pattern
    lf = pl.scan_parquet(src)
    lf.filter(expr).sink_parquet(out)

    # Creating two lazy sinks and collecting together does not seem to privide
    # a speedup but increase memory usage due to multiple buffers being open.
    # For now keeping as is, worth considering in case time issues arise.
    if excluded is not None:
        lf.filter(~expr).sink_parquet(excluded)
