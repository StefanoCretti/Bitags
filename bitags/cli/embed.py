import click
import polars as pl

from bitags._manipulation import embed_barcode, to_unpaired
from bitags._typing import ReadType

from . import cli


@cli.command()
@click.argument("src", type=click.Path(exists=True))
@click.argument("out", type=click.Path(exists=False))
@click.option(
    "--source-read",
    required=True,
    type=click.Choice(["r1", "r2"]),
    help="Read whose barcode is embedded into the target read name.",
)
@click.option(
    "--target-read",
    required=True,
    type=click.Choice(["r1", "r2"]),
    help="Read whose name receives the embedded barcode.",
)
@click.option(
    "--separator",
    default="::",
    show_default=True,
    help="String separating the original name and the barcode.",
)
@click.option(
    "--no-remove-source",
    is_flag=True,
    default=False,
    help="Keep source read columns instead of dropping them.",
)
def embed(
    src: str,
    out: str,
    source_read: ReadType,
    target_read: ReadType,
    separator: str,
    no_remove_source: bool,
) -> None:
    """Embed the barcode from source-read into the read name of target-read.

    SRC is the parquet file with the barcoded reads.
    OUT is the parquet file in which to save the results (must not exist).
    """
    lf = embed_barcode(pl.scan_parquet(src), source_read, target_read, separator)
    if not no_remove_source:
        lf = to_unpaired(lf, source_read)
    lf.sink_parquet(out)
