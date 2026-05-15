import json

import click
import polars as pl

from bitags._classification import classify_reads
from bitags._typing import ReadType

from . import cli


@cli.command()
@click.argument("src", type=click.Path(exists=True))
@click.argument("out_dir", type=click.Path(exists=False))
@click.option(
    "-r",
    "--regex-json",
    required=True,
    type=click.Path(exists=True),
    help="JSON file mapping category labels to regexes.",
)
@click.option(
    "--read",
    type=click.Choice(["r1", "r2"]),
    default="r1",
    show_default=True,
    help="Read whose tag_type column is used for classification.",
)
def split(src: str, out_dir: str, regex_json: str, read: ReadType) -> None:
    """Classify reads by tag type and write per-category parquet partitions.

    SRC is the parquet file with the barcoded reads.
    OUT_DIR is the directory in which to write the hive-partitioned parquet files.
    """
    split_col = f"tag_type_{read}"
    with open(regex_json) as fh:
        regexes = json.load(fh)

    (
        pl.scan_parquet(src)
        .pipe(classify_reads, regexes, col=split_col, into="category")
        .sink_parquet(pl.PartitionBy(out_dir, key="category"), mkdir=True)
    )
