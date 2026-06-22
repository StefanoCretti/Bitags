import json

import click
import polars as pl

from bitags._classification import classify_reads
from bitags._typing import ReadType

from . import cli


@cli.command()
@click.argument("src", type=click.Path(exists=True))
@click.argument("out", type=click.Path(exists=False))
@click.option(
    "-j",
    "--regex-json",
    default=None,
    type=click.Path(exists=True),
    help="JSON file mapping category labels to regexes.",
)
@click.option(
    "-c",
    "--category",
    "categories",
    multiple=True,
    nargs=2,
    metavar="LABEL REGEX",
    help="Category label and regex. Repeatable; order is preserved.",
)
@click.option(
    "-r",
    "--read",
    type=click.Choice(["r1", "r2"]),
    default="r1",
    show_default=True,
    help="Read whose tag_type column is used for classification.",
)
@click.option(
    "-o",
    "--other",
    default="Other",
    show_default=True,
    help="Label assigned to reads matching no category.",
)
def classify(
    src: str,
    out: str,
    regex_json: str | None,
    categories: tuple[str, ...],
    read: ReadType,
    other: str,
) -> None:
    """Classify reads by tag type and write the result to a parquet file.

    SRC is the parquet file with the barcoded reads.
    OUT is the output parquet file path.

    Provide regexes either as a JSON file (-j) or as inline pairs (-c):

        bitags classify src.pq out.pq -j regexes.json
        bitags classify src.pq out.pq -c DNA <pattern> -c RNA <pattern>
    """
    if regex_json and categories:
        raise click.UsageError("Use either -j/--regex-json or -c/--category, not both.")
    if not regex_json and not categories:
        raise click.UsageError("Provide regexes via -j/--regex-json or -c/--category.")

    if regex_json:
        with open(regex_json) as fh:
            regexes: dict[str, str] = json.load(fh)
    else:
        regexes = {label: pattern for label, pattern in categories}

    split_col = f"tag_type_{read}"
    (
        pl.scan_parquet(src)
        .pipe(classify_reads, regexes, col=split_col, into="category", other=other)
        .sink_parquet(out)
    )
