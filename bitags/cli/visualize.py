import json

import click
import polars as pl

from bitags._visualization import visualize_reads

from . import cli


@cli.command()
@click.argument("src", type=click.Path(exists=True))
@click.option("-n", default=10, show_default=True, help="Number of reads to display.")
@click.option(
    "--color-map",
    default=None,
    type=click.Path(exists=True),
    help="JSON file mapping tag type names to hex color strings.",
)
@click.option(
    "-o",
    "--out",
    default=None,
    help="Output file path (.html or .svg). If omitted, prints to terminal.",
)
def visualize(src: str, n: int, color_map: str | None, out: str | None) -> None:
    """Display reads with tag regions highlighted by type.

    SRC is the parquet file with the barcoded reads.
    """
    if color_map is not None:
        with open(color_map) as fh:
            cmap = json.load(fh)
    else:
        cmap = None

    pl.scan_parquet(src).pipe(visualize_reads, num_rows=n, color_map=cmap, out=out)
