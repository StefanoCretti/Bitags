import click


@click.group()
def cli():
    """Bitap-based barcoding and downstream processing of sequencing reads."""


from . import (  # noqa: E402, F401
    barcode,
    classify,
    embed,
    fastq_to_parquet,
    filter,
    parquet_to_fastq,
    trim,
    visualize,
)
