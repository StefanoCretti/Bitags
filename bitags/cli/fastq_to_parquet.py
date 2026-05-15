import click

from bitags._read_io import scan_fastq

from . import cli


@cli.command("fastq-to-parquet")
@click.argument("r1", type=click.Path(exists=True))
@click.argument("r2", type=click.Path(exists=True), required=False, default=None)
@click.option("-o", "--out", required=True, help="Output parquet file path.")
def fastq_to_parquet(r1: str, r2: str | None, out: str) -> None:
    """Convert single or paired fastq.gz files to a parquet file.

    R1 is the path to the R1 fastq.gz file.
    R2 is the optional path to the R2 fastq.gz file (paired reads only).
    """
    scan_fastq(r1, r2).sink_parquet(out)
