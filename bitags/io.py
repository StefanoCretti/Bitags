import oxbow as ox
import polars as pl

from bitags._typing import ReadType, TagMerger, TagParser

# Default way to combine columns into a fastq read
_read_merge_expr: pl.Expr = pl.concat_str(
    [
        pl.format("@{} {}", pl.col("name"), pl.col("description")),
        pl.col("sequence"),
        pl.lit("+"),
        pl.col("quality"),
    ],
    separator="\n",
).alias("read")


def _split_tags(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Parse description into tag_seq, tag_type, and tag_pos columns.

    Expects description format: <description>::<tag_seq>|<tag_type>|<pos1>:<pos2>:...
    """
    return (
        lf.with_columns(pl.col("description").str.split_exact("::", 1).alias("_parsed"))
        .with_columns(
            pl.col("_parsed").struct.field("field_0").alias("description"),
            pl.col("_parsed")
            .struct.field("field_1")
            .str.split_exact("|", 2)
            .struct.rename_fields(["tag_seq", "tag_type", "tag_pos"])
            .alias("_tags"),
        )
        .drop("_parsed")
        .unnest("_tags")
        .with_columns(
            pl.when(pl.col("tag_pos").is_null() | (pl.col("tag_pos") == ""))
            .then(pl.lit([], dtype=pl.List(pl.String)))
            .otherwise(pl.col("tag_pos").str.split(":"))
            .cast(pl.List(pl.UInt32))
            .alias("tag_pos")
        )
    )


def _merge_tags(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Reconstruct the description column from tag_seq, tag_type, and tag_pos."""
    return (
        lf.with_columns(
            pl.concat_str(
                [
                    pl.col("tag_seq"),
                    pl.col("tag_type"),
                    pl.col("tag_pos").cast(pl.List(pl.String)).list.join(":"),
                ],
                separator="|",
            )
            .fill_null("")
            .alias("_tag_info")
        )
        .with_columns(
            pl.concat_str(
                [pl.col("description"), pl.col("_tag_info")], separator="::"
            ).alias("description")
        )
        .drop("_tag_info")
    )


def sink_fastq(lf: pl.LazyFrame, out: str, tag_merger: TagMerger | None = None) -> None:
    """Write a LazyFrame as a gzipped fastq file."""
    if tag_merger is not None:
        lf = tag_merger(lf)

    (
        lf.select(_read_merge_expr).sink_csv(
            out,
            include_header=False,
            quote_style="never",
            compression="gzip",
        )
    )


def scan_fastq(src: str, read: ReadType, tag_parser: TagParser | None) -> pl.LazyFrame:
    """Scan a fastq.gz file, apply tag parsing, and suffix all column names."""

    lf = ox.from_fastq(src).to_polars(lazy=True)
    assert isinstance(lf, pl.LazyFrame)

    if tag_parser is not None:
        lf = tag_parser(lf)

    return lf.rename({col: f"{col}_{read}" for col in lf.collect_schema().names()})


def fastq_to_parquet(
    r1: str,
    r2: str | None = None,
    *,
    out: str,
    tag_parser: TagParser | None = None,
) -> None:
    """Convert single or paired fastq.gz files to a single parquet file.

    All columns are suffixed with '_r1' (and '_r2' for paired reads).
    """
    tag_parser = tag_parser or _split_tags
    lf = scan_fastq(r1, "r1", tag_parser)
    if r2 is not None:
        lf = pl.concat([lf, scan_fastq(r2, "r2", tag_parser)], how="horizontal")
    lf.sink_parquet(out)


def parquet_to_fastq(
    src: str,
    r1: str,
    r2: str | None = None,
    *,
    tag_merger: TagMerger | None = None,
) -> None:
    """Convert a parquet file (created by fastq_to_parquet) back to fastq.gz.

    For paired reads, provide the r2 output path.
    """
    lf = pl.scan_parquet(src)
    schema_names = lf.collect_schema().names()
    merger = tag_merger or _merge_tags

    def _prepare(suffix: str) -> pl.LazyFrame:
        """Select and strip suffix from columns."""
        cols = [c for c in schema_names if c.endswith(f"_{suffix}")]
        return lf.select(cols).rename({c: c[: -(len(suffix) + 1)] for c in cols})

    sink_fastq(_prepare("r1"), r1, tag_merger=merger)
    if r2 is not None:
        sink_fastq(_prepare("r2"), r2, tag_merger=merger)
