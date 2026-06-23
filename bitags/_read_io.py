import warnings
from typing import Iterable, Literal, overload

import oxbow as ox
import polars as pl

from bitags._schema import get_read_cols, is_paired
from bitags._typing import ReadType, SamTag


@overload
def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: str,
    tags: Iterable[SamTag] | None = ...,
    lazy: Literal[True] = ...,
    keep_description: bool = ...,
) -> tuple[pl.LazyFrame, pl.LazyFrame]: ...
@overload
def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: None = ...,
    tags: Iterable[SamTag] | None = ...,
    lazy: Literal[True] = ...,
    keep_description: bool = ...,
) -> pl.LazyFrame: ...
@overload
def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: str | None,
    tags: Iterable[SamTag] | None,
    lazy: Literal[False],
    keep_description: bool = ...,
) -> None: ...


def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: str | None = None,
    tags: Iterable[SamTag] | None = None,
    lazy: bool = False,
    keep_description: bool = False,
) -> None | pl.LazyFrame | tuple[pl.LazyFrame, pl.LazyFrame]:
    """Write a LazyFrame with _r1/_r2 suffixed columns as gzipped fastq files.

    Always writes r1. Writes r2 only if r2 is provided; emits a warning instead
    of raising if the frame is paired but r2 is omitted.

    tags: SAM optional fields appended to the R1 description as tab-separated
    (tag_name, tag_type, column_name) tuples. column_name must match the full
    column name in the frame (e.g. "tag_seq_r1"). By default only the tags are
    written to the description; set keep_description=True to prepend the original
    description (note: this may interfere with alignment tools).
    """

    def _tag_str(tag: SamTag) -> pl.Expr:
        """Convert a tag into its sam/bam tag string representation."""
        tag_name, tag_type, column = tag
        return pl.format(f"{tag_name}:{tag_type}:{{}}", pl.col(column))

    def _description_with_tags(col: str, tags: Iterable[SamTag], keep: bool) -> pl.Expr:
        """Concatenate the bam/sam tags to the read description, if any."""
        tag_exprs = [_tag_str(tag) for tag in tags]
        parts = ([pl.col(col).replace("", None)] if keep else []) + tag_exprs
        if not parts:
            return pl.lit("")
        return pl.concat_str(parts, separator="\t", ignore_nulls=True).fill_null("")

    def _prepare_read(
        lf: pl.LazyFrame,
        read: ReadType,
        tags: Iterable[SamTag],
        keep: bool = True,
    ) -> pl.LazyFrame:
        """Create a string column containing the info to sink to fastq."""

        cols = get_read_cols(read)
        description = _description_with_tags(cols.description, tags, keep)
        return lf.select(
            pl.concat_str(
                [
                    pl.format("@{} {}", pl.col(cols.name), description),
                    pl.col(cols.seq),
                    pl.lit("+"),
                    pl.col(cols.qual),
                ],
                separator="\n",
            ).alias("read")
        )

    paired = is_paired(lf)
    if paired and r2 is None:
        warnings.warn(
            "Data is paired but no r2 path was provided; r2 will not be written.",
            UserWarning,
        )
    if not paired and r2:
        raise ValueError("LazyFrame is unpaired but r2 output path was provided.")

    sink_kwargs = {
        "include_header": False,
        "quote_style": "never",
        "compression": "gzip",
        "lazy": True,
    }

    sinks = [
        _prepare_read(lf, "r1", tags or (), keep=keep_description).sink_csv(
            r1, **sink_kwargs
        )
    ]
    if r2 is not None:
        sinks.append(_prepare_read(lf, "r2", ()).sink_csv(r2, **sink_kwargs))

    if lazy:
        return sinks[0] if len(sinks) == 1 else tuple(sinks)

    pl.collect_all(sinks)


def scan_fastq(r1: str, r2: str | None = None) -> pl.LazyFrame:
    """Scan a (paired) fastq.gz into a polars lazyframe."""

    def _scan_read(fq: str, read: ReadType) -> pl.LazyFrame:
        lf = ox.from_fastq(fq).to_polars(lazy=True)
        assert isinstance(lf, pl.LazyFrame)
        return lf.rename({col: f"{col}_{read}" for col in lf.collect_schema().names()})

    lf = _scan_read(r1, "r1")
    if r2 is not None:
        lf = pl.concat([lf, _scan_read(r2, "r2")], how="horizontal")

    return lf
