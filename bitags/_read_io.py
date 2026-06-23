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
) -> tuple[pl.LazyFrame, pl.LazyFrame]: ...
@overload
def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: None = ...,
    tags: Iterable[SamTag] | None = ...,
    lazy: Literal[True] = ...,
) -> pl.LazyFrame: ...
@overload
def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: str | None,
    tags: Iterable[SamTag] | None,
    lazy: Literal[False],
) -> None: ...


def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: str | None = None,
    tags: Iterable[SamTag] | None = None,
    lazy: bool = False,
) -> None | pl.LazyFrame | tuple[pl.LazyFrame, pl.LazyFrame]:
    """Write a LazyFrame with _r1/_r2 suffixed columns as gzipped fastq files.

    Always sinks r1. Sinks r2 only if a path is given and _r2 columns are present.

    tags is a list of (tag_name, tag_type, column_name) tuples appended to the R1
    description as tab-separated SAM optional fields. column_name must be the full
    column name as it appears in the frame (e.g. "tag_seq_r1" or a bare column).
    """

    def _tag_str(tag: SamTag) -> pl.Expr:
        """Convert a tag into its sam/bam tag string representation."""
        tag_name, tag_type, column = tag
        return pl.format(f"{tag_name}:{tag_type}:{{}}", pl.col(column))

    def _description_with_tags(col: str, tags: Iterable[SamTag]) -> pl.Expr:
        """Concatenate the bam/sam tags to the read description, id any."""
        return pl.concat_str(
            [pl.col(col).replace("", None)] + [_tag_str(tag) for tag in tags],
            separator="\t",
            ignore_nulls=True,
        ).fill_null("")

    def _prepare_read(
        lf: pl.LazyFrame,
        read: ReadType,
        tags: Iterable[SamTag],
    ) -> pl.LazyFrame:
        """Create a string column containing the info to sink to fastq."""

        cols = get_read_cols(read)
        description = _description_with_tags(cols.description, tags)
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
        raise ValueError("LazyFrame is paired but no r2 output path was provided.")
    if not paired and r2:
        raise ValueError("LazyFrame is unpaired but r2 output path was provided.")

    sink_kwargs = {
        "include_header": False,
        "quote_style": "never",
        "compression": "gzip",
        "lazy": True,
    }

    sinks = [_prepare_read(lf, "r1", tags or ()).sink_csv(r1, **sink_kwargs)]
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
