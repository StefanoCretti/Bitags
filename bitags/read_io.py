from typing import Literal, overload

import oxbow as ox
import polars as pl

from bitags._typing import ReadType


@overload
def sink_fastq(
    lf: pl.LazyFrame, *, r1: str, r2: str, lazy: Literal[True] = ...
) -> tuple[pl.LazyFrame, pl.LazyFrame]: ...
@overload
def sink_fastq(
    lf: pl.LazyFrame, *, r1: str, r2: None = ..., lazy: Literal[True] = ...
) -> pl.LazyFrame: ...
@overload
def sink_fastq(
    lf: pl.LazyFrame, *, r1: str, r2: str | None, lazy: Literal[False]
) -> None: ...


def sink_fastq(
    lf: pl.LazyFrame,
    *,
    r1: str,
    r2: str | None = None,
    lazy: bool = False,
) -> None | pl.LazyFrame | tuple[pl.LazyFrame, pl.LazyFrame]:
    """Write a LazyFrame with _r1/_r2 suffixed columns as gzipped fastq files.

    Always sinks r1. Sinks r2 only if a path is given and _r2 columns are present.
    """
    schema_names = lf.collect_schema().names()

    def _prepare(suffix: str) -> pl.LazyFrame:
        cols = [c for c in schema_names if c.endswith(f"_{suffix}")]

        _strip_suffix = {c: c[: -(len(suffix) + 1)] for c in cols}
        _create_read: pl.Expr = pl.concat_str(
            [
                pl.format("@{} {}", pl.col("name"), pl.col("description")),
                pl.col("sequence"),
                pl.lit("+"),
                pl.col("quality"),
            ],
            separator="\n",
        ).alias("read")

        return lf.select(cols).rename(_strip_suffix).select(_create_read)

    has_r2 = any(c.endswith("_r2") for c in schema_names)
    if has_r2 and r2 is None:
        raise ValueError("LazyFrame is paired but no r2 output path was provided.")

    sink_kwargs = {
        "include_header": False,
        "quote_style": "never",
        "compression": "gzip",
    }

    sinks = [_prepare("r1").sink_csv(r1, **sink_kwargs)]
    if r2 is not None:
        sinks.append(_prepare("r2").sink_csv(r2, **sink_kwargs))

    if lazy:
        return sinks[0] if len(sinks) == 0 else tuple(sinks)

    pl.collect_all(sinks)


def scan_fastq(r1: str, r2: str | None) -> pl.LazyFrame:
    """Scan a (paired) fastq.gz into a polars lazyframe."""

    def _scan_read(fq: str, read: ReadType) -> pl.LazyFrame:
        lf = ox.from_fastq(fq).to_polars(lazy=True)
        assert isinstance(lf, pl.LazyFrame)
        return lf.rename({col: f"{col}_{read}" for col in lf.collect_schema().names()})

    lf = _scan_read(r1, "r1")
    if r2 is not None:
        lf = pl.concat([lf, _scan_read(r2, "r2")], how="horizontal")

    return lf
