import json
from typing import Mapping

import polars as pl


def load_regexes(file: str) -> Mapping[str, str]:
    """Load regexes from a json file."""

    with open(file, "r") as handle:
        mapping = json.load(handle)

    assert all(isinstance(k, str) for k in mapping.keys())
    assert all(isinstance(v, str) for v in mapping.values())

    return mapping


def classify_reads(
    lf: pl.LazyFrame,
    regexes: Mapping[str, str],
    *,
    col: str,
    into: str,
    other: str = "Other",
) -> pl.LazyFrame:
    """Classify reads by matching regexes against a column in order.

    Each read is assigned the key of the first matching regex. Reads that
    match none are assigned the value of `other`. The result is stored as
    a Categorical column named `into`.
    """
    expr = pl
    for label, pattern in regexes.items():
        expr = expr.when(pl.col(col).str.contains(pattern)).then(pl.lit(label))
    assert isinstance(expr, pl.Expr)
    return lf.with_columns(
        expr.otherwise(pl.lit(other)).cast(pl.Categorical).alias(into)
    )


def split(
    lf: pl.LazyFrame,
    on: str | list[str],
    out_dir: str,
) -> pl.LazyFrame:
    """Split a LazyFrame into partitioned parquet files based on column values.

    Produces a hive-partitioned directory under `out_dir` (e.g. col=value/data.parquet).
    Uses a single sink pass — no upfront collection of unique values needed.
    """
    cols = [on] if isinstance(on, str) else on
    lf.sink_parquet(pl.PartitionBy(out_dir, key=cols), mkdir=True)
    return pl.scan_parquet(out_dir)
