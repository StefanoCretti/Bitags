from typing import Mapping

import polars as pl


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
    assert isinstance(expr, pl.Expr), "At least one regex is required."
    return lf.with_columns(
        expr.otherwise(pl.lit(other)).cast(pl.Categorical).alias(into)
    )
