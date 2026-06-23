import polars as pl

from bitags._typing import ReadType


def to_unpaired(lf: pl.LazyFrame, read: ReadType) -> pl.LazyFrame:
    """Drop the given read and return an unpaired frame.

    Raises if the frame is not paired (both _r1 and _r2 columns must be present).
    If the dropped read was r1, the remaining r2 columns are renamed to r1.
    """
    cols = lf.collect_schema().names()
    has_r1 = any(c.endswith("_r1") for c in cols)
    has_r2 = any(c.endswith("_r2") for c in cols)
    if not (has_r1 and has_r2):
        raise ValueError("Reads are not paired.")

    lf = lf.drop(c for c in cols if c.endswith(f"_{read}"))

    if read == "r1":
        r2_cols = [c for c in lf.collect_schema().names() if c.endswith("_r2")]
        lf = lf.rename({c: c[:-2] + "r1" for c in r2_cols})

    return lf
