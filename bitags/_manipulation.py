import polars as pl

from bitags._typing import ReadType


def embed_barcode(
    lf: pl.LazyFrame,
    source_read: ReadType,
    target_read: ReadType,
    separator: str = "::",
) -> pl.LazyFrame:
    """Append the barcode from source_read to the read name of target_read.

    Concatenates name_{target_read} and tag_seq_{source_read} using separator
    and writes the result back into name_{target_read}.
    """
    return lf.with_columns(
        pl.concat_str(
            [pl.col(f"name_{target_read}"), pl.col(f"tag_seq_{source_read}")],
            separator=separator,
        ).alias(f"name_{target_read}")
    )


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
