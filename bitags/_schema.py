from typing import NamedTuple

import polars as pl


class ReadCols(NamedTuple):
    name: str
    description: str
    seq: str
    qual: str
    tag_seq: str
    tag_type: str
    tag_pos: str


def get_read_cols(read: str) -> ReadCols:
    """Return column name attributes for the given read suffix."""
    return ReadCols(
        name=f"name_{read}",
        description=f"description_{read}",
        seq=f"sequence_{read}",
        qual=f"quality_{read}",
        tag_seq=f"tag_seq_{read}",
        tag_type=f"tag_type_{read}",
        tag_pos=f"tag_pos_{read}",
    )


def is_paired(lf: pl.LazyFrame) -> bool:
    """Return True if the LazyFrame contains the core read columns for both r1 and r2."""
    cols = set(lf.collect_schema().names())
    return all(set(get_read_cols(read)[:4]).issubset(cols) for read in ("r1", "r2"))
