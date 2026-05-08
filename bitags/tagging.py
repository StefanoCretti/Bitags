import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import polars as pl

from bitags._typing import ReadType


@dataclass
class Tag:
    """Tag descriptor for bitap fuzzy matching."""

    seq: str
    info: str
    max_mism: int = field(default=0)


def load_tags(
    path: str | Path,
    *,
    info_col: int = 0,
    seq_col: int = 2,
    mism_col: int = 3,
) -> list[Tag]:
    """Load tags from a TSV file. Skips lines starting with '#'."""
    tags: list[Tag] = []
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            tags.append(Tag(seq=row[seq_col], info=row[info_col], max_mism=int(row[mism_col])))
    return tags


def barcode_reads(
    lf: pl.LazyFrame,
    tags: Sequence[Tag],
    *,
    read: ReadType,
) -> pl.LazyFrame:
    """Add tag_seq_{read}, tag_type_{read}, tag_pos_{read} columns via bitap matching.

    Applies the Rust plugin to sequence_{read}, finding all tags and resolving
    overlaps. Output columns match the schema used by the rest of the pipeline.
    """
    from . import find_tags

    return (
        lf.with_columns(find_tags(pl.col(f"sequence_{read}"), tags).alias("_tags"))
        .unnest("_tags")
        .rename(
            {
                "tag_seq": f"tag_seq_{read}",
                "tag_type": f"tag_type_{read}",
                "tag_pos": f"tag_pos_{read}",
            }
        )
    )
