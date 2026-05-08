from pathlib import Path

import polars as pl
from polars.plugins import register_plugin_function

from bitags._typing import ReadType

_LIB = Path(__file__).parent


def find_tags(expr: pl.Expr, tags: pl.DataFrame) -> pl.Expr:
    """Find all tags in each sequence string via bitap fuzzy matching.

    Returns a struct column with fields tag_seq (str), tag_type (str), tag_pos (List[u32]).
    """
    return register_plugin_function(
        plugin_path=_LIB,
        function_name="find_tags",
        args=[expr],
        kwargs={"tags": tags.to_dicts()},
        is_elementwise=True,
    )


def barcode_reads(
    lf: pl.LazyFrame,
    tags: pl.DataFrame,
    *,
    read: ReadType,
) -> pl.LazyFrame:
    """Add tag_seq_{read}, tag_type_{read}, tag_pos_{read} columns via bitap matching.

    Applies the Rust plugin to sequence_{read}, finding all tags and resolving
    overlaps. Output columns match the schema used by the rest of the pipeline.
    """
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
