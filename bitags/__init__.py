from pathlib import Path
from typing import Sequence

import polars as pl
from polars.plugins import register_plugin_function

_LIB = Path(__file__).parent


def find_tags(expr: pl.Expr, tags: Sequence) -> pl.Expr:
    """Find all tags in each sequence string via bitap fuzzy matching.

    Returns a struct column with fields tag_seq (str), tag_type (str), tag_pos (List[u32]).
    Each tag must have .seq, .info, and .max_mism attributes.
    """
    return register_plugin_function(
        plugin_path=_LIB,
        function_name="find_tags",
        args=[expr],
        kwargs={
            "tags": [
                {"seq": t.seq, "info": t.info, "max_mism": t.max_mism}
                for t in tags
            ]
        },
        is_elementwise=True,
    )
