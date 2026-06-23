import re

import polars as pl

from bitags._schema import get_read_cols
from bitags._typing import ReadType


def trim_reads(
    lf: pl.LazyFrame,
    regex: str,
    *,
    read: ReadType = "r1",
    tags_only: bool = False,
    fill_empty: bool = True,
) -> pl.LazyFrame:
    """Trim reads using the two capture groups of the regex matched against tag_type.

    Group 1 defines tags to remove from the 5' end, group 2 from the 3' end.
    When tags_only=False (default), sequence and quality are also trimmed and tag
    positions are adjusted relative to the new sequence start. When tags_only=True,
    only tag_seq, tag_type, and tag_pos are trimmed; sequence and quality are untouched.
    When fill_empty=True (default), reads whose trimmed sequence is empty are given
    a single 'N' in the sequence column and a single 'I' in the quality column.
    This prevents empty sequences, which are not valid in FASTQ format.
    """
    if not (regex.startswith("^") and regex.endswith("$")):
        raise ValueError("regex must have start (^) and end ($) anchors")
    if re.compile(regex).groups != 2:
        raise ValueError("regex must have exactly two capturing groups")

    cols = get_read_cols(read)
    original_cols = lf.collect_schema().names()

    def _tag_bounds(lf: pl.LazyFrame, regex: str) -> pl.LazyFrame:
        """Match regex against tag_type and add _first_tag/_last_tag columns.

        _first_tag: index of the first kept tag (left-inclusive).
        _last_tag:  index past the last kept tag (right-exclusive).
        """

        def tag_count(c: str) -> pl.Expr:
            return (
                pl.when(pl.col(c).is_null() | (pl.col(c) == ""))
                .then(pl.lit(0, dtype=pl.UInt32))
                .otherwise(pl.col(c).str.split(":").list.len())
            )

        return (
            lf.with_columns(
                pl.col(cols.tag_type)
                .str.extract_groups(regex)
                .struct.rename_fields(["_left", "_right"])
                .alias("_groups")
            )
            .unnest("_groups")
            .with_columns(
                tag_count("_left").alias("_first_tag"),
                (tag_count(cols.tag_type) - tag_count("_right")).alias("_last_tag"),
            )
            .drop("_left", "_right")
        )

    def _trim_tag_cols(lf: pl.LazyFrame) -> pl.LazyFrame:
        """Slice tag_seq, tag_type, tag_pos in-place using _first_tag/_last_tag."""
        tag_slice = (
            pl.col("_first_tag"),
            pl.col("_last_tag") - pl.col("_first_tag"),
        )
        return lf.with_columns(
            pl.col(cols.tag_pos).list.slice(*tag_slice),
            pl.col(cols.tag_seq).str.split(":").list.slice(*tag_slice).list.join(":"),
            pl.col(cols.tag_type).str.split(":").list.slice(*tag_slice).list.join(":"),
        )

    def _seq_bounds(lf: pl.LazyFrame) -> pl.LazyFrame:
        """Add _first_base/_last_base using the original (pre-trim) tag columns.

        _first_base: sequence position right after the last left-trimmed tag ends.
        _last_base:  sequence position where the first right-trimmed tag starts.
        """
        return lf.with_columns(
            (
                pl.when(pl.col("_first_tag") == 0)
                .then(pl.lit(0, dtype=pl.UInt32))
                .otherwise(
                    pl.col(cols.tag_pos).list.get(
                        pl.col("_first_tag") - 1, null_on_oob=True
                    )
                    + pl.col(cols.tag_seq)
                    .str.split(":")
                    .list.get(pl.col("_first_tag") - 1, null_on_oob=True)
                    .str.len_chars()
                )
            ).alias("_first_base"),
            (
                pl.when(pl.col("_last_tag") >= pl.col(cols.tag_pos).list.len())
                .then(pl.col(cols.seq).str.len_chars())
                .otherwise(
                    pl.col(cols.tag_pos).list.get(pl.col("_last_tag"), null_on_oob=True)
                )
            ).alias("_last_base"),
        )

    def _trim_seq_cols(lf: pl.LazyFrame) -> pl.LazyFrame:
        """Slice sequence and quality using _first_base/_last_base, adjusting tag positions."""
        seq_slice = pl.col("_first_base"), pl.col("_last_base") - pl.col("_first_base")
        return lf.with_columns(
            pl.col(cols.tag_pos) - pl.col("_first_base"),
            pl.col(cols.seq).str.slice(*seq_slice),
            pl.col(cols.qual).str.slice(*seq_slice),
        )

    def _fill_empty_seq(lf: pl.LazyFrame) -> pl.LazyFrame:
        empty = pl.col(cols.seq).str.len_chars() == 0
        return lf.with_columns(
            pl.when(empty)
            .then(pl.lit("N"))
            .otherwise(pl.col(cols.seq))
            .alias(cols.seq),
            pl.when(empty)
            .then(pl.lit("I"))
            .otherwise(pl.col(cols.qual))
            .alias(cols.qual),
        )

    # Note: order is a bit finicky since the expressions are dependent on each other
    # - _seq_bounds must run before _trim_tag_cols since it reads the original tag positions
    # - _trim_seq_cols must run after _trim_tag_cols since it adjusts tag position
    if tags_only:
        return _tag_bounds(lf, regex).pipe(_trim_tag_cols).select(original_cols)
    result = (
        _tag_bounds(lf, regex)
        .pipe(_seq_bounds)
        .pipe(_trim_tag_cols)
        .pipe(_trim_seq_cols)
        .select(original_cols)
    )
    if fill_empty:
        result = result.pipe(_fill_empty_seq)
    return result
