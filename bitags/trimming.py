from types import SimpleNamespace

from bitags._typing import ReadType

import polars as pl


def _get_read_cols(read: str) -> SimpleNamespace:
    """Return column name attributes for the given read suffix."""
    return SimpleNamespace(
        seq=f"sequence_{read}",
        qual=f"quality_{read}",
        tag_seq=f"tag_seq_{read}",
        tag_type=f"tag_type_{read}",
        tag_pos=f"tag_pos_{read}",
    )


def _add_trim_counts(lf: pl.LazyFrame, regex: str, col: str) -> pl.LazyFrame:
    """Match a two-group regex against col and add _n_left/_n_right tag-count columns."""

    def _count(c: str) -> pl.Expr:
        """Count colon-separated items in a capture group; returns 0 if null/empty."""
        return (
            pl.when(pl.col(c).is_null() | (pl.col(c) == ""))
            .then(pl.lit(0, dtype=pl.UInt32))
            .otherwise(pl.col(c).str.split(":").list.len())
        )

    return (
        lf.with_columns(
            pl.col(col)
            .str.extract_groups(regex)
            .struct.rename_fields(["_left", "_right"])
            .alias("_groups")
        )
        .unnest("_groups")
        .with_columns(_n_left=_count("_left"), _n_right=_count("_right"))
        .drop("_left", "_right")
    )


def extract_barcode(
    lf: pl.LazyFrame,
    regex: str,
    *,
    read: ReadType = "r2",
    into: str = "barcode",
) -> pl.LazyFrame:
    """Extract the barcode portion of tag_seq using the same regex format as trim_reads.

    Capture group 1 = left tags to skip, capture group 2 = right tags to skip.
    The remaining colon-joined tag_seq values are stored in `into`.
    """
    cols = _get_read_cols(read)
    split_seq = pl.col(cols.tag_seq).str.split(":")
    return (
        _add_trim_counts(lf, regex, cols.tag_type)
        .with_columns(
            split_seq.list.slice(
                pl.col("_n_left"),
                split_seq.list.len() - pl.col("_n_left") - pl.col("_n_right"),
            )
            .list.join(":")
            .alias(into)
        )
        .drop("_n_left", "_n_right")
    )


def trim_reads(
    lf: pl.LazyFrame,
    regex: str,
    *,
    read: ReadType = "r1",
) -> pl.LazyFrame:
    """Trim reads by removing tags matched by the two capture groups of the regex.

    The regex is matched against tag_type. The first capture group defines tags
    to remove from the 5' end, the second from the 3' end. Either group may be
    absent (no trimming on that side). Sequence, quality, and tag positions are
    sliced and adjusted accordingly.
    """
    cols = _get_read_cols(read)

    def _seq_slice_bounds(lf: pl.LazyFrame) -> pl.LazyFrame:
        """Add _seq_start and _seq_end columns defining the region to keep.

        _seq_start: position right after the last left-trimmed tag ends.
        _seq_end:   position where the first right-trimmed tag starts.
        """
        return lf.with_columns(
            (
                pl.when(pl.col("_n_left") == 0)
                .then(pl.lit(0, dtype=pl.UInt32))
                .otherwise(
                    pl.col(cols.tag_pos).list.get(
                        pl.col("_n_left") - 1, null_on_oob=True
                    )
                    + pl.col(cols.tag_seq)
                    .str.split(":")
                    .list.get(pl.col("_n_left") - 1, null_on_oob=True)
                    .str.len_chars()
                )
            ).alias("_seq_start"),
            (
                pl.when(pl.col("_n_right") == 0)
                .then(pl.col(cols.seq).str.len_chars())
                .otherwise(
                    pl.col(cols.tag_pos).list.get(
                        -pl.col("_n_right").cast(pl.Int32), null_on_oob=True
                    )
                )
            ).alias("_seq_end"),
        )

    def _apply_slices(lf: pl.LazyFrame) -> pl.LazyFrame:
        """Slice sequence, quality, and tag columns, adjusting tag positions accordingly."""
        tag_slice = (
            pl.col("_n_left"),
            pl.col(cols.tag_pos).list.len() - pl.col("_n_left") - pl.col("_n_right"),
        )
        seq_slice = pl.col("_seq_start"), pl.col("_seq_end") - pl.col("_seq_start")

        return lf.with_columns(
            pl.col(cols.tag_pos).list.slice(*tag_slice) - pl.col("_seq_start"),
            pl.col(cols.tag_seq).str.split(":").list.slice(*tag_slice).list.join(":"),
            pl.col(cols.tag_type).str.split(":").list.slice(*tag_slice).list.join(":"),
            pl.col(cols.seq).str.slice(*seq_slice),
            pl.col(cols.qual).str.slice(*seq_slice),
        )

    original_get_read_cols = lf.collect_schema().names()

    return (
        _add_trim_counts(lf, regex, cols.tag_type)
        .pipe(_seq_slice_bounds)
        .pipe(_apply_slices)
        .select(original_get_read_cols)
    )
