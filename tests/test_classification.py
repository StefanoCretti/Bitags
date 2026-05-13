import polars as pl
import pytest

from bitags._classification import classify_reads


@pytest.mark.parametrize("read", ["r1", "r2"])
def test_classify_read_type(paired_tagged_reads, classification_regexes, read):
    result = classify_reads(
        paired_tagged_reads,
        classification_regexes,
        col=f"tag_type_{read}",
        into="category",
    ).collect()
    expected = paired_tagged_reads.collect()[f"category_{read}"].to_list()
    assert result["category"].to_list() == expected


def test_classify_custom_other_label(paired_tagged_reads, classification_regexes):
    result = classify_reads(
        paired_tagged_reads,
        classification_regexes,
        col="tag_type_r1",
        into="category",
        other="Unknown",
    ).collect()
    assert "Unknown" in result["category"]
    assert "Other" not in result["category"]


def test_classify_result_is_categorical(paired_tagged_reads, classification_regexes):
    result = classify_reads(
        paired_tagged_reads,
        classification_regexes,
        col="tag_type_r1",
        into="category",
    ).collect()
    assert result["category"].dtype == pl.Categorical


def test_classify_preserves_existing_columns(
    paired_tagged_reads, classification_regexes
):
    original_cols = paired_tagged_reads.collect_schema().names()
    result = classify_reads(
        paired_tagged_reads,
        classification_regexes,
        col="tag_type_r1",
        into="category",
    ).collect()
    assert all(c in result.columns for c in original_cols)


def test_classify_into_column_name(paired_tagged_reads, classification_regexes):
    result = classify_reads(
        paired_tagged_reads,
        classification_regexes,
        col="tag_type_r1",
        into="read_type",
    ).collect()
    assert "read_type" in result.columns
    assert "category" not in result.columns


def test_classify_empty_regexes_raises(paired_tagged_reads):
    with pytest.raises(ValueError):
        classify_reads(
            paired_tagged_reads, {}, col="tag_type_r1", into="category"
        ).collect()
