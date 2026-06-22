import pytest

from bitags._trimming import trim_reads


@pytest.mark.parametrize(
    "key",
    [
        "trim_nothing",
        "trim_one_5'",
        "trim_one_3'",
        "trim_two_3'",
        "trim_if_both_5'3'",
        "trim_if_either_5'3'",
        "trim_two_3'_one_optional",
    ],
)
def test_trim(paired_tagged_reads, trimming_cases, key):
    regex, overrides = trimming_cases[key]
    result = trim_reads(paired_tagged_reads, regex).collect()
    base = paired_tagged_reads.collect()
    expected = base.with_columns(overrides[col] for col in overrides.columns)
    assert result.equals(expected)


@pytest.mark.parametrize(
    "key",
    [
        "0_groups",
        "1_groups",
        "3_groups",
        "^_missing",
        "$_missing",
        "^$_missing",
        "empty_string",
    ],
)
def test_trim_malformed_regex_raises(
    paired_tagged_reads, malformed_trimming_regexes, key
):
    with pytest.raises(ValueError):
        trim_reads(paired_tagged_reads, malformed_trimming_regexes[key])


def test_trim_fill_empty_false_leaves_empty_sequences(paired_tagged_reads):
    result = trim_reads(
        paired_tagged_reads, r"^(TagA).*?(TagDNA)$", fill_empty=False
    ).collect()
    assert (result["sequence_r1"].str.len_chars() == 0).any()
