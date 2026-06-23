import pytest

from bitags._manipulation import to_unpaired


def test_to_unpaired_drop_r2(paired_tagged_reads):
    result = to_unpaired(paired_tagged_reads, "r2").collect()
    assert not any(c.endswith("_r2") for c in result.columns)
    assert any(c.endswith("_r1") for c in result.columns)


def test_to_unpaired_drop_r1_renames_to_r1(paired_tagged_reads):
    result = to_unpaired(paired_tagged_reads, "r1").collect()
    assert not any(c.endswith("_r2") for c in result.columns)
    assert any(c.endswith("_r1") for c in result.columns)


def test_to_unpaired_drop_r1_preserves_values(paired_tagged_reads):
    result = to_unpaired(paired_tagged_reads, "r1").collect()
    expected = paired_tagged_reads.collect()["sequence_r2"].to_list()
    assert result["sequence_r1"].to_list() == expected


def test_to_unpaired_raises_if_not_paired(unpaired_tagged_reads):
    with pytest.raises(ValueError, match="not paired"):
        to_unpaired(unpaired_tagged_reads, "r2").collect()
