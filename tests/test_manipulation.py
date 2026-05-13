import pytest

from bitags._manipulation import embed_barcode, to_unpaired


@pytest.mark.parametrize(
    "source,target",
    [
        ("r1", "r1"),
        ("r1", "r2"),
        ("r2", "r1"),
        ("r2", "r2"),
    ],
)
def test_embed_barcode_combinations(paired_tagged_reads, source, target):
    result = embed_barcode(
        paired_tagged_reads, source_read=source, target_read=target
    ).collect()
    df = paired_tagged_reads.collect()
    expected = (df[f"name_{target}"] + "::" + df[f"tag_seq_{source}"]).to_list()
    assert result[f"name_{target}"].to_list() == expected


@pytest.mark.parametrize("sep", ["|", "--"])
def test_embed_barcode_separator(paired_tagged_reads, sep):
    result = embed_barcode(
        paired_tagged_reads, source_read="r1", target_read="r2", separator=sep
    ).collect()
    df = paired_tagged_reads.collect()
    expected = (df["name_r2"] + sep + df["tag_seq_r1"]).to_list()
    assert result["name_r2"].to_list() == expected


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
