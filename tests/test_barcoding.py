import pytest

from bitags._barcoding import barcode_reads

# Very light testing as there are some on the rust side


@pytest.mark.parametrize("read", ["r1", "r2"])
def test_barcode_read(paired_reads, paired_tagged_reads, tags, read):
    result = barcode_reads(paired_reads, tags, read=read).collect()
    added = set(result.columns) - set(paired_reads.collect_schema().names())
    assert added == {f"tag_seq_{read}", f"tag_type_{read}", f"tag_pos_{read}"}
    expected = paired_tagged_reads.collect().select(result.columns)
    assert result.equals(expected)


def test_barcode_with_no_tags(paired_reads, tags):
    empty_tags = tags.clear()
    result = barcode_reads(paired_reads, empty_tags, read="r1").collect()
    n = len(result)
    assert result["tag_seq_r1"].to_list() == [""] * n
    assert result["tag_type_r1"].to_list() == [""] * n
    assert result["tag_pos_r1"].to_list() == [[]] * n
