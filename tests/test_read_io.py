import polars as pl
import polars.selectors as cs
import pytest

from bitags._read_io import scan_fastq, sink_fastq

# --- scan_fastq ---


@pytest.mark.parametrize("read,fastq_fixture", [("r1", "fastq_r1"), ("r2", "fastq_r2")])
def test_scan_fastq_unpaired(request, read, fastq_fixture, paired_reads):
    fq = request.getfixturevalue(fastq_fixture)
    result = scan_fastq(fq).collect()
    expected = paired_reads.select(cs.ends_with(f"_{read}")).collect()
    result = result.rename({c: c.replace("_r1", f"_{read}") for c in result.columns})
    assert result.equals(expected)


def test_scan_fastq_paired(fastq_r1, fastq_r2, paired_reads):
    result = scan_fastq(fastq_r1, fastq_r2).collect()
    assert result.equals(paired_reads.collect())


# --- sink_fastq ---


def test_sink_fastq_raises_when_paired_but_no_r2_path(paired_reads, tmp_path):
    with pytest.raises(ValueError, match="paired"):
        sink_fastq(paired_reads, r1=str(tmp_path / "out.fastq.gz"))


def test_sink_fastq_raises_when_unpaired_but_r2_path(unpaired_reads, tmp_path):
    with pytest.raises(ValueError, match="unpaired"):
        sink_fastq(
            unpaired_reads,
            r1=str(tmp_path / "out_r1.fastq.gz"),
            r2=str(tmp_path / "out_r2.fastq.gz"),
        )


def test_sink_fastq_roundtrip_unpaired(unpaired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    sink_fastq(unpaired_reads, r1=out_r1)
    result = scan_fastq(out_r1).collect()
    assert result.equals(unpaired_reads.collect())


def test_sink_fastq_roundtrip_paired(paired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    out_r2 = str(tmp_path / "out_r2.fastq.gz")
    sink_fastq(paired_reads, r1=out_r1, r2=out_r2)
    result = scan_fastq(out_r1, out_r2).collect()
    assert result.equals(paired_reads.collect())


# --- lazy sinks ---


def test_sink_fastq_lazy_unpaired_returns_lazyframe(unpaired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    result = sink_fastq(unpaired_reads, r1=out_r1, lazy=True)
    assert isinstance(result, pl.LazyFrame)


def test_sink_fastq_lazy_paired_returns_tuple(paired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    out_r2 = str(tmp_path / "out_r2.fastq.gz")
    result = sink_fastq(paired_reads, r1=out_r1, r2=out_r2, lazy=True)
    assert isinstance(result, tuple)
    assert len(result) == 2


def test_sink_fastq_lazy_unpaired_roundtrip(unpaired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    sink_fastq(unpaired_reads, r1=out_r1, lazy=True).collect()
    result = scan_fastq(out_r1).collect()
    assert result.equals(unpaired_reads.collect())


def test_sink_fastq_lazy_paired_roundtrip(paired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    out_r2 = str(tmp_path / "out_r2.fastq.gz")
    s1, s2 = sink_fastq(paired_reads, r1=out_r1, r2=out_r2, lazy=True)
    pl.collect_all([s1, s2])
    result = scan_fastq(out_r1, out_r2).collect()
    assert result.equals(paired_reads.collect())


# --- tags ---


def test_sink_fastq_tags_unpaired(unpaired_tagged_reads, tmp_path):
    out = str(tmp_path / "out_r1.fastq.gz")
    sink_fastq(unpaired_tagged_reads, r1=out, tags=[("CB", "Z", "tag_seq_r1")])
    base = unpaired_tagged_reads.collect()
    result = scan_fastq(out).collect()
    expected = (base["description_r1"] + "\tCB:Z:" + base["tag_seq_r1"]).to_list()
    assert result["description_r1"].to_list() == expected


def test_sink_fastq_tags_paired(paired_tagged_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    out_r2 = str(tmp_path / "out_r2.fastq.gz")
    sink_fastq(
        paired_tagged_reads, r1=out_r1, r2=out_r2, tags=[("CB", "Z", "tag_seq_r1")]
    )
    base = paired_tagged_reads.collect()
    assert (
        scan_fastq(out_r1).collect()["description_r1"].to_list()
        == (base["description_r1"] + "\tCB:Z:" + base["tag_seq_r1"]).to_list()
    )
    assert (
        scan_fastq(out_r2).collect()["description_r1"].to_list()
        == base["description_r2"].to_list()
    )


def test_sink_fastq_tags_empty_description(unpaired_tagged_reads, tmp_path):
    out = str(tmp_path / "out_r1.fastq.gz")
    lf = unpaired_tagged_reads.with_columns(pl.lit("").alias("description_r1"))
    sink_fastq(lf, r1=out, tags=[("CB", "Z", "tag_seq_r1")])
    base = unpaired_tagged_reads.collect()
    result = scan_fastq(out).collect()
    expected = ("CB:Z:" + base["tag_seq_r1"]).to_list()
    assert result["description_r1"].to_list() == expected
