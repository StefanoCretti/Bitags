import polars as pl
import pytest
from click.testing import CliRunner

from bitags._read_io import scan_fastq
from bitags.cli import cli


def test_fastq_to_parquet_unpaired(fastq_r1, unpaired_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(cli, ["fastq-to-parquet", fastq_r1, "-o", out])
    assert result.exit_code == 0
    assert pl.read_parquet(out).equals(unpaired_reads.collect())


def test_fastq_to_parquet_paired(fastq_r1, fastq_r2, paired_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli, ["fastq-to-parquet", fastq_r1, fastq_r2, "-o", out]
    )
    assert result.exit_code == 0
    assert pl.read_parquet(out).equals(paired_reads.collect())


def test_parquet_to_fastq_unpaired(parquet_unpaired_reads, unpaired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    result = CliRunner().invoke(
        cli, ["parquet-to-fastq", parquet_unpaired_reads, "--r1", out_r1]
    )
    assert result.exit_code == 0
    assert scan_fastq(out_r1).collect().equals(unpaired_reads.collect())


def test_parquet_to_fastq_unpaired_with_r2_path_exits_nonzero(
    parquet_unpaired_reads, tmp_path
):
    result = CliRunner().invoke(
        cli,
        [
            "parquet-to-fastq",
            parquet_unpaired_reads,
            "--r1",
            str(tmp_path / "out_r1.fastq.gz"),
            "--r2",
            str(tmp_path / "out_r2.fastq.gz"),
        ],
    )
    assert result.exit_code != 0


def test_parquet_to_fastq_paired_without_r2_path_exits_nonzero(
    parquet_paired_reads, tmp_path
):
    result = CliRunner().invoke(
        cli,
        [
            "parquet-to-fastq",
            parquet_paired_reads,
            "--r1",
            str(tmp_path / "out_r1.fastq.gz"),
        ],
    )
    assert result.exit_code != 0


def test_parquet_to_fastq_paired(parquet_paired_reads, paired_reads, tmp_path):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    out_r2 = str(tmp_path / "out_r2.fastq.gz")
    result = CliRunner().invoke(
        cli, ["parquet-to-fastq", parquet_paired_reads, "--r1", out_r1, "--r2", out_r2]
    )
    assert result.exit_code == 0
    assert scan_fastq(out_r1, out_r2).collect().equals(paired_reads.collect())


@pytest.mark.parametrize("read", ["r1", "r2"])
def test_barcode_single_read(
    parquet_paired_reads, tags_tsv, paired_tagged_reads, tmp_path, read
):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli, ["barcode", parquet_paired_reads, out, "-t", tags_tsv, "-r", read]
    )
    assert result.exit_code == 0
    df = pl.read_parquet(out)
    assert df.equals(paired_tagged_reads.collect().select(df.columns))


def test_barcode_both_reads(
    parquet_paired_reads, tags_tsv, paired_tagged_reads, tmp_path
):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli, ["barcode", parquet_paired_reads, out, "-t", tags_tsv]
    )
    assert result.exit_code == 0
    df = pl.read_parquet(out)
    assert df.equals(paired_tagged_reads.collect().select(df.columns))


@pytest.mark.parametrize("source,target", [("r1", "r2"), ("r2", "r1")])
def test_embed(
    parquet_paired_tagged_reads, paired_tagged_reads, tmp_path, source, target
):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "embed",
            parquet_paired_tagged_reads,
            out,
            "--source-read",
            source,
            "--target-read",
            target,
        ],
    )
    assert result.exit_code == 0
    df = pl.read_parquet(out)
    base = paired_tagged_reads.collect()
    expected_names = (
        base[f"name_{target}"] + "::" + base[f"tag_seq_{source}"]
    ).to_list()
    assert df["name_r1"].to_list() == expected_names
    assert not any(c.endswith("_r2") for c in df.columns)


def test_embed_no_remove_source(
    parquet_paired_tagged_reads, paired_tagged_reads, tmp_path
):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "embed",
            parquet_paired_tagged_reads,
            out,
            "--source-read",
            "r1",
            "--target-read",
            "r2",
            "--no-remove-source",
        ],
    )
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    expected = base.with_columns(
        (base["name_r2"] + "::" + base["tag_seq_r1"]).alias("name_r2")
    )
    assert pl.read_parquet(out).equals(expected)


def test_trim(
    parquet_paired_tagged_reads, paired_tagged_reads, trimming_cases, tmp_path
):
    regex, overrides = trimming_cases["trim_one_5'"]
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli, ["trim", parquet_paired_tagged_reads, out, "-r", regex]
    )
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    expected = base.with_columns(overrides[col] for col in overrides.columns)
    assert pl.read_parquet(out).equals(expected)


def test_trim_invalid_regex_exits_nonzero(parquet_paired_tagged_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli, ["trim", parquet_paired_tagged_reads, out, "-r", "^TagA.*?$"]
    )
    assert result.exit_code != 0


def test_filter_exact(parquet_paired_tagged_reads, paired_tagged_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "filter",
            parquet_paired_tagged_reads,
            out,
            "-c",
            "tag_type_r1",
            "-p",
            "TagA:TagDNA",
        ],
    )
    assert result.exit_code == 0
    expected = paired_tagged_reads.filter(
        pl.col("tag_type_r1") == "TagA:TagDNA"
    ).collect()
    assert pl.read_parquet(out).equals(expected)


def test_filter_exact_no_match_substring(parquet_paired_tagged_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "filter",
            parquet_paired_tagged_reads,
            out,
            "-c",
            "tag_type_r1",
            "-p",
            "TagDNA",
        ],
    )
    assert result.exit_code == 0
    assert len(pl.read_parquet(out)) == 0


def test_filter_exact_with_excluded(
    parquet_paired_tagged_reads, paired_tagged_reads, tmp_path
):
    out = str(tmp_path / "out.parquet")
    excluded = str(tmp_path / "excluded.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "filter",
            parquet_paired_tagged_reads,
            out,
            "-c",
            "tag_type_r1",
            "-p",
            "TagA:TagDNA",
            "-e",
            excluded,
        ],
    )
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    assert pl.read_parquet(out).equals(
        base.filter(pl.col("tag_type_r1") == "TagA:TagDNA")
    )
    assert pl.read_parquet(excluded).equals(
        base.filter(pl.col("tag_type_r1") != "TagA:TagDNA")
    )


def test_filter_regex(parquet_paired_tagged_reads, paired_tagged_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "filter",
            parquet_paired_tagged_reads,
            out,
            "-c",
            "tag_type_r1",
            "-p",
            "TagDNA",
            "--regex",
        ],
    )
    assert result.exit_code == 0
    expected = paired_tagged_reads.filter(
        pl.col("tag_type_r1").str.contains("TagDNA")
    ).collect()
    assert pl.read_parquet(out).equals(expected)


@pytest.mark.parametrize("read", ["r1", "r2"])
def test_classify(
    parquet_paired_tagged_reads,
    paired_tagged_reads,
    classification_json,
    tmp_path,
    read,
):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "classify",
            parquet_paired_tagged_reads,
            out,
            "-r",
            classification_json,
            "--read",
            read,
        ],
    )
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    expected = base.with_columns(
        pl.col(f"category_{read}").cast(pl.Categorical).alias("category")
    )
    assert pl.read_parquet(out).sort("name_r1").equals(expected.sort("name_r1"))
