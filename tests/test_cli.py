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


def test_parquet_to_fastq_tags(
    parquet_paired_tagged_reads, paired_tagged_reads, tmp_path
):
    out_r1 = str(tmp_path / "out_r1.fastq.gz")
    out_r2 = str(tmp_path / "out_r2.fastq.gz")
    result = CliRunner().invoke(
        cli,
        [
            "parquet-to-fastq",
            parquet_paired_tagged_reads,
            "--r1",
            out_r1,
            "--r2",
            out_r2,
            "-t",
            "CB",
            "Z",
            "tag_seq_r1",
        ],
    )
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    assert (
        scan_fastq(out_r1).collect()["description_r1"].to_list()
        == (base["description_r1"] + "\tCB:Z:" + base["tag_seq_r1"]).to_list()
    )


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


def test_trim_min_length(
    parquet_paired_tagged_reads, paired_tagged_reads, trimming_cases, tmp_path
):
    regex, overrides = trimming_cases["trim_one_5'"]
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli, ["trim", parquet_paired_tagged_reads, out, "-r", regex, "-m", "10"]
    )
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    expected = base.with_columns(overrides[col] for col in overrides.columns).filter(
        pl.col("sequence_r1").str.len_chars() >= 10
    )
    assert pl.read_parquet(out).equals(expected)


def test_trim_min_length_with_excluded(
    parquet_paired_tagged_reads, paired_tagged_reads, trimming_cases, tmp_path
):
    regex, overrides = trimming_cases["trim_one_5'"]
    out = str(tmp_path / "out.parquet")
    excluded = str(tmp_path / "excluded.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "trim",
            parquet_paired_tagged_reads,
            out,
            "-r",
            regex,
            "-m",
            "10",
            "-e",
            excluded,
        ],
    )
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    trimmed = base.with_columns(overrides[col] for col in overrides.columns)
    assert pl.read_parquet(out).equals(
        trimmed.filter(pl.col("sequence_r1").str.len_chars() >= 10)
    )
    assert pl.read_parquet(excluded).equals(
        trimmed.filter(pl.col("sequence_r1").str.len_chars() < 10)
    )


def test_trim_excluded_without_min_length_exits_nonzero(
    parquet_paired_tagged_reads, tmp_path
):
    result = CliRunner().invoke(
        cli,
        [
            "trim",
            parquet_paired_tagged_reads,
            str(tmp_path / "out.parquet"),
            "-r",
            "^(TagA).*?()$",
            "-e",
            str(tmp_path / "excluded.parquet"),
        ],
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
            "-j",
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


@pytest.mark.parametrize("read", ["r1", "r2"])
def test_classify_inline_categories(
    parquet_paired_tagged_reads,
    paired_tagged_reads,
    classification_regexes,
    tmp_path,
    read,
):
    out = str(tmp_path / "out.parquet")
    args = ["classify", parquet_paired_tagged_reads, out, "--read", read]
    for label, pattern in classification_regexes.items():
        args += ["-c", label, pattern]
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0
    base = paired_tagged_reads.collect()
    expected = base.with_columns(
        pl.col(f"category_{read}").cast(pl.Categorical).alias("category")
    )
    assert pl.read_parquet(out).sort("name_r1").equals(expected.sort("name_r1"))


def test_classify_inline_and_json_mutually_exclusive(
    parquet_paired_tagged_reads, classification_json, tmp_path
):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "classify",
            parquet_paired_tagged_reads,
            out,
            "-j",
            classification_json,
            "-c",
            "DNA",
            ".*",
        ],
    )
    assert result.exit_code != 0
    assert "not both" in result.output


def test_classify_no_regexes_fails(parquet_paired_tagged_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(cli, ["classify", parquet_paired_tagged_reads, out])
    assert result.exit_code != 0


def test_classify_custom_other_label(
    parquet_paired_tagged_reads, classification_json, tmp_path
):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli,
        [
            "classify",
            parquet_paired_tagged_reads,
            out,
            "-j",
            classification_json,
            "-o",
            "Unassigned",
        ],
    )
    assert result.exit_code == 0
    categories = pl.read_parquet(out)["category"].cast(pl.String).unique().to_list()
    assert "Other" not in categories


def test_classify_category_missing_pattern_fails(parquet_paired_tagged_reads, tmp_path):
    out = str(tmp_path / "out.parquet")
    result = CliRunner().invoke(
        cli, ["classify", parquet_paired_tagged_reads, out, "-c", "DNA"]
    )
    assert result.exit_code != 0
