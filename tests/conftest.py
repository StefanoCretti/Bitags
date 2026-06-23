import gzip
import json
from pathlib import Path

import polars as pl
import polars.selectors as cs
import pytest


@pytest.fixture
def paired_tagged_reads() -> pl.LazyFrame:
    path = Path(__file__).parent / "data" / "paired_tagged_reads.tsv"
    return (
        pl.read_csv(path, separator="\t", infer_schema_length=0)
        .with_columns(pl.col(pl.String).fill_null(""))
        .with_columns(
            pl.col("tag_pos_r1")
            .str.split(" ")
            .list.eval(pl.element().filter(pl.element() != "").cast(pl.UInt32)),
            pl.col("tag_pos_r2")
            .str.split(" ")
            .list.eval(pl.element().filter(pl.element() != "").cast(pl.UInt32)),
        )
        .lazy()
    )


@pytest.fixture
def unpaired_tagged_reads(paired_tagged_reads) -> pl.LazyFrame:
    return paired_tagged_reads.drop(cs.ends_with("_r2"))


@pytest.fixture
def paired_reads(paired_tagged_reads) -> pl.LazyFrame:
    return paired_tagged_reads.drop(cs.starts_with("tag") | cs.starts_with("category"))


@pytest.fixture
def unpaired_reads(unpaired_tagged_reads) -> pl.LazyFrame:
    return unpaired_tagged_reads.drop(
        cs.starts_with("tag") | cs.starts_with("category")
    )


def _lf_to_fastq(lf: pl.LazyFrame, read: str) -> str:
    df = lf.collect()
    lines = []
    for row in df.iter_rows(named=True):
        lines.append(f"@{row[f'name_{read}']} {row[f'description_{read}']}")
        lines.append(row[f"sequence_{read}"])
        lines.append("+")
        lines.append(row[f"quality_{read}"])
    return "\n".join(lines) + "\n"


@pytest.fixture
def fastq_r1(paired_reads, tmp_path) -> str:
    path = tmp_path / "r1.fastq.gz"
    with gzip.open(path, "wt") as f:
        f.write(_lf_to_fastq(paired_reads, "r1"))
    return str(path)


@pytest.fixture
def fastq_r2(paired_reads, tmp_path) -> str:
    path = tmp_path / "r2.fastq.gz"
    with gzip.open(path, "wt") as f:
        f.write(_lf_to_fastq(paired_reads, "r2"))
    return str(path)


@pytest.fixture
def tags() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "seq": ["AACAA", "CCCCC", "GGGGG", "CCCGG"],
            "info": ["TagA", "TagDNA", "TagRNA", "TagB"],
            "max_mism": [1, 0, 0, 0],
        }
    )


@pytest.fixture
def classification_regexes() -> dict[str, str]:
    """To match, reads must end with the appropriate tag."""
    return {
        "DNA": r"^.*TagDNA$",
        "RNA": r"^.*TagRNA$",
        "Any": r"^.*Tag(DNA|RNA)$",  # Here just to check that the first match is chosen
    }


@pytest.fixture
def trimming_cases() -> dict[str, tuple[str, pl.DataFrame]]:
    path = Path(__file__).parent / "data" / "trimming_cases.tsv"
    df = (
        pl.read_csv(path, separator="\t", infer_schema_length=0)
        .with_columns(pl.col(pl.String).fill_null(""))
        .with_columns(
            pl.col("tag_pos_r1")
            .str.split(" ")
            .list.eval(pl.element().filter(pl.element() != "").cast(pl.UInt32))
        )
    )
    return {
        key[0]: (group["regex"][0], group.drop("key", "regex"))
        for key, group in df.group_by("key", maintain_order=True)
    }


@pytest.fixture
def parquet_paired_reads(paired_reads, tmp_path) -> str:
    path = str(tmp_path / "paired_reads.parquet")
    paired_reads.sink_parquet(path)
    return path


@pytest.fixture
def parquet_unpaired_reads(unpaired_reads, tmp_path) -> str:
    path = str(tmp_path / "unpaired_reads.parquet")
    unpaired_reads.sink_parquet(path)
    return path


@pytest.fixture
def parquet_unpaired_tagged_reads(unpaired_tagged_reads, tmp_path) -> str:
    path = str(tmp_path / "unpaired_tagged_reads.parquet")
    unpaired_tagged_reads.sink_parquet(path)
    return path


@pytest.fixture
def parquet_paired_tagged_reads(paired_tagged_reads, tmp_path) -> str:
    path = str(tmp_path / "paired_tagged_reads.parquet")
    paired_tagged_reads.sink_parquet(path)
    return path


@pytest.fixture
def tags_tsv(tags, tmp_path) -> str:
    path = str(tmp_path / "tags.tsv")
    tags.write_csv(path, separator="\t")
    return path


@pytest.fixture
def classification_json(classification_regexes, tmp_path) -> str:
    path = str(tmp_path / "classification.json")
    with open(path, "w") as f:
        json.dump(classification_regexes, f)
    return path


@pytest.fixture
def malformed_trimming_regexes() -> dict[str, str]:
    """Should fail validation and raise value error."""
    return {
        "0_groups": r"^TagA.*?$",
        "1_groups": r"^(TagA).*?$",
        "3_groups": r"^(TagA).*?()()$",
        "^_missing": r"(TagA).*?()$",
        "$_missing": r"^(TagA).*?()",
        "^$_missing": r"(TagA).*?()",
        "empty_string": r"",
    }
