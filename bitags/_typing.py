from typing import Callable, Literal

import polars as pl

type MoleculeType = Literal["dna", "rna"]
type ReadType = Literal["r1", "r2"]
type TagParser = Callable[[pl.LazyFrame], pl.LazyFrame]
type TagMerger = Callable[[pl.LazyFrame], pl.LazyFrame]
