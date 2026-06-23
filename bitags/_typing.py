from typing import Literal

type MoleculeType = Literal["dna", "rna"]
type ReadType = Literal["r1", "r2"]
type SamTag = tuple[str, str, str]  # (tag_name, tag_type, column_name)
