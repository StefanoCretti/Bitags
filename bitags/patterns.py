from bitags._typing import MoleculeType

_R1_TAG = {"O": "Odd_r1", "E": "Even_r1", "T": "Term_r1"}
_R2_TAG = {"O": "Odd_r2", "E": "Even_r2", "T": "Term_r2"}
_LIGTAG_R1 = {"dna": "Dpm_r1", "rna": "Rpm_r1"}
_LIGTAG_R2 = {"dna": "Dpm_r2", "rna": "Rpm_r2"}


def _parse_schema(schema: str) -> list[str]:
    """Validate schema and return middle O/E tags followed by T, in schema order.

    Schema format: L[O|E]*T — starts with ligtag, optional alternating Odd/Even, ends with terminal.
    Empty middle (schema 'LT') is valid and represents a fully skipped round product.
    """
    if schema[0] != "L" or schema[-1] != "T":
        raise ValueError("Schema must start with 'L' (ligtag) and end with 'T' (terminal).")
    middle = schema[1:-1]
    if middle and not all(c in "OE" for c in middle):
        raise ValueError("Schema middle must consist only of 'O' and 'E'.")
    return [*middle, "T"]


def _nested_optional(tags: list[str], *, first_sep: bool = True) -> str:
    """Build a nested optional regex chain: (:t1(:t2(:t3)?)?)?

    With first_sep=False, the first tag has no leading colon: (t1(:t2(:t3)?)?)?
    """
    if not tags:
        return ""
    sep = ":" if first_sep else ""
    return f"({sep}{tags[0]}{_nested_optional(tags[1:])})" + "?"


def _capturing_chain(tags: list[str]) -> str:
    """Build a single capturing group over a nested optional chain.

    (t1(?::t2(?::t3)?)?) captures 't1', 't1:t2', or 't1:t2:t3'.
    Used by trim regexes so _count_tags can count the captured tags.
    """
    chain = tags[0]
    for tag in tags[1:]:
        chain = f"{chain}(?::{tag}"
    chain += ")?" * (len(tags) - 1)
    return f"({chain})"


def r2_regex(schema: str, molecule: MoleculeType) -> str:
    """Return a regex matching the expected tag_type sequence for R2.

    R2 reads from the barcode stack inward, so the order is:
    (Adap_l) → Term_r2 → [reversed Odd/Even] → ligtag → (DpmVar_l → DpmVar_r → Adap_r, DNA only).
    """
    tags = list(reversed(_parse_schema(schema)))
    core = ":".join([*[_R2_TAG[t] for t in tags], _LIGTAG_R2[molecule]])
    suffix = (
        _nested_optional(["DpmVar_l", "DpmVar_r", "Adap_r"])
        if molecule == "dna"
        else _nested_optional(["Adap_r"])
    )
    return f"^(Adap_l:)?{core}{suffix}$"


def r2_barcode_regex(schema: str, molecule: MoleculeType) -> str:
    """Return a regex for barcode extraction on R2, compatible with extract_barcode.

    Group 1 captures Adap_l (left to remove).
    Group 2 captures the DpmVar/Adap suffix chain (right to remove).
    The middle — Term_r2 through the ligtag — is kept.
    """
    tags = list(reversed(_parse_schema(schema)))
    core = ":".join([*[_R2_TAG[t] for t in tags], _LIGTAG_R2[molecule]])
    right = (
        "(?::(DpmVar_l(?::DpmVar_r(?::Adap_r)?)?))?$"
        if molecule == "dna"
        else "(?::(Adap_r))?$"
    )
    return f"^(Adap_l)?(?::)?{core}{right}"


def r2_trim_regex(schema: str, molecule: MoleculeType) -> str:
    """Return a regex for trimming all barcode tags from R2, compatible with trim_reads.

    Group 1 captures all barcode tags (Adap_l + core + optional suffix) for left removal.
    Group 2 is always empty — the R2 barcode sits entirely to the left of the insert.
    """
    tags = list(reversed(_parse_schema(schema)))
    core = ":".join([*[_R2_TAG[t] for t in tags], _LIGTAG_R2[molecule]])
    suffix = (
        "(?::DpmVar_l(?::DpmVar_r(?::Adap_r)?)?)?"
        if molecule == "dna"
        else "(?::Adap_r)?"
    )
    return f"^((?:Adap_l:)?{core}{suffix})()$"


def r1_trim_regex(schema: str, molecule: MoleculeType) -> str:
    """Return a regex for trimming barcode tags from R1, compatible with trim_reads.

    Group 1 = left tags to remove: Adap_l + DpmVar_l for DNA, Adap_l for RNA.
    Group 2 = right tags to remove: any barcode spillover past the insert.
    """
    tags = _parse_schema(schema)
    tag_names = [_R1_TAG[t] for t in tags]
    ligtag = _LIGTAG_R1[molecule]

    if molecule == "dna":
        right = _capturing_chain(["DpmVar_r", ligtag, *tag_names, "Adap_r"])
        return f"^((?:Adap_l:)?DpmVar_l)(?::{right})?$"
    else:
        right = _capturing_chain([ligtag, *tag_names, "Adap_r"])
        return f"^(Adap_l)?(?::)?{right}?$"


def r2_classification(schema: str) -> dict[str, str]:
    """Return a label→regex mapping for R2 classification given a schema.

    Includes proper DNA/RNA reads and skipped-round variants (pairs of missing
    Odd/Even tags). Labels follow the convention: 'DNA', 'RNA', 'DNA_skip1',
    'RNA_skip1', 'DNA_skip2', 'RNA_skip2', etc.
    """
    middle = "".join(_parse_schema(schema)[:-1])  # strip trailing "T", rejoin to string
    result: dict[str, str] = {
        "DNA": r2_regex(schema, "dna"),
        "RNA": r2_regex(schema, "rna"),
    }
    for k in range(1, len(middle) // 2 + 1):
        reduced = f"L{middle[:-2 * k] if 2 * k < len(middle) else ''}T"
        label = f"_skip{k}"
        result[f"DNA{label}"] = r2_regex(reduced, "dna")
        result[f"RNA{label}"] = r2_regex(reduced, "rna")
    for k in range(1, len(middle) + 2):
        remaining = middle[k - 1:]
        tags = [_R2_TAG["T"]] + [_R2_TAG[c] for c in reversed(remaining)]
        result[f"Eroded_{k}"] = f"^(Adap_l:)?{':'.join(tags)}$"
    return result


def r1_regex(schema: str, molecule: MoleculeType) -> str:
    """Return a regex matching the expected tag_type sequence for R1.

    R1 reads from the insert outward. DpmVar_l is mandatory for DNA, nothing for RNA.
    All tags after the insert are optional and appear in schema order.
    """
    tags = _parse_schema(schema)
    tag_names = [_R1_TAG[t] for t in tags]
    ligtag = _LIGTAG_R1[molecule]

    if molecule == "dna":
        tail = _nested_optional(["DpmVar_r", ligtag, *tag_names, "Adap_r"])
        return f"^(Adap_l:)?DpmVar_l{tail}$"
    else:
        optional = _nested_optional([ligtag, *tag_names, "Adap_r"], first_sep=False)
        return f"^(Adap_l:)?{optional}$"
