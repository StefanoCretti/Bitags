from rich.console import Console
from rich.text import Text

_TAG_COLORS: dict[str, str] = {
    "Adap_l": "#555555",
    "Adap_r": "#555555",
    "DpmVar_l": "#4488CC",
    "DpmVar_r": "#4488CC",
    "Dpm_r1": "#44AA44",
    "Dpm_r2": "#44AA44",
    "Rpm_r1": "#AAAA00",
    "Rpm_r2": "#AAAA00",
    "Odd_r1": "#CC4444",
    "Odd_r2": "#CC4444",
    "Even_r1": "#AA44AA",
    "Even_r2": "#AA44AA",
    "Term_r1": "#44AAAA",
    "Term_r2": "#44AAAA",
}


def _highlighted_seq(
    seq: str,
    tag_types: list[str],
    tag_seqs: list[str],
    tag_pos: list[int],
    color_map: dict[str, str],
) -> Text:
    """Build a Rich Text with tag regions colored by type."""
    text = Text()
    cursor = 0
    for t_type, t_seq, pos in zip(tag_types, tag_seqs, tag_pos):
        if pos > cursor:
            text.append(seq[cursor:pos])
        text.append(seq[pos : pos + len(t_seq)], style=f"on {color_map.get(t_type, 'white')}")
        cursor = pos + len(t_seq)
    text.append(seq[cursor:])
    return text


def _aligned_tag_names(
    tag_types: list[str],
    tag_seqs: list[str],
    tag_pos: list[int],
    color_map: dict[str, str],
) -> Text:
    """Build a Rich Text with tag names positioned to align with their sequence regions."""
    text = Text()
    cursor = 0
    for t_type, t_seq, pos in zip(tag_types, tag_seqs, tag_pos):
        if pos > cursor:
            text.append(" " * (pos - cursor))
        width = len(t_seq)
        text.append(t_type[:width].center(width), style=color_map.get(t_type, "white"))
        cursor = pos + width
    return text


def render_read_pair(
    row: dict,
    *,
    console: Console,
    color_map: dict[str, str] | None = None,
) -> None:
    """Render a single read pair to the console with tag-highlighted sequences."""
    color_map = color_map or _TAG_COLORS

    console.print(row["name_r1"])

    for read in ("r1", "r2"):
        tag_type = row.get(f"tag_type_{read}") or ""
        tag_seq = row.get(f"tag_seq_{read}") or ""
        tag_pos = row.get(f"tag_pos_{read}") or []
        seq = row.get(f"sequence_{read}") or ""

        tag_types = tag_type.split(":") if tag_type else []
        tag_seqs = tag_seq.split(":") if tag_seq else []

        prefix = f"{read.upper()} 5' "

        tag_line = Text(" " * len(prefix))
        tag_line.append_text(_aligned_tag_names(tag_types, tag_seqs, tag_pos, color_map))

        seq_line = Text()
        seq_line.append(f"{read.upper()} ", style="bold")
        seq_line.append("5' ")
        seq_line.append_text(_highlighted_seq(seq, tag_types, tag_seqs, tag_pos, color_map))
        seq_line.append(" 3'")

        console.print(tag_line)
        console.print(seq_line)
