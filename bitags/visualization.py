import io

import polars as pl
import seaborn as sns
from rich.console import Console
from rich.text import Text


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
        text.append(
            seq[pos : pos + len(t_seq)], style=f"on {color_map.get(t_type, 'white')}"
        )
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


def _make_color_map(rows: list[dict]) -> dict[str, str]:
    """Assign divergent palette colors to the unique tag types found in rows."""
    tag_types = sorted(
        {
            t
            for row in rows
            for read in ("r1", "r2")
            for t in (row.get(f"tag_type_{read}") or "").split(":")
            if t
        }
    )
    palette = sns.color_palette("RdBu", len(tag_types))
    return {
        t: "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
        for t, (r, g, b) in zip(tag_types, palette)
    }


def _export_console(console: Console, out: str | None) -> None:
    """Export a recorded console to a file or print to terminal if out is None.

    Supported file formats: .html, .svg (inferred from extension).
    """
    if out is None:
        Console().print(console.export_text())
        return

    ext = out.rsplit(".", 1)[-1].lower()
    if ext == "html":
        content = console.export_html()
    elif ext == "svg":
        content = console.export_svg()
    else:
        raise ValueError(f"Unsupported format '.{ext}'. Use .html or .svg")
    with open(out, "w") as f:
        f.write(content)


def visualize_reads(
    lf: pl.LazyFrame,
    num_rows: int,
    *,
    color_map: dict[str, str] | None = None,
    out: str | None = None,
) -> None:
    """Render a list of read-pair rows to the console with tag-highlighted sequences."""

    def render_read_pair(
        row: dict,
        *,
        console: Console,
        color_map: dict[str, str],
    ) -> None:
        """Render a single read pair to the console with tag-highlighted sequences."""

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
            tag_line.append_text(
                _aligned_tag_names(tag_types, tag_seqs, tag_pos, color_map)
            )

            seq_line = Text()
            seq_line.append(f"{read.upper()} ", style="bold")
            seq_line.append("5' ")
            seq_line.append_text(
                _highlighted_seq(seq, tag_types, tag_seqs, tag_pos, color_map)
            )
            seq_line.append(" 3'")

            console.print(tag_line)
            console.print(seq_line)

    rows = lf.head(num_rows).collect().to_dicts()
    console = Console(record=True, file=io.StringIO())
    color_map = color_map or _make_color_map(rows)

    for row in rows:
        render_read_pair(row, console=console, color_map=color_map)
        console.rule()

    _export_console(console, out)
