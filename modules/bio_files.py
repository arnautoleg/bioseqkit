import re
import os
from typing import Iterable
from pathlib import Path


def _default_out_path(input_fasta: str) -> str:
    """
    Build a default output path for a one-line-per-record FASTA.

    If the input filename has an extension, insert ".oneline" before it.
    Otherwise, append ".oneline.fasta".

    Examples:
        "/data/genome.fasta"   -> "/data/genome.oneline.fasta"
        "reads.fa"             -> "reads.oneline.fa"
        "/tmp/contigs"         -> "/tmp/contigs.oneline.fasta"

    Args:
        input_fasta: Path to the input FASTA file.

    Returns:
        The suggested output file path with the ".oneline" suffix.
    """
    root, ext = os.path.splitext(input_fasta)
    if ext:
        return f"{root}.oneline{ext}"
    return f"{input_fasta}.oneline.fasta"


def flush_sequence(seq_buf: list[str], fout, *, final: bool = False) -> None:
    """
    Write the accumulated sequence as a single line, then clear the buffer.

    Joins all strings in `seq_buf`, removes spaces and tabs, and writes the
    result to `fout`. If `final` is False (default), a trailing newline
    is added.
    If `final` is True, no newline is appended—useful for avoiding a newline
    at the very end of the FASTA file. The buffer is cleared regardless.

    Notes:
      - If `seq_buf` is empty or the joined string is empty after stripping
        spaces/tabs, nothing is written.
      - Only spaces (' ') and tabs ('\\t') are removed; other characters
      are left as-is.
      - `fout` must be a text-like object with a `.write(str) -> int` method.

    Args:
        seq_buf: List of sequence fragments (e.g., lines) to be joined and written.
        fout:    An open text file-like object to write into.
        final:   When True, do not append a trailing newline.
    """
    if not seq_buf:
        return
    seq_line = "".join(seq_buf).replace(" ", "").replace("\t", "")
    if seq_line:
        if final:
            fout.write(seq_line)  # no newline at the very end of the file
        else:
            fout.write(seq_line + "\n")
    seq_buf.clear()


_rx_num = re.compile(r"(\d+)")


def _normalize_targets(genes: str | Iterable[str]) -> set[str]:
    """
    Normalize the user-specified target genes to a lowercase, deduplicated set.

    Args:
        genes: Gene name (str) or a collection of gene names (Iterable[str]).

    Returns:
        A non-empty set of normalized gene names.

    Examples:
        >>> _normalize_targets("  GltA ")
        {'glta'}
        >>> _normalize_targets(["gltA", "  SDHA", " "])
        {'glta', 'sdha'}
    """
    items = [genes] if isinstance(genes, str) else list(genes)
    targets = {str(g).casefold().strip() for g in items if g and str(g).strip()}
    if not targets:
        raise ValueError(
            "Список генов интереса пуст или содержит только пустые значения"
        )
    return targets


def _parse_location(loc_text: str) -> tuple[int, int, int]:
    """Parse GenBank location into (start, end, strand).
    Supports complement/join/fuzzy; strand ∈ {+1, -1}."""
    txt = loc_text.strip()
    strand = -1 if "complement" in txt else 1
    nums = [int(m.group(1)) for m in _rx_num.finditer(txt)]
    if not nums:
        return (1, 1, strand)
    return (min(nums), max(nums), strand)


def _matches_anchor(qual: dict[str, list[str]], targets_norm: set[str]) -> bool:
    """Return True if feature qualifiers match targets
    by gene/locus_tag/protein_id (exact) or product (substring)."""
    for key in ("gene", "locus_tag", "protein_id"):
        values = qual.get(key, [])
        if values and values[0].casefold() in targets_norm:
            return True
    product = (qual.get("product", [""])[0]).casefold()
    if product:
        for t in targets_norm:
            if t and t in product:
                return True
    return False


def _feature_header(
    rec_id: str, start: int, end: int, strand: int, qual: dict[str, list[str]]
) -> str:
    """Build a FASTA header with record id, coordinates,
    strand, and key qualifiers."""
    strand_chr = "+" if strand >= 0 else "-"
    parts = [f"rec={rec_id}", f"pos={start}..{end}", f"strand={strand_chr}"]
    for key in ("gene", "locus_tag", "protein_id", "product"):
        vals = qual.get(key, [])
        if vals and vals[0]:
            parts.append(f"{key}={vals[0]}")
    return "|".join(parts)


def _wrap60(text: str) -> list[str]:
    """Split a string into lines of length 60 (last line may be shorter)."""
    return [text[i : i + 60] for i in range(0, len(text), 60)]


def _flush_current_feature(
    cds_list: list[dict], current_feat: dict | None
) -> dict | None:
    """
    If current_feat exists, append it to cds_list and return None (cleared).
    Otherwise, return current_feat unchanged.
    """
    if current_feat is not None:
        cds_list.append(current_feat)
        return None
    return current_feat


def _finish_record(
    rec_id: str | None,
    cds_list: list[dict],
    targets_norm: set[str],
    n_before: int,
    n_after: int,
    seen_keys: set[tuple[str, int, int, int]],
    fout,
) -> int:
    """
    Select neighbors around anchors and write their translations to FASTA.
    Excludes anchors themselves. Clears cds_list before returning.
    Returns the number of sequences written for this record.
    """
    if not cds_list:
        return 0

    anchors = [
        i
        for i, feat in enumerate(cds_list)
        if _matches_anchor(feat["qual"], targets_norm)
    ]
    if not anchors:
        cds_list.clear()
        return 0

    selected_idx: set[int] = set()
    for idx in anchors:
        if n_before > 0:
            a = max(0, idx - n_before)
            selected_idx.update(range(a, idx))
        if n_after > 0:
            b = min(len(cds_list), idx + 1 + n_after)
            selected_idx.update(range(idx + 1, b))

    selected_idx.difference_update(anchors)

    written_here = 0
    rec = rec_id or ""
    for j in sorted(selected_idx):
        feat = cds_list[j]
        qual = feat["qual"]
        translations = qual.get("translation", [])
        if not translations:
            continue
        aa = translations[0].replace(" ", "").replace("\n", "")
        if not aa:
            continue

        key = (rec, feat["start"], feat["end"], feat["strand"])
        if key in seen_keys:
            continue
        seen_keys.add(key)

        header = _feature_header(rec, feat["start"], feat["end"], feat["strand"], qual)
        fout.write(f">{header}\n")
        for chunk in _wrap60(aa):
            fout.write(chunk + "\n")
        written_here += 1

    cds_list.clear()
    return written_here
