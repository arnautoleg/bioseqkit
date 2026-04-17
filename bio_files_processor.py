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


def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str | None = None) -> str:
    """
    Convert a standard (non-compressed) FASTA so that each record's
    sequence is a single line.

    For every header line starting with '>', all subsequent non-empty sequence
    lines are concatenated (whitespace within sequences is removed) and written
    as one line.
    A trailing newline is NOT added after the final sequence in the output
    file.

    Args:
        input_fasta: Path to the input FASTA file (text, non-compressed).
        output_fasta: Optional path for the resulting FASTA. If omitted,
            a path is auto-derived by inserting ".oneline" before the original
            extension (or appending ".oneline.fasta" if none).

    Returns:
        The path to the written one-line-per-record FASTA.

    Raises:
        FileNotFoundError: If `input_fasta` does not exist.
    """
    if output_fasta is None:
        output_fasta = _default_out_path(input_fasta)

    out_dir = os.path.dirname(os.path.abspath(output_fasta))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    with open(input_fasta, "r", encoding="utf-8") as fin, \
         open(output_fasta, "w", encoding="utf-8") as fout:

        seq_buf: list[str] = []
        seen_header = False

        for raw in fin:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header:
                    flush_sequence(seq_buf, fout)
                fout.write(line + "\n")
                seen_header = True
            else:
                seq_buf.append(line)

        flush_sequence(seq_buf, fout, final=True)

    return output_fasta


def parse_blast_best(input_file: str, output_file: str) -> list[str]:
    """
    Extract the single “best hit” (top Description line) from each
    “Sequences producing significant alignments:” block of a BLAST
    plain-text report and save unique hits to a one-column file.

    Parsing rules:
      - Locate each block starting with the line:
        "Sequences producing significant alignments:"
      - Skip the column header line(s) and any ellipsis-only lines.
      - Take ONLY the first non-header Description line from that block.
      - Trim trailing dots/ellipsis and extra columns.
      - Collect unique descriptions across the whole file.
      - Sort case-insensitively once at the end and write to `output_file`.

    Args:
        input_file: Path to a BLAST TXT report
        (e.g., output from blastp/blastn -outfmt 0).
        output_file: Path to a text file that will contain
        one description per line.

    Returns:
        Sorted list of unique best-hit descriptions.

    Example:
        >>> hits = parse_blast_best("blast_report.txt", "best_hits.txt")
        >>> hits[:3]
        ['Escherichia coli K12 protein X',
         'Hypothetical protein ABC',
         'Aldolase class I']
    """
    header_words = (
        "Description", "Accession", "Max score", "Total score",
        "Score", "Query cover", "Value", "Ident", "Scientific",
        "Common", "Acc."
    )

    best_hits: set[str] = set()

    in_hit_block = False
    took_top_hit = False

    with open(input_file, "rt", encoding="utf-8", errors="replace") as fin:
        for line in fin:
            stripped_line = line.rstrip("\r\n")

            # Enter the hit-list block
            if not in_hit_block:
                if "Sequences producing significant alignments:" in stripped_line:
                    in_hit_block = True
                    took_top_hit = False
                continue

            # Exit conditions for the block
            if not stripped_line.strip() or "Alignments:" in stripped_line:
                in_hit_block = False
                continue

            if any(header_word in stripped_line for header_word in header_words):
                continue
            if stripped_line.strip().startswith("..."):
                continue

            if not took_top_hit:
                match = re.match(r'^(?P<desc>.*?)\s+\.{3,}', stripped_line)
                if match:
                    description = match.group("desc").strip()
                else:
                    description = re.sub(r'\s{2,}\S.*$', '', stripped_line).strip()

                description = re.sub(r'\s*\.\.\.\s*$', '', description)
                description = re.sub(r'^(.*\])\s+[A-Za-z].*$', r'\1', description)

                if description:
                    best_hits.add(description)
                took_top_hit = True

    unique_sorted_hits = sorted(best_hits, key=lambda x: x.lower())

    with open(output_file, "wt", encoding="utf-8", newline="\n") as fout:
        for description in unique_sorted_hits:
            fout.write(description + "\n")

    return unique_sorted_hits


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: str | Iterable[str],
    output_fasta: str,
    n_before: int = 1,
    n_after: int = 1,
) -> int:
    """
    Extract neighbor CDS protein sequences around target genes from a GenBank file (no Biopython)
    and write them to a FASTA file. Anchor CDS (the targets themselves) are excluded.

    For each record:
      1) Parse FEATURES → CDS in file order.
      2) Mark CDS as an anchor if its qualifiers match `genes` by:
         - exact (case-insensitive) match in gene/locus_tag/protein_id, or
         - substring (case-insensitive) match in product.
      3) Select up to `n_before` CDS immediately before and `n_after` CDS immediately after
         each anchor (indices relative to the CDS list within the same record).
      4) Write only neighbors that have a `translation` qualifier, wrapping sequences to 60 chars/line.
         The FASTA header includes rec id, coordinates, strand, and key qualifiers.

    Args:
        input_gbk: Path to the input GenBank (.gb/.gbk) file (text).
        genes:     Target genes (string or iterable of strings) to anchor by.
        output_fasta: Path to the output FASTA file to write neighbors into.
        n_before:  Number of CDS to include before each anchor (>= 0). Default: 1.
        n_after:   Number of CDS to include after each anchor (>= 0).  Default: 1.

    Returns:
        The number of FASTA records written (unique neighbors; anchors excluded).

    """
    targets_norm = _normalize_targets(genes)

    written_total = 0
    seen_keys: set[tuple[str, int, int, int]] = set()

    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)

    with open(input_gbk, "rt", encoding="utf-8", errors="replace") as fin, \
         open(output_fasta, "wt", encoding="utf-8", newline="\n") as fout:

        rec_id: str | None = None
        in_features = False

        cds_list: list[dict] = []
        current_feat: dict | None = None
        collecting_multiline = False
        current_qual_key = ""

        for raw in fin:
            line = raw.rstrip("\n")

            if line.startswith("LOCUS"):
                if rec_id is not None:
                    current_feat = _flush_current_feature(cds_list, current_feat)
                    written_total += _finish_record(
                        rec_id, cds_list, targets_norm, n_before, n_after, seen_keys, fout
                    )
                    in_features = False
                parts = line.split()
                rec_id = parts[1] if len(parts) > 1 else ""
                current_feat = None
                collecting_multiline = False
                current_qual_key = ""
                cds_list.clear()
                continue

            if line.strip() == "//":
                current_feat = _flush_current_feature(cds_list, current_feat)
                written_total += _finish_record(
                    rec_id, cds_list, targets_norm, n_before, n_after, seen_keys, fout
                )
                in_features = False
                rec_id = None
                current_feat = None
                collecting_multiline = False
                current_qual_key = ""
                cds_list.clear()
                continue

            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN"):
                current_feat = _flush_current_feature(cds_list, current_feat)
                in_features = False
                collecting_multiline = False
                current_qual_key = ""
                continue

            if not in_features:
                continue

            if line.startswith("     "):
                body = line[5:]
                if body and not body.startswith(" "):
                    current_feat = _flush_current_feature(cds_list, current_feat)
                    current_feat = None
                    collecting_multiline = False
                    current_qual_key = ""

                    tokens = body.split()
                    if not tokens:
                        continue
                    ftype = tokens[0]
                    if ftype != "CDS":
                        continue

                    loc_text = body[len(ftype):].strip()
                    start, end, strand = _parse_location(loc_text)
                    current_feat = {"start": start, "end": end, "strand": strand, "qual": {}}
                    continue

                if current_feat is not None:
                    text = body.strip()
                    if text.startswith("/"):
                        collecting_multiline = False
                        current_qual_key = ""
                        if "=" in text:
                            key, val = text[1:].split("=", 1)
                            key = key.strip()
                            val = val.strip()
                            if val.startswith('"') and not val.endswith('"'):
                                val = val[1:]
                                current_feat["qual"].setdefault(key, []).append(val)
                                collecting_multiline = True
                                current_qual_key = key
                            else:
                                if val.startswith('"') and val.endswith('"'):
                                    val = val[1:-1]
                                current_feat["qual"].setdefault(key, []).append(val)
                        else:
                            key = text[1:].strip()
                            current_feat["qual"].setdefault(key, []).append("")
                    else:
                        if collecting_multiline and current_qual_key:
                            cont = text
                            if cont.endswith('"'):
                                cont = cont[:-1]
                                collecting_multiline = False
                                current_feat["qual"][current_qual_key][-1] += cont
                                current_qual_key = ""
                            else:
                                current_feat["qual"][current_qual_key][-1] += cont

        if rec_id is not None:
            current_feat = _flush_current_feature(cds_list, current_feat)
            written_total += _finish_record(
                rec_id, cds_list, targets_norm, n_before, n_after, seen_keys, fout
            )

    return None
