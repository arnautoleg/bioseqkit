import os
from modules.types_alliases import bounds, seq_dict


def _normalize_bounds(values: bounds) -> tuple[float, float]:
    """
    Converts bounds to a pair (low, high).
    - Number x → (0, x)
    - Pair (a, b) → (min(a, b), max(a, b))
    """
    if isinstance(values, (int, float)):
        low, high = 0, float(values)
    else:
        if not (isinstance(values, (tuple, list)) and len(values) == 2):
            raise ValueError("bounds must be a number or a 2-tuple")
        limit1, limit2 = values
        low, high = float(min(limit1, limit2)), float(max(limit1, limit2))

    return low, high


def gc_content_percent(seq: str) -> float:
    """
    GC content (as a percentage of the read length).
    Calculated over the full read length: (G + C) / len(seq) * 100.
    Register is ignored.
    """
    if not seq:
        return 0.0
    gc = 0
    for ch in seq.upper():
        if ch == "G" or ch == "C":
            gc += 1
    return 100.0 * gc / len(seq)


def mean_phred33(qual: str) -> float:
    """
    Average quality according to the Phred+33 scale across
    all characters of the qual string.
    """
    if not qual:
        return 0.0
    total = 0
    for ch in qual:
        total += ord(ch) - 33
    return total / len(qual)


def is_gc_within_bounds(seq: str, gc_bounds: bounds) -> bool:
    """
    Check that the GC content (as a percentage of the read length)
    is within the valid range.
    """
    low, high = _normalize_bounds(gc_bounds)
    gc = gc_content_percent(seq)
    return low <= gc <= high


def is_length_within_bounds(seq: str, length_bounds: bounds) -> bool:
    """
    Check that the sequence length is within the valid range.
    """

    low, high = _normalize_bounds(length_bounds)
    length = len(seq)
    return low <= length <= high


def is_quality_pass(qual: str, threshold: float) -> bool:
    """
    Check that the read quality according to the Phred+33
    scale is within the valid range.
    """
    return mean_phred33(qual) >= float(threshold)


def is_valid_pair(seq: str, qual: str) -> bool:
    """
    Check that the lengths of the sequence and the quality string match.
    """
    return len(seq) == len(qual)


def read_fastq_to_dict(input_fastq) -> seq_dict:
    """
    Reads a standard (uncompressed) FASTQ file and returns a dictionary mapping
    read identifiers to their corresponding sequence and quality strings.

    Returns:
    dict[str, tuple[str, str]]: A dictionary where each key is a read ID,
    and each value is a tuple (sequence, quality).
    """
    seqs = {}
    if not os.path.isfile(input_fastq):
        raise FileNotFoundError(f"Файл не найден: {input_fastq}")

    with open(input_fastq, "r", encoding="utf-8", errors="replace") as fin:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not (header and seq and plus and qual):
                break
            if not header.startswith("@") or not plus.startswith("+"):
                continue
            read_id = header[1:].strip()
            seqs[read_id] = (seq.strip(), qual.strip())

    return seqs


def write_fastq_from_dict(seqs, output_fastq) -> str:
    """
    Writes a dictionary {id: (sequence, quality)} to a FASTQ file
    inside the 'filtered' folder. Creates the folder if it does not exist.
    """
    # Force saving into the "filtered" folder
    folder = "filtered"
    if not os.path.exists(folder):
        os.makedirs(folder)

    # Keep only the file name, not any path provided
    file_name = os.path.basename(output_fastq)
    output_path = os.path.join(folder, file_name)

    # Write FASTQ
    with open(output_path, "w", encoding="utf-8") as fout:
        for read_id, (seq, qual) in seqs.items():
            fout.write(f"@{read_id}\n{seq}\n+\n{qual}\n")

    print(f"File saved to: {output_path}")
    return output_path


def filter_fastq_dict(
    seqs: seq_dict,
    gc_bounds: bounds = (0, 100),
    length_bounds: bounds = (0, 2**32),
    quality_threshold: float = 0.0,
) -> seq_dict:
    """
    Filters a dict {id: (sequence, quality)}
    by GC-content, length, and mean Phred+33 quality.

    Arguments:
        seqs: dict[str, tuple[str, str]]
            Dictionary of reads in the form {id: (sequence, quality)}.
        gc_bounds: float | tuple[float, float], default (0, 100)
            Lower and upper bounds for GC-content in percent.
            If a single number is provided, it is treated
            as the upper bound (0, x).
        length_bounds: float | tuple[float, float], default (0, 2**32)
            Lower and upper bounds for read length.
            If a single number is provided, it is treated
            as the upper bound (0, x).
        quality_threshold: float, default 0.0
            Minimum allowed average Phred+33 quality (inclusive).

    Returns:
        dict[str, tuple[str, str]]
            Dictionary containing only reads that passed all filters.

    Raises:
        TypeError:
            If inputs are not in the expected types or formats.
        ValueError:
            If provided bounds are invalid or incorrectly formatted.

    """
    result: seq_dict = {}

    for read_id, pair in seqs.items():
        if not (isinstance(pair, tuple) and len(pair) == 2):
            continue
        seq, qual = pair
        if not is_valid_pair(seq, qual):
            continue

        if (
            is_gc_within_bounds(seq, gc_bounds)
            and is_length_within_bounds(seq, length_bounds)
            and is_quality_pass(qual, quality_threshold)
        ):
            result[read_id] = (seq, qual)
    return result
