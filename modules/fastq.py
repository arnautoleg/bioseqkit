from modules.types_alliases import bounds


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
        total += (ord(ch) - 33)
    return total / len(qual)


def is_gc_within_bounds(seq: str, gc_bounds: bounds) -> bool:
    """
    Check that the GC content (as a percentage of the read length)
    is within the valid range.
    """
    low, high = _normalize_bounds(
        gc_bounds
    )
    gc = gc_content_percent(seq)
    return low <= gc <= high


def is_length_within_bounds(seq: str, length_bounds: bounds) -> bool:
    """
    Check that the sequence length is within the valid range.
    """

    low, high = _normalize_bounds(
        length_bounds
    )
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
