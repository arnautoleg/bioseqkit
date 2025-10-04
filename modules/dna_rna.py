DNA_ALPHABET = set("ATCGatcg")
RNA_ALPHABET = set("AUCGaucg")

DNA_COMPLEMENT_MAP = str.maketrans("ATCGatcg", "TAGCtagc")
RNA_COMPLEMENT_MAP = str.maketrans("AUCGaucg", "UAGCuagc")


def _is_subset(seq: str, alphabet: set[str]) -> bool:
    """True if all characters of seq belong to the given alphabet 
    (case-insensitive via mixed-case alphabet)."""
    return set(seq).issubset(alphabet)


def is_nucleic_acid(seq: str) -> bool:
    """True if seq is a valid DNA or RNA (case-insensitive)."""

    set_seq = set(seq)
    return set_seq.issubset(DNA_ALPHABET) or set_seq.issubset(RNA_ALPHABET)


def validate_sequence(seq: str) -> None:
    """
    Validate that seq is DNA or RNA.
    Raises:
        TypeError  – if seq is not a string
        ValueError – if seq contains invalid characters
    """
    if not isinstance(seq, str):
        raise TypeError("Sequence must be a string.")
    if not is_nucleic_acid(seq):
        raise ValueError("Invalid nucleic acid sequence.")


def transcribe(seq: str) -> str:
    """
    Transcribe DNA → RNA (T/t → U/u). If already RNA, return as-is.
    """
    validate_sequence(seq)
    if _is_subset(seq, DNA_ALPHABET):
        return seq.replace("T", "U").replace("t", "u")
    return seq  # already RNA


def reverse(seq: str) -> str:
    """Return the reversed sequence."""
    validate_sequence(seq)
    return seq[::-1]


def complement(seq: str) -> str:
    """"Returns the complementary sequence (DNA or RNA)."""
    validate_sequence(seq)
    if _is_subset(seq, DNA_ALPHABET):
        table = DNA_COMPLEMENT_MAP
    else:
        table = RNA_COMPLEMENT_MAP

    return seq.translate(table)


def reverse_complement(seq: str) -> str:
    """Returns the reverse complementary sequence."""
    return reverse(complement(seq))


def _ensure_strings(items) -> None:
    """Checks that all elements are strings."""
    for it in items:
        if not isinstance(it, str):
            raise TypeError("All sequences must be strings.")
        