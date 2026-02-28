from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar, Iterator, Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


Index = Union[int, slice]
Bounds = Union[int, float, tuple[Union[int, float], Union[int, float]]]


# Task 1: Abstract sequences (OOP)
@dataclass(frozen=True)
class BiologicalSequence(ABC):
    """
    Base abstract class for biological sequences.
    Provides: len, indexing/slicing, iteration, pretty printing.
    Requires: check_alphabet().
    """
    sequence: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "sequence", self.sequence.upper())

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, idx: Index) -> Union[str, "BiologicalSequence"]:
        if isinstance(idx, slice):
            return self.__class__(self.sequence[idx])
        return self.sequence[idx]

    def __iter__(self) -> Iterator[str]:
        return iter(self.sequence)

    def __str__(self) -> str:
        cls = self.__class__.__name__
        seq = self.sequence
        preview = seq if len(seq) <= 60 else seq[:57] + "..."
        return f"{cls}(length={len(self)}, sequence='{preview}')"

    def __repr__(self) -> str:
        return str(self)

    @abstractmethod
    def check_alphabet(self) -> None:
        """Raise ValueError if sequence contains invalid symbols."""
        raise NotImplementedError


class NucleicAcidSequence(BiologicalSequence):
    """
    Base class for DNA/RNA.
    Implements: check_alphabet, complement, reverse, reverse_complement.
    Uses polymorphism via subclass-defined _alphabet and _complement_map.
    """

    _alphabet: ClassVar[Optional[frozenset[str]]] = None
    _complement_map: ClassVar[Optional[dict[str, str]]] = None

    def _ensure_concrete(self) -> None:
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is abstract. Use DNASequence or RNASequence."
            )

    def check_alphabet(self) -> None:
        self._ensure_concrete()

        if self._alphabet is None:
            raise NotImplementedError("Alphabet not defined in subclass.")

        invalid = sorted(set(self.sequence) - self._alphabet)
        if invalid:
            raise ValueError(
                f"Invalid symbols for {self.__class__.__name__}: {invalid}. "
                f"Allowed symbols: {sorted(self._alphabet)}"
            )

    def complement(self) -> "NucleicAcidSequence":
        self._ensure_concrete()
        self.check_alphabet()

        if self._complement_map is None:
            raise NotImplementedError("Complement map not defined in subclass.")

        complemented = "".join(self._complement_map[ch] for ch in self.sequence)
        return self.__class__(complemented)

    def reverse(self) -> "NucleicAcidSequence":
        self._ensure_concrete()
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    _alphabet: ClassVar[frozenset[str]] = frozenset({"A", "T", "G", "C"})
    _complement_map: ClassVar[dict[str, str]] = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self) -> "RNASequence":
        self.check_alphabet()
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    _alphabet: ClassVar[frozenset[str]] = frozenset({"A", "U", "G", "C"})
    _complement_map: ClassVar[dict[str, str]] = {"A": "U", "U": "A", "G": "C", "C": "G"}


class AminoAcidSequence(BiologicalSequence):
    _alphabet: ClassVar[frozenset[str]] = frozenset(set("ACDEFGHIKLMNPQRSTVWY"))

    def check_alphabet(self) -> None:
        invalid = sorted(set(self.sequence) - self._alphabet)
        if invalid:
            raise ValueError(
                f"Invalid amino acids: {invalid}. Allowed symbols: {sorted(self._alphabet)}"
            )

    def fraction_of(self, residues: str) -> float:
        """Fraction of residues (e.g., hydrophobic fraction)."""
        self.check_alphabet()
        rset = set(residues.upper())
        if not rset:
            raise ValueError("Residues must not be empty.")
        return 0.0 if len(self) == 0 else sum(aa in rset for aa in self.sequence) / len(self)


# Optional helper
def is_nucleic_acid(seq: str) -> bool:
    s = seq.upper()
    if "T" in s and "U" in s:
        return False
    try:
        (RNASequence(seq) if "U" in s else DNASequence(seq)).check_alphabet()
        return True
    except ValueError:
        return False



# Task 2: FastQ filtering with Biopython
def _normalize_bounds(values: Bounds) -> tuple[float, float]:
    if isinstance(values, (int, float)):
        return 0.0, float(values)

    if not (isinstance(values, (tuple, list)) and len(values) == 2):
        raise ValueError("Bounds must be a number or a 2-tuple")

    a, b = values
    return float(min(a, b)), float(max(a, b))


def _gc_percent(record: SeqRecord) -> float:
    seq_str = str(record.seq)
    return 0.0 if not seq_str else 100.0 * gc_fraction(seq_str)


def _mean_phred(record: SeqRecord) -> float:
    quals = record.letter_annotations.get("phred_quality")
    return 0.0 if not quals else sum(quals) / len(quals)


def _passes_filters(
    record: SeqRecord,
    gc_bounds: Bounds,
    length_bounds: Bounds,
    quality_threshold: float,
) -> bool:
    gc_low, gc_high = _normalize_bounds(gc_bounds)
    len_low, len_high = _normalize_bounds(length_bounds)

    length = len(record.seq)
    gc = _gc_percent(record)
    q = _mean_phred(record)

    return (gc_low <= gc <= gc_high) and (len_low <= length <= len_high) and (q >= float(quality_threshold))


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Bounds = (0, 100),
    length_bounds: Bounds = (0, 2**32),
    quality_threshold: float = 0.0,
    out_dir: str = "filtered",
) -> str:
    """
    Filter FASTQ by GC%, length, and mean Phred quality using Biopython.
    Writes output into ./filtered (or out_dir). Returns written file path.
    """
    in_path = Path(input_fastq)
    if not in_path.is_file():
        raise FileNotFoundError(f"File not found: {input_fastq}")

    out_folder = Path(out_dir)
    out_folder.mkdir(parents=True, exist_ok=True)
    out_path = out_folder / Path(output_fastq).name

    with in_path.open("rt", encoding="utf-8", errors="replace") as fin, out_path.open("wt", encoding="utf-8") as fout:
        records = SeqIO.parse(fin, "fastq")
        passed = (
            rec
            for rec in records
            if _passes_filters(rec, gc_bounds=gc_bounds, length_bounds=length_bounds, quality_threshold=quality_threshold)
        )
        n_written = SeqIO.write(passed, fout, "fastq")

    print(f"Written {n_written} reads -> {out_path}")
    return str(out_path)


# Demo/tests
if __name__ == "__main__":
    # FASTQ demo
    out = filter_fastq(
        input_fastq="example_data/example_fastq.fastq",
        output_fastq="example_filtered.fastq",
        gc_bounds=(40, 60),
        length_bounds=(50, 300),
        quality_threshold=20,
    )
    print("FASTQ output:", out)

    # Sequence demos (print plain sequences like in your examples)
    print(DNASequence("ATG").transcribe().sequence)          # AUG
    print(DNASequence("ATG").reverse().sequence)             # GTA
    print(DNASequence("AtG").complement().sequence)          # TAC
    print(DNASequence("ATc").reverse_complement().sequence)  # GAT
    print(is_nucleic_acid("TTUU"))                           # False