from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioseqkit import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    _normalize_bounds,
    filter_fastq,
    is_nucleic_acid,
)
from bio_files_processor import (
    _normalize_targets,
    convert_multiline_fasta_to_oneline,
    parse_blast_best,
)


def _write_fastq(path: Path, records: list[SeqRecord]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "fastq")


def test_dna_transcribe_returns_rna_sequence() -> None:
    result = DNASequence("ATGC").transcribe()
    assert isinstance(result, RNASequence)
    assert result.sequence == "AUGC"


def test_reverse_complement_dna() -> None:
    result = DNASequence("ATGC").reverse_complement()
    assert result.sequence == "GCAT"


def test_amino_acid_fraction_of() -> None:
    seq = AminoAcidSequence("AVILGG")
    assert seq.fraction_of("AVL") == pytest.approx(3 / 6)


def test_invalid_dna_alphabet_raises_value_error() -> None:
    with pytest.raises(ValueError):
        DNASequence("ATBX").check_alphabet()


def test_is_nucleic_acid_false_when_t_and_u_mixed() -> None:
    assert is_nucleic_acid("TTUU") is False


def test_normalize_bounds_single_number() -> None:
    assert _normalize_bounds(10) == (0.0, 10.0)


def test_convert_multiline_fasta_to_oneline_writes_file(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    input_fasta.write_text(">seq1\nATG\nCC\n>seq2\nAAA\nTTT\n", encoding="utf-8")

    written_path = convert_multiline_fasta_to_oneline(str(input_fasta), str(output_fasta))

    assert written_path == str(output_fasta)
    assert output_fasta.read_text(encoding="utf-8") == ">seq1\nATGCC\n>seq2\nAAATTT"


def test_filter_fastq_filters_by_quality_and_creates_output(tmp_path: Path) -> None:
    input_fastq = tmp_path / "reads.fastq"
    out_dir = tmp_path / "filtered"

    good = SeqRecord(Seq("GCGC"), id="good", description="")
    good.letter_annotations["phred_quality"] = [40, 40, 40, 40]

    bad = SeqRecord(Seq("ATAT"), id="bad", description="")
    bad.letter_annotations["phred_quality"] = [0, 0, 0, 0]

    _write_fastq(input_fastq, [good, bad])

    output_path = filter_fastq(
        input_fastq=str(input_fastq),
        output_fastq="filtered.fastq",
        gc_bounds=(0, 100),
        length_bounds=(0, 10),
        quality_threshold=20,
        out_dir=str(out_dir),
    )

    written_records = list(SeqIO.parse(output_path, "fastq"))
    assert Path(output_path).exists()
    assert len(written_records) == 1
    assert written_records[0].id == "good"


def test_parse_blast_best_extracts_unique_top_hits(tmp_path: Path) -> None:
    blast_text = (
        "Sequences producing significant alignments:\n"
        "Description    Score\n"
        "Protein kinase ABC   ...   100\n"
        "Other hit            ...   90\n"
        "\n"
        "Sequences producing significant alignments:\n"
        "Description    Score\n"
        "Protein kinase ABC   ...   100\n"
        "\n"
    )
    input_file = tmp_path / "blast.txt"
    output_file = tmp_path / "best_hits.txt"
    input_file.write_text(blast_text, encoding="utf-8")

    hits = parse_blast_best(str(input_file), str(output_file))

    assert hits == ["Protein kinase ABC"]
    assert output_file.read_text(encoding="utf-8") == "Protein kinase ABC\n"


def test_normalize_targets_empty_raises_value_error() -> None:
    with pytest.raises(ValueError):
        _normalize_targets([" ", ""])
