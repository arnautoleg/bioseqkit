from modules.dna_rna import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement,
    _ensure_strings,
)

from modules.fastq import read_fastq_to_dict, write_fastq_from_dict, filter_fastq_dict


from modules.types_alliases import bounds


def run_dna_rna_tools(*args):
    """
    Runs the specified DNA/RNA procedure on one or more input sequences.

    Arguments:
        *args:
            One or more strings followed by a procedure name (str).
            The last argument must be one of:
              - "is_nucleic_acid"
              - "transcribe"
              - "reverse"
              - "complement"
              - "reverse_complement"

    Returns:
        str | bool | list[str] | list[bool]
        Single value for one input sequence; list of
        results for multiple inputs.

    Raises:
        ValueError:
            If fewer than two arguments are provided or
            if the procedure name is unknown.
        TypeError:
            If the last argument is not a string or if
            any of the sequences are not strings.

    Notes:
        Prints the result to the console before returning it.
    """

    if len(args) < 2:
        raise ValueError("Provide at least one sequence and a procedure")

    *sequences, procedure = args

    if not isinstance(procedure, str):
        raise TypeError("The last argument is a name of the procedure (str).")

    _ensure_strings(sequences)

    procedures = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    func = procedures.get(procedure)
    if func is None:
        raise ValueError(f"Unknown procedure: {procedure}")

    results = [func(seq) for seq in sequences]

    to_print = results[0] if len(results) == 1 else results
    print(to_print)

    return to_print


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: bounds = (0, 100),
    length_bounds: bounds = (0, 2**32),
    quality_threshold: float = 0.0,
) -> str:
    """
    Reads FASTQ from `input_fastq`, filters reads by GC%, length,
    and mean Phred+33 quality, then writes the passed reads to `output_fastq`
    (saved into ./filtered/ by your writer).
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
        str: Path to the written FASTQ (as returned by write_fastq_from_dict).
        In my opinion is comfortable for users
    """

    seqs = read_fastq_to_dict(input_fastq)
    filtered = filter_fastq_dict(
        seqs,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )
    out_path = write_fastq_from_dict(filtered, output_fastq)
    return out_path
