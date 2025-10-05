from modules.dna_rna import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement,
    _ensure_strings
)

from modules.fastq import (
    is_valid_pair,
    is_gc_within_bounds,
    is_length_within_bounds,
    is_quality_pass,
)


from modules.types_alliases import seq_dict, bounds


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


def filter_fastq(seqs: seq_dict,
                 gc_bounds: bounds = (0, 100),
                 length_bounds: bounds = (0, 2**32),
                 quality_threshold: float = 0.0) -> seq_dict:
    """
    Filters FASTQ reads by GC-content, sequence length,
    and average quality (Phred+33).

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

