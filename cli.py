from __future__ import annotations

import argparse
from typing import Sequence

from bioseqkit import filter_fastq
from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    parse_blast_best,
    select_genes_from_gbk_to_fasta,
)


import logging


logging.basicConfig(
    filename="bioseqkit.log",
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)


def _parse_pair(raw: Sequence[str], cast=float) -> tuple[float, float]:
    if len(raw) != 2:
        raise argparse.ArgumentTypeError("Expected exactly two values.")
    try:
        left = cast(raw[0])
        right = cast(raw[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Both values must be numeric.") from exc
    return left, right


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioseqkit",
        description="CLI for sequence utilities, FASTQ filtering, and file processing.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fastq_parser = subparsers.add_parser(
        "filter-fastq",
        help="Filter FASTQ reads by GC percent, length, and mean quality.",
    )
    fastq_parser.add_argument("--input-fastq", required=True, help="Path to input FASTQ.")
    fastq_parser.add_argument("--output-fastq", required=True, help="Name of output FASTQ file.")
    fastq_parser.add_argument(
        "--gc-bounds",
        nargs=2,
        metavar=("MIN", "MAX"),
        default=(0.0, 100.0),
        help="Lower and upper GC%% bounds.",
    )
    fastq_parser.add_argument(
        "--length-bounds",
        nargs=2,
        metavar=("MIN", "MAX"),
        default=(0, 2**32),
        help="Lower and upper read-length bounds.",
    )
    fastq_parser.add_argument(
        "--quality-threshold",
        type=float,
        default=0.0,
        help="Minimum mean Phred quality.",
    )
    fastq_parser.add_argument(
        "--out-dir",
        default="filtered",
        help="Directory for filtered FASTQ output.",
    )

    fasta_parser = subparsers.add_parser(
        "fasta-oneline",
        help="Convert a multiline FASTA into one-line-per-record FASTA.",
    )
    fasta_parser.add_argument("--input-fasta", required=True, help="Path to input FASTA.")
    fasta_parser.add_argument("--output-fasta", help="Optional path to output FASTA.")

    blast_parser = subparsers.add_parser(
        "parse-blast",
        help="Extract the best hits from a BLAST plain-text report.",
    )
    blast_parser.add_argument("--input-file", required=True, help="Path to BLAST report.")
    blast_parser.add_argument("--output-file", required=True, help="Path to output text file.")

    gbk_parser = subparsers.add_parser(
        "gbk-neighbors",
        help="Extract translated CDS neighbors around target genes from GenBank.",
    )
    gbk_parser.add_argument("--input-gbk", required=True, help="Path to input GenBank file.")
    gbk_parser.add_argument(
        "--genes",
        nargs="+",
        required=True,
        help="Target gene names (space-separated).",
    )
    gbk_parser.add_argument("--output-fasta", required=True, help="Path to output FASTA.")
    gbk_parser.add_argument("--n-before", type=int, default=1, help="Number of CDS before anchor.")
    gbk_parser.add_argument("--n-after", type=int, default=1, help="Number of CDS after anchor.")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.command == "fasta-oneline":
            logger.info(f"Running fasta-oneline on {args.input_fasta}")

            output_path = convert_multiline_fasta_to_oneline(
                args.input_fasta,
                args.output_fasta,
            )

            logger.info(f"Output written to {output_path}")
            print(output_path)

        elif args.command == "parse-blast":
            logger.info(f"Parsing BLAST file {args.input_file}")

            hits = parse_blast_best(
                args.input_file,
                args.output_file,
            )

            logger.info(f"Extracted {len(hits)} hits")
            print(f"Extracted {len(hits)} unique best hits.")

        elif args.command == "filter-fastq":
            logger.info(f"Filtering FASTQ: {args.input_fastq}")

            gc_bounds = _parse_pair(args.gc_bounds, cast=float)
            length_bounds = _parse_pair(args.length_bounds, cast=int)

            output = filter_fastq(
                input_fastq=args.input_fastq,
                output_fastq=args.output_fastq,
                gc_bounds=gc_bounds,
                length_bounds=length_bounds,
                quality_threshold=args.quality_threshold,
                out_dir=args.out_dir,
            )

            logger.info(f"Filtered FASTQ saved to {output}")
            print(output)

        elif args.command == "gbk-neighbors":
            logger.info(f"Processing GBK: {args.input_gbk}")

            n = select_genes_from_gbk_to_fasta(
                input_gbk=args.input_gbk,
                genes=args.genes,
                output_fasta=args.output_fasta,
                n_before=args.n_before,
                n_after=args.n_after,
            )

            logger.info(f"Written {n} FASTA records")
            print(f"Written {n} FASTA records.")

        return 0

    except Exception as e:
        logger.exception("Error occurred during execution")
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
