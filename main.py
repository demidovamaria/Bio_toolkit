from typing import Dict, Tuple, Union
import os
import argparse
import logging

from bio_utils.bio_files_processor import convert_multiline_fasta_to_oneline
from bio_utils.parse_blast_output import parse_blast_output
from abc import ABC, abstractmethod
from typing import Iterator
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.
    Fixes the common interface.
    """

    def __init__(self, sequence: str):
        self._sequence = sequence.upper()
        self._check_alphabet()

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, item):
        """
        Supports indexing and slicing.
        """
        return self.__class__(self._sequence[item])

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequence)

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({self._sequence})"

    @abstractmethod
    def _check_alphabet(self) -> None:
        """
        Checks that the sequence alphabet is valid.
        """
        pass

class NucleicAcidSequence(BiologicalSequence, ABC):
    """
    Abstract class for DNA and RNA sequences.
    """

    _alphabet: set[str]
    _complement_map: dict[str, str]

    def complement(self):
        """
        Returns complementary sequence.
        """
        if not hasattr(self, "_complement_map"):
            raise NotImplementedError

        return self.__class__(
            "".join(self._complement_map[n] for n in self._sequence)
        )

    def reverse(self):
        """
        Returns reversed sequence.
        """
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self):
        """
        Returns reverse-complement sequence.
        """
        return self.complement().reverse()

    def _check_alphabet(self) -> None:
        invalid = set(self._sequence) - self._alphabet
        if invalid:
            raise ValueError(f"Invalid characters: {invalid}")
        
class DNASequence(NucleicAcidSequence):
    """
    DNA sequence.
    """

    _alphabet = {"A", "T", "G", "C"}
    _complement_map = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }

    def transcribe(self):
        """
        Transcribes DNA to RNA.
        """
        return RNASequence(self._sequence.replace("T", "U"))
    
class RNASequence(NucleicAcidSequence):
    """
    RNA sequence.
    """

    _alphabet = {"A", "U", "G", "C"}
    _complement_map = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
    }

class AminoAcidSequence(BiologicalSequence):
    """
    Amino acid (protein) sequence.
    """

    _alphabet = set("ACDEFGHIKLMNPQRSTVWY")

    def molecular_weight(self) -> int:
        """
        Returns approximate molecular weight.
        """
        return len(self._sequence) * 110

    def _check_alphabet(self) -> None:
        invalid = set(self._sequence) - self._alphabet
        if invalid:
            raise ValueError(f"Invalid amino acids: {invalid}")
        
def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] = (0.0, 1.0),
    length_bounds: tuple[int, int] = (0, 2**32),
    quality_threshold: float = 0.0,
) -> None:
    """
    Filters FASTQ reads using Biopython.

    Parameters:
        input_fastq: path to input FASTQ file
        output_fastq: path to output FASTQ file
        gc_bounds: allowed GC fraction range (0–1)
        length_bounds: allowed read length range
        quality_threshold: minimal mean Phred quality
    """

    try:
        passed_records = []

        for record in SeqIO.parse(input_fastq, "fastq"):
            seq_length = len(record.seq)
            gc = gc_fraction(record.seq)
            qualities = record.letter_annotations["phred_quality"]
            mean_quality = sum(qualities) / len(qualities)

            if not (length_bounds[0] <= seq_length <= length_bounds[1]):
                continue
            if not (gc_bounds[0] <= gc <= gc_bounds[1]):
                continue
            if mean_quality < quality_threshold:
                continue

            passed_records.append(record)

        logging.info(f"Total reads passed filter: {len(passed_records)}")

        SeqIO.write(passed_records, output_fastq, "fastq")

    except Exception as e:
        logging.error(f"Error in filter_fastq: {str(e)}")
        raise
    
    logging.info(f"Total reads passed filter: {len(passed_records)}")
    SeqIO.write(passed_records, output_fastq, "fastq")
    

# filter_fastq(
#     input_fastq="data/example_fastq.fastq",
#     output_fastq="filtered_out.fastq",
#     gc_bounds=(0.15, 1.0),
#     length_bounds=(10, 2**32),
# )

# convert_multiline_fasta_to_oneline(
#     input_fasta='data/example_multiline_fasta.fasta',
# )

# parse_blast_output(
#     input_file='data/example_blast_results.txt',
#     output_file='parsed_blast_result.txt'
# )


def main():
    """
    Entry point for the Bio Utils command-line interface (CLI).

    This function parses command-line arguments using argparse and dispatches
    execution to one of the supported subcommands.

    Supported subcommands:

    1. fastq
        Filters FASTQ reads using the following criteria:
        - GC content range
        - Sequence length range
        - Mean Phred quality score

        Arguments:
            --input     (str) Path to input FASTQ file (required)
            --output    (str) Path to output FASTQ file (required)
            --min_gc    (float) Minimum GC fraction (default: 0.0)
            --max_gc    (float) Maximum GC fraction (default: 1.0)
            --min_len   (int) Minimum sequence length (default: 0)
            --max_len   (int) Maximum sequence length (default: large number)
            --quality   (float) Minimum mean quality (default: 0.0)

    2. fasta
        Converts a multiline FASTA file into a single-line FASTA format.

        Arguments:
            --input     (str) Path to input FASTA file (required)
            --output    (str) Path to output FASTA file (optional)

    3. blast
        Parses BLAST output and extracts unique top-hit protein names,
        sorted alphabetically.

        Arguments:
            --input     (str) Path to input BLAST file (required)
            --output    (str) Path to output file (required)

    Logging:
        - Logs the selected command (INFO level)
        - Logs execution steps for each subcommand
        - Logs errors if execution fails (ERROR level)

    Raises:
        Exception: Re-raises any exception after logging it.

    Notes:
        This function is intended to be used as a script entry point and should
        be executed via the command line, not imported and called directly.
    """

    logging.basicConfig(
    filename="logs/bio_tool.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(description="Bioinformatics Toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # FASTQ
    fastq_parser = subparsers.add_parser("fastq")
    fastq_parser.add_argument("--input", required=True)
    fastq_parser.add_argument("--output", required=True)
    fastq_parser.add_argument("--min_gc", type=float, default=0.0)
    fastq_parser.add_argument("--max_gc", type=float, default=1.0)
    fastq_parser.add_argument("--min_len", type=int, default=0)
    fastq_parser.add_argument("--max_len", type=int, default=2**32)
    fastq_parser.add_argument("--quality", type=float, default=0.0)

    # FASTA
    fasta_parser = subparsers.add_parser("fasta")
    fasta_parser.add_argument("--input", required=True)
    fasta_parser.add_argument("--output")

    # BLAST
    blast_parser = subparsers.add_parser("blast")
    blast_parser.add_argument("--input", required=True)
    blast_parser.add_argument("--output", required=True)

    args = parser.parse_args()

    logging.info(f"Command received: {args.command}")

    try:
        if args.command == "fastq":
            logging.info(f"Running FASTQ filtering: {args.input} -> {args.output}")

            filter_fastq(
                input_fastq=args.input,
                output_fastq=args.output,
                gc_bounds=(args.min_gc, args.max_gc),
                length_bounds=(args.min_len, args.max_len),
                quality_threshold=args.quality,
            )

        elif args.command == "fasta":
            logging.info(f"Running FASTA conversion: {args.input}")

            convert_multiline_fasta_to_oneline(
                input_fasta=args.input,
                output_fasta=args.output,
            )

        elif args.command == "blast":
            logging.info(f"Running BLAST parsing: {args.input}")

            parse_blast_output(
                input_file=args.input,
                output_file=args.output,
            )

    except Exception as e:
        logging.error(f"Error occurred: {str(e)}")
        raise

if __name__ == "__main__":
    main()