"""
Module for FASTQ sequence processing.
Includes functions for filtering FASTQ sequences by GC content, length, and quality.
"""

from typing import Dict, Tuple, Union, Iterator, TextIO


def calculate_gc_content(seq: str) -> float:
    """
    Calculates the GC content of a sequence in percentage.

    Arguments:
        seq: Input DNA/RNA sequence.

    Returns:
        GC content as a percentage.
    """
    if not seq:
        return 0.0
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100


def calculate_mean_quality(quality_str: str) -> float:
    """
    Calculates the mean Phred33 quality score of a FASTQ read.

    Arguments:
        quality_str: Quality string in Phred33 format.

    Returns:
        Mean quality score.
    """
    if not quality_str:
        return 0.0
    total_quality = sum(ord(char) - 33 for char in quality_str)
    return total_quality / len(quality_str)


def is_within_bounds(value: float, bounds: Union[float, Tuple[float, float]]) -> bool:
    """
    Checks if a value is within the specified bounds.

    Arguments:
        value: Value to check.
        bounds: Single number (upper bound) or tuple (lower, upper bounds).

    Returns:
        True if value is within bounds, False otherwise.
    """
    if isinstance(bounds, (int, float)):
        return value <= bounds
    return bounds[0] <= value <= bounds[1]


def read_fastq_to_dict(input_fastq: str) -> Dict[str, Tuple[str, str]]:
    """
    Read a FASTQ file and return a dictionary of sequences.

    Arguments:
        input_fastq: Path to the FASTQ file.

    Returns:
        Dictionary with sequence names as keys and (sequence, quality) tuples as values.
    """
    seqs = {}
    file = open(input_fastq, 'r')
    while True:
        name = file.readline().strip()
        if not name:
            break  # End of file
        if not name.startswith('@'):
            raise ValueError("Invalid FASTQ format: Expected '@' at start of name")

        seq = file.readline().strip()
        plus = file.readline().strip()
        qual = file.readline().strip()

        if not plus.startswith('+') or len(seq) != len(qual):
            raise ValueError("Invalid FASTQ format: Readline error")

        seqs[name] = (seq, qual)
    return seqs

def write_fastq(output_path: str, seqs: Dict[str, Tuple[str, str]]) -> None:
    """
    Write sequences from a dictionary to a FASTQ file.

    Arguments:
        output_path: Path to the output FASTQ file.
        seqs: Dictionary with sequence names as keys and (sequence, quality) tuples as values.
    """
    with open(output_path, 'w') as file:
        for name, (seq, qual) in seqs.items():
            file.write(f"{name}\n{seq}\n+\n{qual}\n")
