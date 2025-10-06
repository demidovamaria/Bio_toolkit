"""
Module for FASTQ sequence processing.
Includes functions for filtering FASTQ sequences by GC content, length, and quality.
"""

from typing import Dict, Tuple, Union

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

def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[float, Tuple[float, float]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: float = 0
) -> Dict[str, Tuple[str, str]]:
    """
    Filters FASTQ sequences based on GC content, length, and quality thresholds.

    Arguments:
        seqs: Dictionary of sequences {name: (sequence, quality)}.
        gc_bounds: GC content range in percentage (single number or tuple).
        length_bounds: Sequence length range (single number or tuple).
        quality_threshold: Minimum average Phred33 quality score.

    Returns:
        Dictionary with filtered sequences.
    """
    filtered_seqs = {}

    for name, (seq, qual) in seqs.items():
        # Check sequence length
        seq_length = len(seq)
        if not is_within_bounds(seq_length, length_bounds):
            continue

        # Check GC content
        gc_content = calculate_gc_content(seq)
        if not is_within_bounds(gc_content, gc_bounds):
            continue

        # Check mean quality
        mean_quality = calculate_mean_quality(qual)
        if mean_quality < quality_threshold:
            continue

        # If all checks pass, add to result
        filtered_seqs[name] = (seq, qual)
    return filtered_seqs