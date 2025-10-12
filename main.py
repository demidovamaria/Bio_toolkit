from typing import Dict, Tuple, Union
import os

from bio_utils.bio_files_processor import convert_multiline_fasta_to_oneline
from bio_utils.fastq import is_within_bounds, calculate_gc_content, calculate_mean_quality, read_fastq_to_dict, \
    write_fastq
from bio_utils.parse_blast_output import parse_blast_output

from bio_utils.rna_dna_utils import is_nucleic_acid, transcribe, reverse, complement, reverse_complement


def run_dna_rna_tools(*args):
    if not args or len(args) < 2:
        raise Exception("Input error")
    sequences = args[:-1]
    procedure = args[-1]

    procedure_map = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    if procedure not in procedure_map:
        raise Exception("Input error")

    results = [procedure_map[procedure](sequence) for sequence in sequences]
    return results[0] if len(results) == 1 else results


def filter_fastq(
        input_fastq: str,
        output_fastq: str,
        gc_bounds: Union[float, Tuple[float, float]] = (0, 100),
        length_bounds: Union[int, Tuple[int, int]] = (0, 2 ** 32),
        quality_threshold: float = 0
) -> None:
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

    # Create filtered directory if it doesn't exist
    output_dir = "filtered"
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, output_fastq)

    # Read FASTQ file into dictionary
    seqs = read_fastq_to_dict(input_fastq)

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

    write_fastq(
        output_path, filtered_seqs
    )


filter_fastq(
    'data/example_fastq.fastq',
    'filtered_out.fastq',
    (15, 100),
    (10, 2 ** 32)
)

convert_multiline_fasta_to_oneline(
    input_fasta='data/example_multiline_fasta.fasta',
)

parse_blast_output(
    input_file='data/example_blast_results.txt',
    output_file='parsed_blast_result.txt'
)