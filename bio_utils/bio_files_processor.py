import os
from typing import Optional, Dict


def convert_multiline_fasta_to_oneline(
        input_fasta: str,
        output_fasta: Optional[str] = None,
) -> Dict[str, str]:
    """
    Reads a multiline FASTA file, concatenates each sequence into a single line,
    and optionally writes to an output FASTA file.

    Arguments:
        input_fasta: Path to the input FASTA file.
        output_fasta: Optional path to the output FASTA file.

    Returns:
        Dictionary with sequence names as keys and concatenated sequences as values.
    """
    sequences = {}
    current_name = None
    current_seq = []

    # Read the input FASTA file
    with open(input_fasta, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            if line.startswith('>'):
                # Save the previous sequence if it exists
                if current_name and current_seq:
                    sequences[current_name] = ''.join(current_seq)
                    current_seq = []
                # Start a new sequence
                current_name = line[1:]  # Remove '>'
            else:
                # Append sequence line
                current_seq.append(line)

        # Save the last sequence
        if current_name and current_seq:
            sequences[current_name] = ''.join(current_seq)

    # Write to output file if specified
    if output_fasta:
        with open(output_fasta, 'w') as file:
            for name, seq in sequences.items():
                file.write(f">{name}\n{seq}\n")
    else:
        default_out_name = 'oneline_fasta.fasta'
        output_dir = "processed"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, default_out_name)
        with open(output_path, 'w') as file:
            for name, seq in sequences.items():
                file.write(f">{name}\n{seq}\n")

    return sequences