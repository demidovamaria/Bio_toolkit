"""
Module for DNA/RNA sequence processing.
Includes functions for validating and manipulating nucleic acid sequences (reverse, complement, transcription, reverse-complement).
"""

def is_nucleic_acid(sequence: str) -> bool:
    """
    Validates whether a sequence is a DNA or RNA sequence.

    Checks if the sequence contains only valid nucleotide bases (A, T, U, C, G)
    and ensures it does not mix T (DNA) and U (RNA).

    Arguments:
        sequence: The input sequence to validate (case-insensitive).

    Returns:
        bool: True if the sequence is a valid DNA or RNA sequence, False otherwise.

    Examples:
        >>> is_nucleic_acid("ATCG")
        True
        >>> is_nucleic_acid("AUCG")
        True
        >>> is_nucleic_acid("ATUCG")
        False
        >>> is_nucleic_acid("ATXG")
        False
    """
    sequence = sequence.upper()
    bases = set("ATUCG")
    if "T" in sequence and "U" in sequence:
        return False
    return set(sequence).issubset(bases)

def transcribe(sequence: str) -> str:
    """
    Transcribes a DNA sequence to RNA by replacing T with U.

    The function is case-preserving: lowercase 't' becomes 'u', uppercase 'T' becomes 'U'.
    Non-T bases are returned unchanged.

    Arguments:
        sequence: The input DNA sequence (case-insensitive).

    Returns:
        str: The transcribed RNA sequence.

    Examples:
        >>> transcribe("ATCG")
        'AUCG'
        >>> transcribe("atcg")
        'aucg'
    """
    return "".join("U" if base.upper() == "T" else base for base in sequence)

def reverse(sequence: str) -> str:
    """
    Reverses the input sequence.

    Arguments:
        sequence: The input DNA/RNA sequence.

    Returns:
        str: The reversed sequence.

    Examples:
        >>> reverse("ATCG")
        'GCTA'
        >>> reverse("AUCG")
        'GCUA'
    """
    return sequence[::-1]

def complement(sequence: str) -> str:
    """
    Computes the complement of a DNA or RNA sequence.

    For DNA: A<->T, C<->G. For RNA: A<->U, C<->G. The function preserves the case
    of the input (e.g., lowercase 'a' complements to 't' or 'u'). The sequence type
    (DNA or RNA) is determined by the presence of 'U' (RNA) or absence (DNA).

    Arguments:
        sequence: The input DNA/RNA sequence (case-insensitive).

    Returns:
        str: The complementary sequence.

    Examples:
        >>> complement("ATCG")
        'TAGC'
        >>> complement("AUCG")
        'UAGC'
        >>> complement("atcg")
        'tagc'
    """
    seq = sequence.upper()
    is_rna = "U" in seq
    complement_map = (
        str.maketrans("AUCGaucg", "UAGCuagc")
        if is_rna
        else str.maketrans("ATCGatcg", "TAGCtagc")
    )
    return sequence.translate(complement_map)

def reverse_complement(sequence: str) -> str:
    """
    Computes the reverse complement of a DNA or RNA sequence.

    First computes the complement of the sequence (A<->T/U, C<->G), then reverses it.
    Preserves the case of the input sequence.

    Arguments:
        sequence: The input DNA/RNA sequence (case-insensitive).

    Returns:
        str: The reverse complement sequence.

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("AUCG")
        'CGAU'
        >>> reverse_complement("atcg")
        'cgat'
    """
    return reverse(complement(sequence))