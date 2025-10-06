def is_nucleic_acid(sequence):
    sequence = sequence.upper()
    bases = set("ATUCG")
    if "T" in sequence and "U" in sequence:
        return False
    return all(base in bases for base in sequence)


def transcribe(sequence):
    if not is_nucleic_acid(sequence):
        raise Exception(
            "Error, only nucleic acids is available for processing"
        )
    return "".join("U" if base.upper() == "T" else base for base in sequence)


def reverse(sequence):
    if not is_nucleic_acid(sequence):
        raise Exception(
            "Error, only nucleic acids is available for processing"
        )
    return sequence[::-1]


def complement(sequence):
    if not is_nucleic_acid(sequence):
        raise Exception(
            "Error, only nucleic acids is available for processing"
        )
    seq = sequence.upper()
    is_rna = "U" in seq
    complement_map = (
        str.maketrans("AUCGaucg", "UAGCuagc")
        if is_rna
        else str.maketrans("ATCGatcg", "TAGCtagc")
    )
    return sequence.translate(complement_map)


def reverse_complement(sequence):
    if not is_nucleic_acid(sequence):
        raise Exception(
            "Error, only nucleic acids is available for processing"
        )
    return reverse(complement(sequence))


