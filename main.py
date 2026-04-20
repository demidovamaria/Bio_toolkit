from typing import Dict, Tuple, Union
import os

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

    SeqIO.write(passed_records, output_fastq, "fastq")

filter_fastq(
    input_fastq="data/example_fastq.fastq",
    output_fastq="filtered_out.fastq",
    gc_bounds=(0.15, 1.0),
    length_bounds=(10, 2**32),
)

convert_multiline_fasta_to_oneline(
    input_fasta='data/example_multiline_fasta.fasta',
)

parse_blast_output(
    input_file='data/example_blast_results.txt',
    output_file='parsed_blast_result.txt'
)