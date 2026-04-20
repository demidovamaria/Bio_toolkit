# Bio Utils

A lightweight Python package for bioinformatics tasks, designed for educational purposes.
The project demonstrates clean Python code, object-oriented programming, and practical usage
of third-party bioinformatics libraries.

## Features

### Biological Sequences (OOP)
- Object-oriented representation of biological sequences
- Abstract base class `BiologicalSequence` defining a common interface
- DNA and RNA sequences implemented via inheritance and polymorphism
- Protein (amino acid) sequences with domain-specific methods

Supported operations:
- Length calculation (`len(sequence)`)
- Indexing and slicing
- Pretty string representation
- Alphabet validation
- DNA/RNA complement, reverse, and reverse-complement
- DNA transcription into RNA

### FASTQ Filtering (Biopython)
- FASTQ filtering using Biopython (`SeqIO`, `SeqRecord`, `SeqUtils`)
- Filtering by:
  - Sequence length
  - GC content (fraction)
  - Mean Phred quality score

### FASTA Processing
- Convert multiline FASTA files into single-line FASTA format

### BLAST Output Parsing
- Parse BLAST text output
- Extract unique top-hit protein names
- Sort results alphabetically

---

## Requirements

- Python 3.9+
- Biopython

Install dependencies with:

```bash
pip install -r requirements.txt