Bio Utils
A lightweight Python package for bioinformatics tasks, designed for educational purposes. It provides tools for processing DNA/RNA sequences, filtering FASTQ sequences, converting multiline FASTA files to single-line format, and parsing BLAST output files. The package relies solely on the Python standard library, ensuring ease of use without external dependencies.
Features

DNA/RNA Processing: Validate, reverse, complement, transcribe, and compute reverse-complement of DNA/RNA sequences with case-preserving operations.
FASTQ Filtering: Filter FASTQ reads by GC content (percentage bounds), sequence length (bounds), and average Phred33 quality scores.
FASTA Processing: Convert multiline FASTA files to single-line format for easier downstream processing.
BLAST Output Parsing: Extract and sort unique top-match protein names from BLAST output text files.
Modular design: Core logic is organized into separate modules (rna_dna_utils, fastq, bio_files_processor, parse_blast_output) for reusability.
Type annotations and comprehensive docstrings for improved code readability and IDE support.
Flexible bounds: Support for single-value (upper limit) or tuple-based (range) bounds for filtering.

Usage
The package is organized into several modules, each handling specific bioinformatics tasks. Below are examples of how to use the core functionalities provided by the main.py script and its dependencies.
DNA/RNA Sequence Processing
The rna_dna_utils module provides functions to validate and manipulate nucleic acid sequences. These can be accessed via the run_dna_rna_tools function in main.py.
from bio_utils.main import run_dna_rna_tools

# Validate a DNA sequence
result = run_dna_rna_tools("ATCG", "is_nucleic_acid")
print(result)  # True

# Transcribe DNA to RNA
result = run_dna_rna_tools("ATCG", "transcribe")
print(result)  # AUCG

# Get reverse complement of a sequence
result = run_dna_rna_tools("ATCG", "reverse_complement")
print(result)  # CGAT

FASTQ Filtering
The fastq module supports filtering FASTQ files based on GC content, sequence length, and quality scores. Use the filter_fastq function in main.py to apply these filters.
from bio_utils.main import filter_fastq

# Filter FASTQ file with GC content between 15-100% and minimum length of 10
filter_fastq(
    input_fastq='data/example_fastq.fastq',
    output_fastq='filtered_out.fastq',
    gc_bounds=(15, 100),
    length_bounds=(10, 2**32),
    quality_threshold=0
)

Output is written to filtered/filtered_out.fastq.
FASTA File Conversion
The bio_files_processor module converts multiline FASTA files to single-line format using the convert_multiline_fasta_to_oneline function.
from bio_utils.bio_files_processor import convert_multiline_fasta_to_oneline

# Convert multiline FASTA to single-line
sequences = convert_multiline_fasta_to_oneline(
    input_fasta='data/example_multiline_fasta.fasta',
    output_fasta='processed/oneline_fasta.fasta'
)

Output is written to processed/oneline_fasta.fasta if output_fasta is not specified, or to the specified output_fasta path.
BLAST Output Parsing
The parse_blast_output module extracts and sorts unique top-match protein names from BLAST output text files.
from bio_utils.parse_blast_output import parse_blast_output

# Parse BLAST output to extract sorted protein names
parse_blast_output(
    input_file='data/example_blast_results.txt',
    output_file='parsed_blast_result.txt'
)

Output is written to parsed_blast_result.txt.
Module Overview

rna_dna_utils.py: Functions for validating and manipulating DNA/RNA sequences (is_nucleic_acid, transcribe, reverse, complement, reverse_complement).
fastq.py: Tools for reading, filtering, and writing FASTQ files based on GC content, length, and quality scores.
bio_files_processor.py: Converts multiline FASTA files to single-line format.
parse_blast_output.py: Parses BLAST output files to extract and sort unique top-match protein names.
main.py: Orchestrates the above functionalities with high-level functions for sequence processing and file handling.