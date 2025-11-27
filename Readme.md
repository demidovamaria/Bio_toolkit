# Bio Utils

A lightweight Python package for bioinformatics tasks, designed for educational purposes. It provides tools for processing DNA/RNA sequences (reversing, complementing, transcribing) and filtering FASTQ sequences based on GC content, length, and quality scores. The package uses only the Python standard library, making it easy to use without external dependencies.

## Features
- **DNA/RNA Processing**: Reverse, complement, and transcribe sequences.
- **FASTQ Filtering**: Filter reads by GC percentage (bounds), length (bounds), and average Phred33 quality.
- Modular design: Core logic separated into modules for reusability.
- Type annotations and docstrings for better code readability and IDE support.
- Supports flexible bounds: single values (upper limit) or tuples (range).

## Installation
1. Clone the repository from GitHub:
   ```bash
   git clone demidovamaria/Bio_toolkit.git
   cd Bio_toolkit