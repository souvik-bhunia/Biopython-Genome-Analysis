# Biopython-Based Genome Sequence Analysis

A Python-based workflow for parsing FASTA sequences, calculating nucleotide composition and GC content, validating sequence symbols, and exporting structured sequence statistics to CSV.

## Project Overview

This project implements a sequence-level analysis workflow using Python and Biopython.

The workflow processes multiple FASTA files and extracts quantitative statistics for each sequence record. The calculated results are organised into a structured Pandas DataFrame and exported as a CSV file for downstream analysis.

## Workflow

```text
FASTA files
     │
     ▼
Biopython SeqIO
     │
     ▼
Sequence-level analysis
     │
     ├── Sequence length
     ├── A/T/G/C/N counts
     ├── GC content
     └── DNA symbol validation
     │
     ▼
GC calculation verification
     │
     ▼
Pandas DataFrame
     │
     ▼
CSV output
