![Dayhoff's Unzipper Logo](logo.jpg)
<h1 align="center">Dayhoff's Unzipper</h1>

A lightweight first-pass genome annotation tool written in Python.

It processes nucleotide sequences from FASTA files and performs basic structural and regulatory feature detection.

**Note:** This is an educational tool useful for quick inspection of sequences.

## Features

- Reads multi-FASTA files using Biopython
- Detects open reading frames (ORFs) in all three forward frames (ATG → stop codon)
- Searches for common promoter, regulatory, and translation initiation motifs using an efficient Rabin-Karp implementation
- Recognizes the following motifs:
  - TATA box (TATAAA)
  - CAAT box (CCAAT)
  - GC box (GGGCGG)
  - Pribnow box / -10 element (TATAAT)
  - AP-1 site (TGACTCA)
  - CRE site (TGACGTCA)
  - Octamer motif (ATGCAAAT)
  - Kozak sequence (GCCACCATGG)
  - Poly(A) signal (AATAAA)
  - Shine-Dalgarno sequence (AGGAGG)
- Outputs annotations in GFF3 format
- Includes `##sequence-region` pragma for each input sequence
- Processes files silently (no console output)

## First-Pass Annotation

Dayhoff's Unzipper is designed as a **quick first-pass annotation pipeline**:

1. **Input**: one or more nucleotide sequences in FASTA format
2. **Processing**:
   - Scans all three forward reading frames for potential protein-coding regions (ORFs)
   - Searches the entire sequence for exact matches to well-known regulatory and translation motifs
3. **Output**:
   - GFF3 file containing:
     - `CDS` entries for each detected ORF
     - `motif` entries for each detected regulatory sequence
   - All coordinates are 1-based, consistent with GFF3 specification

This allows visualization in genome browsers (IGV, JBrowse, etc.) or further processing with other bioinformatics tools.

## Requirements

- Python 3.6+
- Biopython (`pip install biopython`)

## Usage

Run the script and provide the path to your FASTA file when prompted:

```bash
python main.py
