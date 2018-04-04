TALON is a Python program for identifying known and novel genes/isoforms
in long read transcriptome data sets. TALON is technology-agnostic in that it
works from mapped SAM files, allowing data from different sequencing platforms
(i.e. PacBio and Oxford Nanopore) to be analyzed side by side. 

## Installation

Requires:
* intervaltree (v2.1.0)

## How to run

```
Usage: talon.py [options]

Options:
  -h, --help           show this help message and exit
  --f=FILE             Comma-delimited list of input SAM files
  -g FILE, --gtf=FILE  GTF annotation containing genes, transcripts, and exons

```
