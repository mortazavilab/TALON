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
  -h, --help           Show help message and exit  
  --f                  Comma-delimited dataset config file providing sam files for TALON to run on, as well as metadata that will be tracked in the dataset table. The required format is: dataset name, sample description, platform, sam file (full path).  
  -a, --annot          TALON database. Created using build_talon_annotation.py  
  -b BUILD, --build    Genome build to use. Note: must be in the TALON database.  
  --o=FILE              Outfile name  
  --encode              If this option is set, TALON will require novel transcripts to be corroborated by at least one other                  dataset in order to be included in the output abundance file.  
```
