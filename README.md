TALON is a Python program for identifying known and novel genes/isoforms
in long read transcriptome data sets. TALON is technology-agnostic in that it
works from mapped SAM files, allowing data from different sequencing platforms
(i.e. PacBio and Oxford Nanopore) to be analyzed side by side. 

## Installation
TALON is designed to be run with Python version 2.7.

Requires:
* intervaltree (v2.1.0)

To install TALON, simply download the files using Github's "Download ZIP" button, then unzip them in the directory where you would like to install the program. Alternately, you can download a specific version of the program from the Releases tab. The TALON scripts can now be run directly from the command line- just include the path.

## How to run

### Initializing a TALON database
The first step in using TALON is to initialize a SQLite database from the GTF annotation of your choice (i.e. GENCODE). This step is done using initialize_talon_database.py, and only needs to be performed once.  

```
Usage: initialize_talon_database.py [options]

Options:
  -h, --help           Show help message and exit
  --f                  GTF annotation file
  --a                  The name of the annotation (for metadata purposes)
  --g                  The name of the reference genome build that the annotation describes. Use a short and memorable name since you will need to specify the genome build when you run TALON later.
  --o                  Output prefix for the database
```

### Running TALON
Now that you've initilialized your database, you're ready to annotate long read datasets using TALON. 

```
Usage: talon.py [options]

Options:
  -h, --help           Show help message and exit  
  --f                  Comma-delimited dataset config file providing sam files for TALON to run on, as well as metadata that   will be tracked in the dataset table. The required format is: dataset name, sample description, platform, sam file (full path).  
  -a, --annot          TALON database. Created using build_talon_annotation.py  
  -b, --build          Genome build to use. Note: must be in the TALON database.  
  --o                  Outfile name  
  --encode             If this option is set, TALON will require novel transcripts to be corroborated by at least one other        dataset in order to be included in the output abundance file.  
```
