TALON is a Python program for identifying known and novel genes/isoforms
in long read transcriptome data sets. TALON is technology-agnostic in that it
works from mapped SAM files, allowing data from different sequencing platforms
(i.e. PacBio and Oxford Nanopore) to be analyzed side by side. 

# Installation
TALON is designed to be run with Python version 2.7.

Requires:
* intervaltree (v2.1.0)

To install TALON, simply download the files using Github's "Download ZIP" button, then unzip them in the directory where you would like to install the program. Alternately, you can download a specific version of the program from the Releases tab. The TALON scripts can now be run directly from the command line- just include the path.

# How to run

## Initializing a TALON database
The first step in using TALON is to initialize a SQLite database from the GTF annotation of your choice (i.e. GENCODE). This step is done using initialize_talon_database.py, and only needs to be performed once.  

```
python initialize_talon_database.py --h

Usage: initialize_talon_database.py [options]

Options:
  -h, --help           Show help message and exit
  --f                  GTF annotation file
  --a                  The name of the annotation (for metadata purposes)
  --g                  The name of the reference genome build that the annotation describes. Use a short and memorable name since you will need to specify the genome build when you run TALON later.
  --o                  Output prefix for the database
```

## Running TALON
Now that you've initilialized your database, you're ready to annotate long read datasets using TALON. 

```
python talon.py --h

Usage: talon.py [options]

Options:
  -h, --help           Show help message and exit  
  --f                  Comma-delimited dataset config file providing sam files for TALON to run on, as well as metadata that   will be tracked in the dataset table. The required format is: dataset name, sample description, platform, sam file (full path).  
  -a, --annot          TALON database. Created using build_talon_annotation.py  
  -b, --build          Genome build to use. Note: must be in the TALON database.  
  --o                  Outfile name
```

## Post-TALON utilities

### Obtaining an abundance matrix from your TALON database
If you would like to extract an abundance matrix for your TALON-processed datasets, use the script *create_abundance_file_from_database.py* from the post-TALON_tools directory.

```
python post-TALON_tools/create_abundance_file_from_database.py --h

Usage: create_abundance_file_from_database.py [options]

Options:
  -h, --help            show this help message and exit
  --db=FILE             TALON database
  -a ANNOT, --annot=ANNOT
                        Which annotation version to use. Will determine which
                        annotation transcripts are considered known or novel
                        relative to. Note: must be in the TALON database.
  --filter              If this option is set, the transcripts in the
                        database will be filtered prior to GTF creation
                        (for more information, see
                        filter_talon_transcripts.py)
  -p FILE, --pairings=FILE
                        Optional (only relevant if filter = true): A file
                        indicating which datasets should be
                        considered together when filtering
                        novel transcripts (i.e. biological replicates).
                        Format: Each line of the file constitutes a group,
                        with member datasets separated by
                        commas. If no file is provided, then
                        novel transcripts appearing in any
                        two datasets will be accepted.
  --o=FILE              Prefix for output file
```
The columns in the abundance file are as follows:
1. TALON gene ID
2. TALON transcript ID	
3. Gene ID from your annotation of choice. If the gene is novel relative to that annotation, this will be 'NA'.
4. Transcript ID from your annotation of choice. If the transcript is novel relative to that annotation, this will be 'NA'.
5. Gene name from your annotation of choice (makes the file a bit more human-readable!). If the transcript is novel relative to that annotation, this will be 'NA'.
6. Transcript name from your annotation of choice. If the transcript is novel relative to that annotation, this will be 'NA'.	
7. Number of exons in the transcript
8. Gene status (KNOWN or NOVEL)	
9. Transcript status (KNOWN or NOVEL)  
10. One column per dataset, with a count indicating how many times the current transcript was observed in that dataset.

# License
MIT, see LICENSE
