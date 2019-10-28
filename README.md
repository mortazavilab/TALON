<img align="left" width="400" src="diagram.png">
TALON is a Python program for identifying known and novel genes/isoforms
in long read transcriptome data sets. TALON is technology-agnostic in that it
works from mapped SAM files, allowing data from different sequencing platforms
(i.e. PacBio and Oxford Nanopore) to be analyzed side by side.

To learn more, please see our preprint in BioRxiv: https://www.biorxiv.org/content/10.1101/672931v1

# Installation
Newer version of TALON (v4.0+) are designed to be run with Python 3.6+. 

To install TALON, simply download the files using Github's "Download ZIP" button, then unzip them in the directory where you would like to store the program. Alternately, you can download a specific version of the program from the Releases tab.

Go to the directory and run `pip install .`. This will install TALON. You can now run the commands from anywhere.

If the above installation doesn't work, it's likely due to difficulties with MacOS allowing installation of pybedtools. A workaround is to install pybedtools using conda. To do so run `conda install -c bioconda pybedtools`, then try to `pip install .` from the TALON directory again.

NOTE: Talon versions 4.2 and lower are not installable. Check the README of those releases to see how you can run the scripts from the install directory, or visit the wiki here: https://github.com/dewyman/TALON/wiki/Archived-TALON-documentation.

# How to run
For a small, self-contained example with all necessary files included, see https://github.com/dewyman/TALON/tree/master/example

## Initializing a TALON database
For documentation of TALON versions 4.2 and lower, see https://github.com/dewyman/TALON/wiki/Archived-TALON-documentation.

The first step in using TALON is to initialize a SQLite database from the GTF annotation of your choice (i.e. GENCODE). This step is done using talon_initialize_database, and only needs to be performed once. Keep track of the build and annotation names you choose, as these will be used downstream when running TALON and its utilities.

NOTE: The GTF file you use must contain genes, transcripts, and exons. If the file does not contain explicit gene and/or transcript entries, key tables of the database will be empty and you will experience problems in the downstream analysis. We have included a tool, talon_reformat_gtf, that can convert this type of GTF into the proper format.

```
talon_initialize_database --h

Usage: talon_initialize_database [options]

Options:
  -h, --help           Show help message and exit
  --f                  GTF annotation file
  --g                  The name of the reference genome build that the annotation describes. Use a short and memorable name since you will need to specify the genome build when you run TALON later.
  --a                  The name of the annotation (for metadata purposes)
  --l                  Minimum required transcript length (default = 300 bp)
  --idprefix           Prefix for naming novel discoveries in eventual TALON runs (default = 'TALON')
  --5p                 Maximum allowable distance (bp) at the 5' end during annotation (default = 500 bp)
  --3p                 Maximum allowable distance (bp) at the 3' end during annotation (default = 300 bp)
  --o                  Output prefix for the database
```

## Running TALON
Now that you've initilialized your database, you're ready to annotate long read datasets using TALON. The input database is modified in place to track and quantify transcripts in the provided dataset(s). You can add more datasets at any time by creating a config file for them and running this command. Please note that TALON versions 4.4+ can be run in multithreaded fashion for a much faster runtime.

```
talon --h

Usage: talon [options]

Options:
-h, --help            Show help message and exit
--f                   Comma-delimited dataset config file providing sam files for TALON to run on, as well as metadata that   will be tracked in the dataset table. The required format is: dataset name, sample description, platform, sam file (full path).
  --db FILE,            TALON database. Created using talon_initialize_database.
  --build STRING,       Genome build (i.e. hg38) to use. Must be in the
                        database.
  --threads THREADS, -t THREADS
                        Number of threads to run program with. Default = 2.
  --cov, -c             Minimum alignment coverage in order to use a SAM entry. Default = 0.9
  --identity, -i        Minimum alignment identity in order to use a SAM entry. Default = 0
  --o OUTPREFIX         Prefix for output files
```
TALON generates two output files in the course of a run. The QC log (file with suffix 'QC.log') contains the following columns:
1. dataset	
2. read_ID	
3. passed_QC (1/0)	
4. primary_mapped	(1/0)
5. read_length	
6. fraction_aligned	
7. identity
This is useful for tracking why a particular read was or was not included in the TALON analysis.

The second output file (suffix 'read_annot.tsv') appears at the very end of the run and contains a line for every read that was successfully annotated. The columns are as follows:
1. Name of individual read
2. Name of dataset the read belongs to
3. Name of genome build used in TALON run
4. Chromosome
5. Read start position (1-based). This refers to the 5' end start, so for reads on the - strand, this number will be larger than the read end (col 6).
6. Read end position (1-based). This refs to the 3' end stop, so for reads on the - strand, this will be smaller than the read start (col 5).
7. Strand (+ or -)
8. Number of exons in the transcript
9. Read length (soft-clipped bases not included)
10. Gene ID (issued by TALON, integer)
11. Transcript ID (issued by TALON, integer)
12. Annotation gene ID (ie Ensembl). Will start with 'ENCODE' if gene is a novel gene called by TALON pipeline.
13. Annotation transcript ID (ie Ensembl). Will start with 'ENCODE' if transcript is a novel transcript called by TALON pipeline.
14. Annotation gene name (human-readable gene symbol)
15. Annotion transcript name (human-readable transcript symbol)
16. Gene novelty: one of 'Known", "Antisense", or "Intergenic". Called by TALON pipeline.
17. Transcript novelty: one of "Known", "ISM", "NIC", "NNC", "Antisense", "Intergenic", or "Genomic". Called by TALON pipeline.
18. ISM subtype. If transcript novelty is not ISM, this field will be 'None'. If the transcript is an ISM, then this field can be 'Prefix', 'Suffix', 'Both', or 'None'.

Note: It is possible to obtain this file from a TALON database at any time by running this utility
```
talon_fetch_reads --h

Usage: talon_fetch_reads [-h] [--db FILE,] [--build STRING,]
                         [--datasets STRING,] [--o OUTPREFIX]

This utility queries a TALON database in order to get read-specific annotation
information.

optional arguments:
  -h, --help          show this help message and exit
  --db FILE,          TALON database
  --build STRING,     Genome build (i.e. hg38) to use. Must be in the
                      database.
  --datasets STRING,  Optional: Comma-delimited list of datasets to include.
                      Default behavior is to include all datasets in the
                      database.
  --o OUTPREFIX       Prefix for output files
```

## TALON utilities

### Filtering your transcriptome
If you have run TALON on biological replicates or other datasets you would like to leverage for quality control, you might want to obtain a filtered list of transcripts that are 1) known, or 2) reproducible in at least two of your datasets. To get such a list, run the following TALON utility:
```
talon_filter_transcripts --h

Usage: talon_filter_transcripts [options]

Options:
  -h, --help            show this help message and exit
  --db=FILE             TALON database
  -a ANNOT, --annot=ANNOT
                        The name of the annotation version to use.
                        Will determine which annotation transcripts
                        are considered known or novel relative to.
                        Note: must be in the TALON database.
  -p FILE, --pairings=FILE
                        Optional: A file indicating which datasets
                        should be considered together when filtering
                        novel transcripts (i.e. biological replicates).
                        Format: Each line of the file constitutes a group,
                        with member datasets separated by commas.
                        If no file is provided, then novel transcripts
                        appearing in any two datasets will be accepted.
  --o=FILE              Outfile name
```
The columns in the resulting output file are:
1. TALON gene ID (an integer). This is the same type of ID found in column 1 of TALON abundance files.
2. TALON transcript ID (an integer). This is the same type of ID found in column 2 of TALON abundance files.
3. Novelty category designation of transcript.


### Obtaining an abundance matrix from your TALON database
If you would like to extract an abundance matrix for your TALON-processed datasets, use the tool *talon_create_abundance_file_from_database*. You can provide the filtered transcript list you obtained from talon_filter_transcripts if you would like to restrict the abundance file to those transcripts.

```
talon_abundance --h

Usage: talon_abundance [options]

Options:
  -h, --help            show this help message and exit
  --db=FILE             TALON database
  -a ANNOT, --annot=ANNOT
                        Which annotation version to use. Will determine which
                        annotation transcripts are considered known or novel
                        relative to. Note: must be in the TALON database.
  -b BUILD, --build=BUILD
                        Genome build to use. Note: must be in the TALON
                        database.
  --whitelist=FILE      Whitelist file of transcripts to include in the
                        output. First column should be TALON gene ID,
                        second column should be TALON transcript ID.
                        Other columns are ignored.
  -d FILE, --datasets=FILE
                        Optional: A file indicating which datasets should be
                        included (one dataset name per line). Default is to
                        include all datasets.
  --o=FILE              Prefix for output file
```
The columns in the abundance file are as follows:
1. TALON gene ID
2. TALON transcript ID
3. Gene ID from your annotation of choice. If the gene is novel relative to that annotation, this will be 'NA'.
4. Transcript ID from your annotation of choice. If the transcript is novel relative to that annotation, this will be 'NA'.
5. Gene name from your annotation of choice (makes the file a bit more human-readable!). If the transcript is novel relative to that annotation, this will be the TALON-derived name.
6. Transcript name from your annotation of choice. If the transcript is novel relative to that annotation, this will be the TALON-derived name.
7. Number of exons in the transcript
8. Length of transcript model (basepairs)
9. Gene novelty (Known, Antisense, Intergenic)
10. Transcript status (Known, ISM, NIC, NNC, Antisense, Intergenic)
11. ISM subtype (Both, Prefix, Suffix, None)  
**---------------------------- Remaining columns -----------------------------**  
One column per dataset, with a count indicating how many times the current transcript was observed in that dataset.

### Obtaining a custom GTF transcriptome annotation from a TALON database

```
talon_create_GTF --h

Options:
  -h, --help            show this help message and exit
  --db=FILE             TALON database
  -b BUILD, --build=BUILD
                        Genome build to use. Note: must be in the TALON
                        database.
  -a ANNOT, --annot=ANNOT
                        Which annotation version to use. Will determine which
                        annotation transcripts are considered known or novel
                        relative to. Note: must be in the TALON database.
  --whitelist=FILE      Whitelist file of transcripts to include in the
                        output. First column should be TALON gene ID,
                        second column should be TALON transcript ID.
                        Other columns are ignored.
  --observed            If this option is set, the GTF file will only
                        include transcripts that were observed in at least one
                        dataset (redundant if dataset file provided).
  -d FILE, --datasets=FILE
                        Optional: A file indicating which datasets should be
                        included (one dataset name per line). Default is to
                        include                   all datasets.
  --o=FILE              Prefix for output GTF
```

# License
MIT, see LICENSE
