# TALON example

## User note
This example was initially constructed for TALON v4.2. In the packaged versions since, the commands have changed since TALON does not need to be run from a specific path. If you are using a TALON version below 4.3, please see the archived documentation [here](https://github.com/mortazavilab/TALON/wiki/Archived-TALON-Example-Instructions-(v4.2)) for the correct commands. In addition, the full tutorial as of v4.4.2 is available [here](https://github.com/mortazavilab/TALON/wiki/Archived-TALON-Example-Instructions-(v4.4)).

## Tutorial
Use the provided files to try TALON out on spike-in RNA variant (SIRV) reads sequenced on the PacBio Sequel II platform. 

### Database initialization
First, initialize the TALON database from the provided SIRV GTF annotation:

```
talon_initialize_database \
        --f SIRV_annotation.gtf \
        --a SIRV_annot \
        --g SIRV \
        --o example_talon
```

### Internal priming check
Before annotating the reads, we run talon_label_reads on each file in order to compute how likely each read is to be an internal priming product. Since this is a labeling step, no reads are removed- the script simply annotates each SAM read with the fraction of As present in the 20 bases immediately after the end of the alignment. This is why we need the reference genome fasta that the reads were aligned to. 
```
mkdir -p labeled
talon_label_reads --f aligned_reads/SIRV_rep1.sam \
    --g SIRV.fa  \
    --t 1 \
    --ar 20 \
    --deleteTmp \
    --o labeled/SIRV_rep1

talon_label_reads --f aligned_reads/SIRV_rep2.sam \
    --g SIRV.fa  \
    --t 1 \
    --ar 20 \
    --deleteTmp \
    --o labeled/SIRV_rep2
```
Note: for bigger files, we use the --t option to specify the number of threads for parallelized processing. Here, it should not make a big difference.

### Run TALON annotator
Now, use the provided config file and the newly initialized database to annotate and quantify the reads. The database is modified in place.
```
talon \
       --f config.csv \
       --db example_talon.db \
       --build SIRV \
       --o example
```
The file 'example_talon_QC.log' contains a record of each input transcript along with the computed coverage and identity of the alignment. During a longer TALON run, you can count the number of lines in the file in order to track how many reads have been processed so far. You can also specify the number of threads to run in parallel with the --t option.

### Abundance and filtering
Now, we can explore the results of the run. To summarize how many of each transxcript were found (prior to any filtering), run:
```
talon_summarize \
       --db example_talon.db \
       --v \
       --o example
```
You may notice at this point that there are a lot of novel transcripts. This is because no filtering has been performed. We highly recommend processing biological replicates together and then filtering transcript models for reproducibility (as described later in this tutorial).

To create an abundance matrix without filtering (for use computing gene expression), we run the following:
```
talon_abundance \
       --db example_talon.db \
       -a SIRV_annot \
       --build SIRV \
       --o example
```
The output file is named 'example_talon_abundance.tsv'.

We can also generate a whitelist file containing only transcripts pass the TALON filters. With the default settings, a transcript must be either a) known, or b) detected 5 times in both replicates. For condition b, all of the supporting reads must have 50% or fewer As in the 20 bp interval after alignment. This is the information that was recorded by our talon-label_reads step.
```
talon_filter_transcripts \
       --db example_talon.db \
       --datasets SIRV_Rep1,SIRV_Rep2 \
       -a SIRV_annot \
       --maxFracA 0.5 \
       --minCount 5 \
       --minDatasets 2 \
       --o filtered_transcripts.csv
```
Then, we can run the abundance script, using the filtered transcript list as input. This will give us the filtered abundance table appropriate for isoform quantification.
```
talon_abundance \
       --db example_talon.db \
       --whitelist filtered_transcripts.csv \
       -a SIRV_annot \
       --build SIRV \
       --o example
```
The output file is named 'example_talon_abundance_filtered.tsv'.

We can also use the whitelist to generate a custom GTF annotation of the filtered transcripts in our SIRV datasets like this:
```
talon_create_GTF \
       --db example_talon.db \
       --whitelist filtered_transcripts.csv \
       -a SIRV_annot \
       --build SIRV \
       --o example
```
The output file is named 'example_talon.gtf'.
