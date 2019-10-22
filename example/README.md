# TALON example

## User note
This example was initially constructed for TALON v4.2. In the packaged versions since, the commands have changed since TALON does not need to be run from a specific path. If you are using a TALON version below 4.3, please see the archived documentation here for the correct commands: https://github.com/dewyman/TALON/wiki/Archived-TALON-Example-Instructions-(v4.2). 

## Tutorial
Use the provided files to try TALON out on chromosome 21 PacBio SAM reads from two biological replicates of human cell line GM12878. 

First, initialize the TALON database from the provided GTF annotation of chromosome 21:

```
talon_initialize_database \
        --f gencode.v29.chr21.gtf \
        --a gencode_v29 \
        --g hg38 \
        --o example_talon
```

Now, use the provided config file and the newly initialized database to annotate and quantify the reads. The database is modified in place.
```
talon \
       --f config.csv \
       --db example_talon.db \
       --build hg38 \
       --o example
```
The file 'example_talon_QC.log' contains a record of each input transcript along with the computed coverage and identity of the alignment. During a longer TALON run, you can count the number of lines in the file in order to track how many reads have been processed so far.

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
       -a gencode_v29 \
       --build hg38 \
       --o example
```

We can also generate a filtered abundance matrix where we require transcripts to be either a) known, or b) detected in both replicates (for transcript-level expression). First, we run the TALON filtering script as follows:
```
talon_filter_transcripts \
       --db example_talon.db \
       -a gencode_v29 \
       -p pairings.csv \
       --o filtered_transcripts.csv
```
Then, we can run the abundance script, using the filtered transcript list as input.
```
talon_abundance \
       --db example_talon.db \
       -a gencode_v29 \
       --build hg38 \
       --whitelist filtered_transcripts.csv \
       --o example
```


