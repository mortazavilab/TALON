# TALON example

Use the provided files to try TALON out on chromosome 21 PacBio SAM reads from two biological replicates of human cell line GM12878. 

First, initialize the TALON database from the provided GTF annotation of chromosome 21:
```
python ../initialize_talon_database.py \
        --f gencode.v29.chr21.gtf \
        --a gencode_v29 \
        --g hg38 \
        --o example_talon
```



Now, use the provided config file and the newly initialized database to annotate and quantify the reads. The database is modified in place.
```
python ../talon.py \
       --f config.csv \
       --db example_talon.db \
       --build hg38 \
       --o example
```
The file 'example_talon_QC.log' contains a record of each input transcript along with the computed coverage and identity of the alignment. During a longer TALON run, you can count the number of lines in the file in order to track how many reads have been processed so far.

Now, we can explore the results of the run. To summarize how many of each transxcript were found (prior to any filtering), run:
```
python ../post-TALON_tools/summarize_datasets.py \
       --db example_talon.db \
       --v \
       --o example
```

To create an abundance matrix without filtering (for use computing gene expression), we run the following:
```
python ../post-TALON_tools/create_abundance_file_from_database.py \
       --db example_talon.db \
       -a gencode_v29 \
       --build hg38 \
       --o example
```

We can also generate a filtered abundance matrix where we require transcripts to be either a) known, or b) detected in both replicates (for transcript-level expression):
```
python ../post-TALON_tools/create_abundance_file_from_database.py \
       --db example_talon.db \
       -a gencode_v29 \
       --filter \
       -p pairings.csv \
       --build hg38 \
       --o example
```


