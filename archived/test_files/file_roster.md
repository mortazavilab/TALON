# File Roster

## gencode.v24.annotation.chr1.gtf
### Description:
All GENCODE v24 annotations on chromosome 1.

```
awk '{if($1 == "chr1") print $0}' /bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf > gencode.v24.annotation.chr1.gtf
```

### Uses:
Testing TALON on chromosome 1 datasets.

## sorted.GM12878_chr1_canonicalOnly.sam
### Description:
Sorted chromosome 1 transcripts from PacBio library PB14 (GM12878 cell line) PRIOR TO transcriptClean processing. Noncanonical transcripts have been removed. 

```
cp /pub/dwyman/clean_splice_jns/example_files/header.txt GM12878_chr1_canonicalOnly.sam
awk '{if($16 !~ "0") print $0}' /pub/dwyman/clean_splice_jns/example_files/GM12878_chr1.sam >> GM12878_chr1_canonicalOnly.sam
samtools sort GM12878_chr1_canonicalOnly.sam -o sorted.GM12878_chr1_canonicalOnly.sam -O sam
```
### Uses:
4/11/18: Trying to get Tofu to work. I think a bug in TranscriptClean may be affecting the corrected noncanonical transcripts

## GM12878_chr1_canonicalOnly.fq
### Description:
Fastq file derived from sorted.GM12878_chr1_canonicalOnly.sam.
```
samtools fastq -n sorted.GM12878_chr1_canonicalOnly.sam > GM12878_chr1_canonicalOnly.fq
```
### Uses:
4/11/18: Trying to get Tofu to work. I think a bug in TranscriptClean may be affecting the corrected noncanonical transcripts

## GM12878_chr1_clean.sam
### Description:
Chromosome 1 transcripts from PacBio library PB14 (GM12878 cell line) 
processed using TranscriptClean version 1.0 + log bug fix that has yet to be
put into a release.

Created using command:
```
    python TranscriptClean.py \
           --sam example_files/GM12878_chr1.sam \
           --genome test_TranscriptClean/reference_files/chr1.fa \
           --spliceJns example_files/GM12878_SJs_chr1.tab \
           --variants example_files/GM12878_chr1.vcf.gz \
           --outprefix /bio/dwyman/pacbio_f2016/TALON/test_files/GM12878_chr1
```

### Uses:
Good SAM input file that is large enough to contain many different transcript
cases, but small enough to run quickly. I used the uncorrected version of this 
dataset frequently to test and develop TranscriptClean. 

## sorted.GM12878_chr1_clean.sam
Description:
Sorted version of GM12878_chr1_clean.sam

```
samtools sort GM12878_chr1_clean.sam -o sorted.GM12878_chr1_clean.sam -O sam
```

## extraID.sorted.GM12878_chr1_clean.sam
###Description:
Sorted SAM file from TranscriptClean, processed to add "sample|" to the original ID.

Command:
```
awk '{if($1 !~ "@") print "sample|"$0; else print $0}' sorted.GM12878_chr1_clean.sam > extraID.sorted.GM12878_chr1_clean.sam
```

### Uses:
Tofu, since Tofu now requires a two-part read ID (sample|cluster_id)

## extraID.GM12878_chr1_clean.fa
###Description:
Fasta file from TranscriptClean, processed to add "sample|" to the original ID.

Command:
```
awk -F"c" -v OFS="" '{if($1 ~ ">") print $1,"sample|c"$2$3$4; else print $0}' GM12878_chr1_clean.fa > extraID.GM12878_chr1_clean.fa
```

### Uses:
Tofu, since Tofu now requires a two-part read ID (sample|cluster_id)

## test.gtf
### Description:
Small set of chromosome 1 and chromosome 2 gene annotations (with transcripts
and exons) derived from original file
/bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf

### Uses:
Good for GTF parsing development because the small file size allows lines to be
printed without lots of extra exit statements.
