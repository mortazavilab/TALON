# File Roster

## gencode.v24.annotation.chr1.gtf
### Description:
All GENCODE v24 annotations on chromsome 1.

```
awk '{if($1 == "chr1") print $0}' /bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf > gencode.v24.annotation.chr1.gtf
```

### Uses:
Testing TALON on chromosome 1 datasets.

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

## test.gtf
### Description:
Small set of chromosome 1 and chromosome 2 gene annotations (with transcripts
and exons) derived from original file
/bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf

### Uses:
Good for GTF parsing development because the small file size allows lines to be
printed without lots of extra exit statements.
