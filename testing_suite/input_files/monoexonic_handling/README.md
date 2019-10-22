# Monoexonic handling

These test files help address the issue of how TALON annotates single-exon transcripts belonging to single-exon genes. In the past, an unitentional edge case led to the following scenario:

Reference transcript:       ===================================
Long read:       diff > min 5' diff ======================= diff > min 3' diff

--> Outcome: long read gets annotated as 'genomic' even though it is a good match to the reference transcript.

What we want instead is for such a transcript to be annotated as known. However, there could be unintended effects of entirely removing diff constraints. For instance, consider a possible miRNA gene:

Reference miRNA:                          ===
Long read:                      =======================

Ideally, we would not want this read to be assigned as an miRNA transcript

## Cases in monoexon.gtf:
Monoexonic gene from chr1:1000-2000(+)
Monoexonic gene from chr1:1500-1550(+)

## Cases in monoexon_reads.sam
Monoexonic transcript from chr1:1010-1909(+). This transcript overlaps both genes and has a 5' end diff of 10 bp and a 3' end diff of -91. If we run TALON with allowable 5' and 3' end differences of 99 bp at each end, then this transcript will be annotated as known.
Monoexonic transcript from chr1:1100-1599(+). This means that it overlaps both genes, but has a 100 bp 5' end difference and a -401 bp difference at the 3' end from the larger gene. Using the same threshold as for the first transcript, the edge case behovior would be to call this transcript as genomic. 
