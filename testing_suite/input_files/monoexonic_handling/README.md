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
