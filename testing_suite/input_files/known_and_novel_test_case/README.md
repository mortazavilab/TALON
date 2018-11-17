# Known/novel test case for TALON
The GTF file contains annotations for 3 genes from GENCODE v24:
* LINC01128
* NOC2L
* RPS6KA1

The first sam file, known_novel_test_case.sam, contains 4 transcripts. Transcripts
1 and 4 are known (RPS6KA1-008 and NOC2L-001, respectively). Transcripts 2 and 3
are the same novel transcript of LINC01128. The correct TALON behavior for this case
is to create a novel transcript when it sees transcript 2, and then assign 
transcript 3 to that same novel transcript.

The second sam file, novel_test_case.sam contains two novel transcripts. The first 
is identical to transcripts 2 and 3 from before, and the second is an EBV transcript
that should be labeled as novel because EBV is missing from the annotation. 
This is designed to test the post-TALON filtering step. The filter is supposed to
whitelist the novel LINC01128 transcript because it appears in two datasets, but
should reject the EBV transcript because it is only seen once. All known transcripts
should be whitelisted. 
