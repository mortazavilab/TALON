# Toy transcript test case for TALON
Unlike other test GTFs in the testing suite, this test case consists of an articifial transcript annotation. Those transcripts include:

### TG1-001 on chromosome 1 (+)
* Exon 1: 1-100  
* Exon 2: 500-600  
* Exon 3: 900-1000
### TG2-001 on chromosome 2 (+)
Same structure as TG1-001, but it is on chromosome 2 instead of chromosome 1.
* Exon 1: 1-100
* Exon 2: 500-600
* Exon 3: 900-1000
### TG3-001 on chromosome 1 (-)
Shares an exon with TG1-001, but is on the opposite strand
* Exon 1: 1500-2000
* Exon 2: 900-1000

One test intended for this case is to look at how TALON and the post-TALON filtering behave on a transcript containing only exons 2 and 3 (with possible 5' end differences in either direction). Two sam transcripts are included for this purpose. The first matches exons 2 and 3, but starts further downstream than expected. The second starts further upstream than expected.

TG1-001 and TG2-001 are useful for testing the AND properties of TALON queries. When we query a position shared by TG1-001 and TG2-001 (except for the chromosome of course), we should only get one result back if we specified a particular chromosome.
