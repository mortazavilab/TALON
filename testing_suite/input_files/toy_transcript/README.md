# Toy transcript test case for TALON
Unlike other test GTFs in the testing suite, this test case consists of an articifial transcript annotation. The transcript contains three exons with the following position ranges:  
* Exon 1: 1-100  
* Exon 2: 500-600  
* Exon 3: 900-1000

One test intended for this case is to look at how TALON and the post-TALON filtering behave on a transcript containing only exons 2 and 3 (with possible 5' end differences in either direction). Two sam transcripts are included for this purpose. The first matches exons 2 and 3, but starts further downstream than expected. The second starts further upstream than expected.
