# Chromosome 11 tests
These test cases are designed to among other things make sure that the queries used to access information in the TALON database are working as expected.

The database is initialized from the mouse chromosome 11 annotation plus the Tcf3 gene (which is on chromosome 10). The sam files contain reads mapped to the Tcf3, Drg1, and Canx loci. The idea is to annotate the following types of transcripts across multiple datasets, and then try to compute detetion metrics correctly using queries.

# GTF creation
```
awk '{if ($1 == "chr11") print $0}' /bio/dwyman/pacbio_f2016/data/GENCODEvM7/gencode.vM7.annotation.gtf > chr11_and_Tcf3.gtf
grep "Tcf3" /bio/dwyman/pacbio_f2016/data/GENCODEvM7/gencode.vM7.annotation.gtf >> chr11_and_Tcf3.gtf
```

# Read assignments:
Files used:
```
/share/samdata/dwyman/PB65/TC_v1.0.7/IS_primer_bc1017--IS_primer_bc1017/sorted_canonical_noERCC.sam
/share/samdata/dwyman/PB65/TC_v1.0.7/IS_primer_bc1018--IS_primer_bc1018/sorted_canonical_noERCC.sam
/share/samdata/dwyman/D12/TC_v1.0.7/sorted_canonical.sam
```

## PB65 cell BC017 (8 reads)
* m54284_180814_204240/72352410/ccs is an ISM transcript of Canx
* m54284_180814_002203/18677911/ccs is a prefix ISM of Canx
* m54284_180814_002203/19268005/ccs is a genomic transcript of Tcf3
* m54284_180814_002203/18809472/ccs is a suffix ISM of Tcf3
* m54284_180814_002203/49414590/ccs is NIC transcript of Drg1
* m54284_180814_002203/40042763/ccs is an FSM of Drg1
* m54284_180814_002203/8126905/ccs is antisense and overlaps Grb10
* m54284_180814_204240/59310495/ccs is intergenic

## PB65 cell BC018 (3 reads)
* m54284_180814_204240/9633944/ccs is a prefix ISM transcript of Canx
* m54284_180814_204240/44761238/ccs is an FSM transcript of Tcf3
* m54284_180814_002203/49414590/ccs is copied from BC017 (NIC transcript of Drg1) to test that is gets assigned to the same thing this time

## D12 C2C12 (5 reads)
* m54284_181114_004841/51511668/3202_53_CCS is a NNC transcript of Canx
* m54284_181114_004841/17433269/26_2129_CCS maps to the final exon of Canx but is antisense
* m54284_180829_231707/19726894/30_2449_CCS is FSM with novel 3p of Canx.
* m54284_181114_004841/11862379/29_4135_CCS is an NNC of Canx.
* m54284_181114_004841/71172442/1255_244_CCS maps to the final exon of Canx. Right now, TALON treats it as genomic, but I can see how we might want to call it an ISM at some point. 

