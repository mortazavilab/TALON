# Known genes for database test
The GTF file contains annotations for 4 genes from GENCODE v24:
* ACTB (- strand)
* LINC01128 (+ strand)
* RPS6KA1 (+ strand)
* XIST (- strand)


```
grep '"ACTB";' /bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf >> test.gtf
grep '"LINC01128";' /bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf >> test.gtf
grep '"RPS6KA1";' /bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf >> test.gtf
grep '"XIST";' /bio/dwyman/pacbio_f2016/data/GENCODEv24/gencode.v24.annotation.gtf >> test.gtf
```
