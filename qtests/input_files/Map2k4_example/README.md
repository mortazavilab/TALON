# Map2k4 test case
A transcript mapped to the Map2k4 locus is causing TALON to crash. It turns out that the transcript looks a lot like an FSM (visually and to TALON) because it has only known edges. However, there is no match for it in the annotation because it in fact combines attributes of different annotated transcripts- it is an NIC transcript.

```
grep "Map2k4" /bio/dwyman/pacbio_f2016/data/GENCODEvM7/gencode.vM7.annotation.gtf > Map2k4.gtf
```
