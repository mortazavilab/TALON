import pytest
import sys
sys.path.append("..")
from sam_transcript import *

class TestComputeTranscriptEnd(object):
    def test_MN_Only(self):
        """ This example (from transcript c3098/f3p2/3199 in 
            GM12878_chr1_clean.sam) contains only M and N operations. The 
            transcript is on the minus strand, but the sequence and CIGAR in 
            BAM are relative to the forward strand.
        """
        start = 1081827
        cigar = "2557M97N26M1371N135M1126N66M297N96M2755N" + \
                "76M1043N94M425N113M23956N38M"
        assert compute_transcript_end(start, cigar) == 1116097

    def test_deletion(self):
        """ This example (from transcript c2986/f14p15/2748 in 
            GM12878_chr1_clean.sam) contains M, N, and D operations. This is a 
            plus strand transcript. 
        """
        start = 203305518
        cigar = "231M1355N1013M1D1504M"

        assert compute_transcript_end(start, cigar) == 203309621

    def test_softClipping(self):
        """ This example (from transcript c28239/f1p4/3076 in
            GM12878_chr1_clean.sam) contains S operations. This is a
            plus strand transcript.
        """
        start = 167936402
        cigar = "136S114M15291N55M14767N93M8108N186M12479N114M3595N136M" + \
                "1886N215M9041N94M1294N120M543N261M33577N118M4536N116M" + \
                "1444N87M228N328M18393N139M1830N157M630N89M1892N106M6907N420M"
        
        assert compute_transcript_end(start, cigar) == 168075790

    def test_softClipping_atEnd(self):
        """ This example (from transcript m54284_180611_122744/21037936/24_2799_CCS in
            /share/samdata/dwyman/D1/TC_v1.0.5/sorted_canonical.sam) is on the
            minus strand and contains soft clipping at the end"""

        start = 152673671
        cigar = "36M2698N375M671N25M1580N54M242N121M1753N111M9625N6M45I3581N81M21372N133M1810S"

        assert compute_transcript_end(start, cigar) == 152716134 
