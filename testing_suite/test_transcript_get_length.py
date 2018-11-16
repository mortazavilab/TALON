import pytest
import sys
sys.path.append("..")
from transcript import *
from edge import *
@pytest.mark.unit

class TestAddExonToTranscript(object):
    def test1(self):
        """ Create transcript with three exons (covered in other tests) and 
            compute the transcript length (sum of exon lengths). Correct
            answer is 11 + 11 + 6 = 28
        """
        start = 10
        end = 50
        e1_start = 10
        e1_end = 20
        e2_start = 30
        e2_end = 40
        e3_start = 45
        e3_end = 50

        t = Transcript("t1", "chr1", start, end, "+", "gene1", None)
        e1 = Edge("e1", "chr1", e1_start, e1_end, "+", "gene1", "t1", None)
        e2 = Edge("e2", "chr1", e2_start, e2_end, "+", "gene1", "t1", None) 
        e3 = Edge("e3", "chr1", e3_start, e3_end, "+", "gene1", "t1", None)

        t.add_exon(e1)
        t.add_exon(e2)
        t.add_exon(e3)
      
        assert t.get_length() == 28

