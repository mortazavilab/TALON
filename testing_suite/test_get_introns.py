import pytest
import sys
sys.path.append("..")
from sam_transcript import *

class TestGetIntrons(object):

    def test_multiexon_without_jI(self):
        """ This example (from transcript c3098/f3p2/3199 in
            GM12878_chr1_clean.sam) contains 8 introns and 9 exons. The jI
            field is omitted from the input args, so the function must call
            compute_jI(start, cigar).
        """
        jI = [ 1084384,1084480,1084507,1085877,1086013, \
               1087138,1087205,1087501,1087598,1090352,1090429, \
               1091471,1091566,1091990,1092104,1116059 ]
        fields = [ "NH:i:1", "HI:i:1", "NM:i:0", "MD:Z:3201", \
                   "jM:B:c,22,22,22,22,22,22,22,22"]
        start = 1081827
        cigar = "2557M97N26M1371N135M1126N66M297N96M2755N" + \
                "76M1043N94M425N113M23956N38M"
        assert get_introns(fields, start, cigar) == jI

    def test_multiexon_with_jI(self):
        """ This example (from transcript c3098/f3p2/3199 in 
            GM12878_chr1_clean.sam) contains 8 introns and 9 exons.
            The jI field is supplied in the input args to the function, so
            all it has to do is find the field, split, and return it.
        """
        jI = [ 1084384,1084480,1084507,1085877,1086013, \
               1087138,1087205,1087501,1087598,1090352,1090429, \
               1091471,1091566,1091990,1092104,1116059 ]
        fields = [ "NH:i:1", "HI:i:1", "NM:i:0", "MD:Z:3201", \
                   "jM:B:c,22,22,22,22,22,22,22,22", jI ] 
        start = 1081827
        cigar = "2557M97N26M1371N135M1126N66M297N96M2755N" + \
                "76M1043N94M425N113M23956N38M"
        assert get_introns(fields, start, cigar) == jI

def test_compute_jI():
    """ This example (from transcript c3098/f3p2/3199 in
            GM12878_chr1_clean.sam) contains 8 introns and 9 exons."""

    jI = "jI:B:i,1084384,1084480,1084507,1085877,1086013," + \
                   "1087138,1087205,1087501,1087598,1090352,1090429," + \
                   "1091471,1091566,1091990,1092104,1116059"
    start = 1081827
    cigar = "2557M97N26M1371N135M1126N66M297N96M2755N" + \
            "76M1043N94M425N113M23956N38M"
    assert compute_jI(start, cigar) == jI

