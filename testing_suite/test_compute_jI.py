import pytest
from talon import transcript_utils as tu
@pytest.mark.unit

class TestComputejI(object):

    def test_compute_jI_example1(self):
        """ This example (from transcript c3098/f3p2/3199 in
            GM12878_chr1_clean.sam) contains 8 introns and 9 exons."""

        jI = "jI:B:i,1084384,1084480,1084507,1085877,1086013," + \
                   "1087138,1087205,1087501,1087598,1090352,1090429," + \
                   "1091471,1091566,1091990,1092104,1116059"
        start = 1081827
        cigar = "2557M97N26M1371N135M1126N66M297N96M2755N" + \
                "76M1043N94M425N113M23956N38M"
        assert tu.compute_jI(start, cigar) == jI

