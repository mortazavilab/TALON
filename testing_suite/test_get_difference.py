import pytest
import sys
sys.path.append("..")
from transcript_match_tracker import *
@pytest.mark.unit

class TestGetDifference(object):
    def test_1(self):
        """ Example where intervals match exactly. So the answer
            should be diff_3 = diff_5 = 0. 
        """
        a = [ 10, 20 ] 
        b = [ 10, 20 ]
        strand = "+"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == 0 and diff_3 == 0

    def test_2(self):
        """ Example where intervals match exactly. So the answer
            should be diff_3 = diff_5 = 0.
        """
        a = [ 10, 20 ]
        b = [ 10, 20 ]
        strand = "-"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == 0 and diff_3 == 0

    def test_3(self):
        """ Example where a starts upstream of b (+). So the answer
            should be diff_3 = 0 and diff_5 = -5.
        """
        a = [ 5, 20 ]
        b = [ 10, 20 ]
        strand = "+"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == -5 and diff_3 == 0

    def test_4(self):
        """ Example where a starts upstream of b (-). So the answer
            should be diff_3 = 0 and diff_5 = -5.
        """
        a = [ 10, 25 ]
        b = [ 10, 20 ]
        strand = "-"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == -5 and diff_3 == 0

    def test_5(self):
        """ Example where a starts later than b (+). So the answer
            should be diff_3 = 0 and diff_5 = 5.
        """
        a = [ 15, 20 ]
        b = [ 10, 20 ]
        strand = "+"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == 5 and diff_3 == 0

    def test_6(self):
        """ Example where a starts downstream of b (-). So the answer
            should be diff_3 = 0 and diff_5 = 5.
        """
        a = [ 10, 15 ]
        b = [ 10, 20 ]
        strand = "-"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == 5 and diff_3 == 0

    def test_7(self):
        """ Example where a ends downstream of b (+). So the answer
            should be diff_3 = 5 and diff_5 = 0.
        """
        a = [ 10, 25 ]
        b = [ 10, 20 ]
        strand = "+"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == 0 and diff_3 == 5

    def test_8(self):
        """ Example where a ends downstream of b (-). So the answer
            should be diff_3 = 5 and diff_5 = 0.
        """
        a = [ 5, 20 ]
        b = [ 10, 20 ]
        strand = "-"
        diff_5, diff_3 = get_difference(a, b, strand)
        assert diff_5 == 0 and diff_3 == 5
