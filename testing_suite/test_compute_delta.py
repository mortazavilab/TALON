import pytest
import sys
sys.path.append("..")
import talon as talon
@pytest.mark.unit

class TestComputeDelta(object):
    def test_1(self):
        """ Plus-strand upstream
        """
        a = 100 
        b = 96
        strand = "+"
        assert talon.compute_delta(a, b, strand) == -4

    def test_2(self):
        """ Plus-strand downstream
        """
        a = 96
        b = 100
        strand = "+"
        assert talon.compute_delta(a, b, strand) == 4

    def test_3(self):
        """ Minus-strand upstream
        """
        a = 96
        b = 100
        strand = "-"
        assert talon.compute_delta(a, b, strand) == -4

    def test_4(self):
        """ Minus-strand downstream
        """
        a = 100
        b = 96
        strand = "-"
        assert talon.compute_delta(a, b, strand) == 4

