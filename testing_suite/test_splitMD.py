import pytest
import sys
sys.path.append("..")
import transcript_utils as tu
@pytest.mark.unit

class TestSplitMD(object):
    def test_splitMD(self):
        """ Easy case- full match"""
    
        MD = "MD:Z:100"
        ops, cts = tu.splitMD(MD)
        assert ops == ["M"]
        assert cts == [100]

    def test_with_mismatches(self):
        """ MD tag with mismatches in it """
        MD = "MD:Z:48T42G8"
        ops, cts = tu.splitMD(MD)
        assert ops == ["M", "X", "M", "X", "M"]
        assert cts == [48, 1, 42, 1, 8]

    def test_with_deletion(self):
        """ MD tag with deletions in it """
        MD = "MD:Z:56^ACG45"
        ops, cts = tu.splitMD(MD)
        assert ops == ["M", "D", "M"]
        assert cts == [56, 3, 45]
    
    def test_with_deletions_and_mismatches(self):
        """ MD tag with consecutive deletions and mismatches in it """
        MD = "MD:Z:56^ACG0AT45"
        ops, cts = tu.splitMD(MD)
        assert ops == ["M", "D", "M", "X", "M"]
        assert cts == [56, 3, 0, 2, 45] 
