import pytest
import sys
sys.path.append("..")
import transcript_utils as tu
@pytest.mark.unit

def test_splitMD():
    """ Make sure that the splitMD function is working correctly"""
    
    MD = "MD:Z:100"
    ops, cts = tu.splitMD(MD)
    assert ops == ["M"]
    assert cts == [100]
