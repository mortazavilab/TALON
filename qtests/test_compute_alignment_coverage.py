import pytest
import sys
sys.path.append("..")
import transcript_utils as tu
@pytest.mark.unit

def test_compute_alignment_coverage():
    """ Computes what fraction of the read is actually aligned to
        the genome by excluding hard or soft-clipped bases."""
    
    cigar = "5S90M5H"
    assert tu.compute_alignment_coverage(cigar) == 0.9
