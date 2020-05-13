import pytest
from talon import transcript_utils as tu
@pytest.mark.unit

def test_compute_alignment_coverage():
    """ Computes what fraction of the read is actually aligned to
        the genome by excluding hard or soft-clipped bases."""
    
    cigar = "5S90M5H"
    assert tu.compute_alignment_coverage(cigar) == 0.9

    cigar = "5S90=5H"
    assert tu.compute_alignment_coverage(cigar) == 0.9


def test_compute_coverage_with_introns():
    """ Make sure that intron portion of alignment is not being counted """

    cigar = "5S45M1000N45=5H"
    assert tu.compute_alignment_coverage(cigar) == 0.9
