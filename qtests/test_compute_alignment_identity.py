import pytest
import sys
sys.path.append("..")
import transcript_utils as tu
@pytest.mark.unit

def test_compute_alignment_coverage():
    """ Computes what fraction of the read is actually aligned to
        the genome by excluding hard or soft-clipped bases."""
    
    MD = "MD:Z:10A3T0T10"
    SEQ = "GGGGGGGGGGTGGGAATGGGGGGGGGG"
    assert tu.compute_alignment_identity(MD, SEQ) == 23.0/27

    MD = "MD:Z:56^A45"
    SEQ = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    assert tu.compute_alignment_identity(MD, SEQ) == 101.0/102

    MD = "MD:Z:6G4C20G1A5C5A1^C3A15G1G15"
    SEQ = "GGGGGGAGGGGAGGGGGGGGGGGGGGGGGGGCGAGTGGGGGAGGGGGCGGGGTGGGGGGGGGGGGGGGAGAGGGGGGGGGGGGGGG" 

    align_length = len(SEQ) + 1 # incremented by 1 because of single deletion
    assert tu.compute_alignment_identity(MD, SEQ) == 76.0/align_length
