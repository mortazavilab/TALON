from talon import talon_label_reads as tlr

def test_frac_as():
    """ Compute the fraction of As in the sequence, making sure we don't have
        int rounding """

    assert tlr.compute_frac_As("AAAAAA") == 1
    assert tlr.compute_frac_As("AATG") == 0.5
    assert tlr.compute_frac_As("ACTGACTGG") == 2.0/9.0
