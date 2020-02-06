import pysam
from talon import talon_label_reads as tlr

def test_plus_strand():
    """ Plus-strand transcript starts at pos 4 and is 10 bp long. The correct
        end position is 13 """

    read_file = "test_inputs/plus_strand_read.sam"
    with pysam.AlignmentFile(read_file) as sam:
        for record in sam:
            assert tlr.compute_transcript_end(record) == 13
            break

def test_minus_strand():
    """ Minus-strand transcript has SAM start = 4 and is 10 bp long. The correct
        end position is 4 """

    read_file = "test_inputs/minus_strand_read.sam"
    with pysam.AlignmentFile(read_file) as sam:
        for record in sam:
            assert tlr.compute_transcript_end(record) == 4
            break
