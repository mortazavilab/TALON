import pytest
import pysam
from talon import talon

def test_parse_custom_SAM_tags():
    """ Test that custom SAM tags are handled as expected """
   
    sam_file = "input_files/test_parse_custom_SAM_tags/toy_reads.sam"
    with pysam.AlignmentFile(sam_file, "rb") as sam: 
        for sam_record in sam:
            fraction_As, custom_label, allelic_label, \
            start_support, end_support = talon.parse_custom_SAM_tags(sam_record)
            if sam_record.query_name == "read_1":
                assert round(fraction_As,1) == 0.2
                assert custom_label == "yes"
                assert allelic_label == "paternal"
                assert start_support == "yes"
                assert end_support == "no"
            elif sam_record.query_name == "read_4":
                assert fraction_As == custom_label == allelic_label == None
                assert start_support == end_support == None
            else:
                pytest.fail("Did not recognize read name")
