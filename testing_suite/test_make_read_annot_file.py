import pytest
import subprocess
from talon.post import get_read_annotations

@pytest.mark.integration
class TestGetReadAnnot(object):
    """ Make sure that the post-TALON get_read_annotations utility is working """

    def test_fetch_reads(self):
        """ Test if the start and end positions returned for each read are
            correct. These positions and read lengths should not include
            soft-clipped regions. """
        database = "scratch/toy_mod.db"
        build = "toy_build"
        reads = get_read_annotations.fetch_reads(database, build)

        for r in reads:
            read_name, dataset, genome_build, gene_ID, \
            transcript_ID, chrom, read_start, read_end, \
            strand, n_exons, read_length, fraction_As, custom_label, \
            allelic_label = r
           
            if read_name == "read_1":
                assert read_start == 1 
                assert read_end == 1004
                assert read_length == 306
            elif read_name == "read_2":
                assert read_start == 999 
                assert read_end == 900
                assert read_length == 100
            elif read_name == "read_3":
                assert read_start == 1 
                assert read_end == 100
                assert read_length == 100
            elif read_name == "read_4":
                assert read_start == 500 
                assert read_end == 519
            elif read_name == "read_5":
                assert read_start == 109
                assert read_end == 100
            elif read_name == "read_6":
                assert read_start == 200
                assert read_end == 209
            else:
                pytest.fail("Unexpected read ID")

   
