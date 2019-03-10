import pytest
import sys
import sqlite3
sys.path.append("..")
import talon as talon
@pytest.mark.dbunit

class TestAllExonsKnown(object):

    def test_find_true(self):
        """ Example where all of the exons are known.
        """
        # Remember that first pos is first intron, last is last intron
        novelty = [0, 0, 0, 0, 0 ]      

        # Make sure that no match got returned
        assert talon.check_all_exons_known(novelty) == True

    def test_find_true_with_novel_exons(self):
        """ Example where all of the exons are known, but the introns are not.
            Note: This is not necessarily realistic biologically.
        """
        novelty = [1, 0, 1, 0, 1]

        # Make sure that no match got returned
        assert talon.check_all_exons_known(novelty) == True
   
    def test_find_false(self):
        """ Example with novel exons
        """
        novelty = [0, 1, 0, 1, 0]

        # Make sure that no match got returned
        assert talon.check_all_exons_known(novelty) == False 

    def test_monoexonic(self):
        """ Monoexonic known exon """
        novelty = [0]

        assert talon.check_all_exons_known(novelty) == True
