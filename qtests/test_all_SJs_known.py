import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
@pytest.mark.dbunit

class TestAllSJsKnown(object):

    def test_find_true(self):
        """ Example where all of the introns are known.
        """
        novelty = [0, 0, 0, 0, 0,0 ,0]      

        # Make sure that no match got returned
        assert talon.all_SJs_known(novelty) == True

    def test_find_true_with_novel_exons(self):
        """ Example where all of the introns are known, but the exons are not.
            Note: This is not realistic biologically.
        """
        novelty = [1, 0, 1, 0, 1, 0 ,1]

        # Make sure that no match got returned
        assert talon.all_SJs_known(novelty) == True
   
    def test_find_false(self):
        """ Example with novel introns
        """
        novelty = [0, 1, 0, 1, 0, 1, 0]

        # Make sure that no match got returned
        assert talon.all_SJs_known(novelty) == False 


