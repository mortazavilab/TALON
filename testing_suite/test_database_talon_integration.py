import pytest
import sys
sys.path.append("..")

class test_database_talon_integration(object):
    """ Compared to a lot of the other tests in the suite, this one is intended 
        to make sure that the major parts of TALON (database init, dataset 
        addition over time) are working correctly. """

    def test_db_initialization(self):
        """ Initializes a TALON database from a very minimal GTF entry 
        """
              
        assert 

