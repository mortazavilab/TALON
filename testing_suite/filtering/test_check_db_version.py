import os
import pytest
from talon.post import filter_talon_transcripts as filt

class Test(object):

    def test_check_v4_db(self):
        """ check that db with v4 throws an error """

        database = "input_files/example_v4.db"
        

        with pytest.raises(ValueError) as errinfo:
            filt.check_db_version(database)
        assert 'Database version is not compatible with v5.0 filtering.' in str(errinfo.value)