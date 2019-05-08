import pytest
import sys
import sqlite3
from helper_fns import *
sys.path.append("..")
import length_utils as lu
@pytest.mark.unit

class TestComputeTxLen(object):

    def test_compute_exon_len_plus(self):
        """ Plus strand example: chr1:1-100 """
        conn, cursor = get_db_cursor()
        build = "toy_build" 

        # Get exon lengths
        exon_lens = lu.get_all_exon_lengths(cursor, build)

        # Fetch transcript row
        cursor.execute("SELECT * FROM transcripts WHERE transcript_ID = '1'")
        transcript_row = cursor.fetchone()

        assert lu.get_transcript_length(transcript_row, exon_lens) == 100 + 101 + 101

        conn.close()

