import pytest
import sys
import sqlite3
from helper_fns import *
sys.path.append("..")
import length_utils as lu
@pytest.mark.unit

class TestComputeExonLens(object):

    def test_compute_exon_len_plus(self):
        """ Plus strand example: chr1:1-100 """
        conn, cursor = get_db_cursor()
        build = "toy_build" 

        exon_lens = lu.get_all_exon_lengths(cursor, build)
        assert exon_lens[1] == 100

        conn.close()

    def test_compute_exon_len_minus(self):
        """ Minus strand example: chr1:2000-1500 """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        exon_lens = lu.get_all_exon_lengths(cursor, build)
        assert exon_lens[6] == 501

        conn.close()
