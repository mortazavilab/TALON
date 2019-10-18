import pytest
from talon import talon, init_refs
from .helper_fns import get_db_cursor
import sqlite3
@pytest.mark.unit

class TestEdgeDict(object):
    def test_all_edges(self):
        """ Get all edges in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        e_dict = init_refs.make_edge_dict(cursor, build)
     
        conn.close()
        assert len(e_dict.keys()) == 31

    def test_interval_edges(self):
        """ Get all edges in the database within the specified interval """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        e_dict = init_refs.make_edge_dict(cursor, build, chrom = "chr1",
                                          start = 1, end = 1000)

        conn.close()
        assert len(e_dict.keys()) == 6
        assert set(e_dict.keys()) == set([(4, 5, 'intron'),
                                          (2, 3, 'intron'),
                                          (5, 6, 'exon'),
                                          (3, 4, 'exon'),
                                          (1, 2, 'exon'),
                                          (6, 5, 'exon')])


