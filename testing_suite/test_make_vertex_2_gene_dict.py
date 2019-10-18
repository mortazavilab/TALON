import pytest
from talon import talon, init_refs
import sqlite3
from .helper_fns import get_db_cursor
@pytest.mark.unit

class TestVertex2GeneDict(object):
    def test_all_v2g(self):
        """ Get all vertices in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        v2g_dict = init_refs.make_vertex_2_gene_dict(cursor, build)
     
        conn.close()
        assert len(v2g_dict.keys()) == 34

    def test_interval_v2g(self):
        """ Get all vertices in the database within the specified interval """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        v2g_dict = init_refs.make_vertex_2_gene_dict(cursor, build, chrom = "chr1",
                                               start = 1, end = 1000)

        conn.close()
        assert len(v2g_dict.keys()) == 6

