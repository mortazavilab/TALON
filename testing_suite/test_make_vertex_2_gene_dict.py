import pytest
from talon import talon
import sqlite3
@pytest.mark.unit

class TestVertex2GeneDict(object):
    def test_all_v2g(self):
        """ Get all vertices in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        v2g_dict = talon.make_vertex_2_gene_dict(cursor, build)
     
        conn.close()
        assert len(v2g_dict.keys()) == 34

    def test_interval_v2g(self):
        """ Get all vertices in the database within the specified interval """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        v2g_dict = talon.make_vertex_2_gene_dict(cursor, build, chrom = "chr1",
                                               start = 1, end = 1000)

        conn.close()
        assert len(v2g_dict.keys()) == 6


def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor 
