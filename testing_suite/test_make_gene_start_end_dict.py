import pytest
from talon import talon, init_refs
import sqlite3
from .helper_fns import get_db_cursor
@pytest.mark.unit

class TestGeneStartEnd(object):
    def test_all_genes(self):
        """ Get starts and ends for all known genes in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        conn.close()

        # Check starts. Field 1: Gene ID, Field 2: start pos, Field 3: vertex ID 
        assert starts == {1: {1: 1}, 
                          2: {2000: 8}, 
                          3: {5000: 9}, 
                          4: {1: 13}, 
                          5: {1: 19, 800: 23}, 
                          6: {4000: 34}}

        assert ends == {1: {1000: 6}, 
                        2: {900: 5}, 
                        3: {6500: 12}, 
                        4: {1000: 18}, 
                        5: {2200: 30, 2600: 32}, 
                        6: {1000: 33}}

    def test_interval_1_1000(self):
        """ Get starts and ends for genes in the database that overlap the 
            interval chr1:1-1000 """
       
        conn, cursor = get_db_cursor()
        build = "toy_build"

        starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start",
                                                          chrom = "chr1",
                                                          start = 1,
                                                          end = 1000)
        ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end",
                                                          chrom = "chr1",
                                                          start = 1,
                                                          end = 1000)

        conn.close()

        assert starts == {1: {1: 1}}

        assert ends == {1: {1000: 6},
                        2: {900: 5}}

    def test_interval_1_2000(self):
        """ Get starts and ends for genes in the database that overlap the
            interval chr1:1-2000 """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start",
                                                          chrom = "chr1",
                                                          start = 1,
                                                          end = 2000)
        ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end",
                                                          chrom = "chr1",
                                                          start = 1,
                                                          end = 2000)

        conn.close()

        assert starts == {1: {1: 1},
                          2: {2000: 8}}

        assert ends == {1: {1000: 6},
                        2: {900: 5}}
