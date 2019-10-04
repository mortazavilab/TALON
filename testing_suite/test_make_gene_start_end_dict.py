import pytest
from talon import talon
import sqlite3
@pytest.mark.unit

class TestGeneStartEnd(object):
    def test_all_genes(self):
        """ Get starts and ends for all known genes in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        starts, ends = talon.make_gene_start_and_end_dict(cursor, build)
     
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

        starts, ends = talon.make_gene_start_and_end_dict(cursor, build,
                                                          chrom = "chr1",
                                                          start = 1,
                                                          end = 1000)

        conn.close()

        assert starts == {1: {1: 1},
                          2: {2000: 8}}

        assert ends == {1: {1000: 6},
                        2: {900: 5}}

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor 
