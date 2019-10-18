import pytest
from talon import talon, init_refs
import sqlite3
from .helper_fns import get_db_cursor
@pytest.mark.unit

class TestTempNovelGeneTab(object):
    def test_all(self):
        """ Get all genes in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        init_refs.make_temp_novel_gene_table(cursor, build)

        # Count number of entries
        cursor.execute(""" SELECT gene_ID FROM temp_gene """)
        results = [ x[0]  for x in cursor.fetchall()]
        conn.close()
        assert results == [1, 2, 3, 4, 5, 6]

    def test_empty_interval(self):
        """ The specified interval contains no monoexonic transcripts """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        init_refs.make_temp_novel_gene_table(cursor, build, chrom = "chr1",
                                                    start = 10000, end = 20000)
         # Count number of entries
        cursor.execute(""" SELECT gene_ID FROM temp_gene """)
        assert cursor.fetchall() == []

        conn.close()

    def test_non_empty(self):
        """ The specified interval contains one monoexonic transcript """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        init_refs.make_temp_novel_gene_table(cursor, build, chrom = "chr1",
                                                    start = 500, end = 1500)
        cursor.execute(""" SELECT gene_ID FROM temp_gene """)
        results = [ x[0]  for x in cursor.fetchall()]
        conn.close()
        assert results == [1, 2]
