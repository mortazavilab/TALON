import pytest
from talon import talon, init_refs
import sqlite3
from .helper_fns import get_db_cursor
@pytest.mark.unit

class TestTempMonoexonTab(object):
    def test_all(self):
        """ Get all monoexonic in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        init_refs.make_temp_monoexonic_transcript_table(cursor, build)

        # Count number of entries
        cursor.execute(""" SELECT transcript_ID FROM temp_monoexon """)
        results = [ x[0]  for x in cursor.fetchall()]
        conn.close()
        assert results == [7]

    def test_empty_interval(self):
        """ The specified interval contains no monoexonic transcripts """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        init_refs.make_temp_monoexonic_transcript_table(cursor, build, chrom = "chr1",
                                                    start = 1, end = 1000)
         # Count number of entries
        cursor.execute(""" SELECT transcript_ID FROM temp_monoexon """)
        assert cursor.fetchall() == []

        conn.close()

    def test_non_empty(self):
        """ The specified interval contains one monoexonic transcript """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        init_refs.make_temp_monoexonic_transcript_table(cursor, build, chrom = "chr4",
                                                    start = 2000, end = 3000)
         # Count number of entries
        cursor.execute(""" SELECT transcript_ID FROM temp_monoexon """)
        results = [ x[0]  for x in cursor.fetchall()]
        conn.close()
        assert results == [7]
