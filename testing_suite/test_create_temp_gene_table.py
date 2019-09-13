import pytest
from talon import talon

from .helper_fns import get_db_cursor


@pytest.mark.dbunit

class TestMakeTempTable(object):

    def test_create_temp_table(self):
        """ Create the table and make sure it is accessible even if it is 
            empty, and make sure it doesn't clash with the TALON database
        """
        # Open TALON database at the same time
        conn, cursor = get_db_cursor()
        build = "toy_build"

        # Now run the temp table creation
        talon.make_temp_novel_gene_table(cursor, build)
        try:
            query = """SELECT * FROM temp_gene"""
            cursor.execute(query)
            results = cursor.fetchall()
        except:
            pytest.fail("Something went wrong with temp table query")
        
        conn.close()
