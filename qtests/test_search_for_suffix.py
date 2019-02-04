import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
@pytest.mark.dbunit

class TestSearchForSuffix(object):

    def test_find_no_match(self):
        """ Example where the toy transcript database contains no matches
            for the edge set.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)
        conn.close()

        edges = [ 9, 11 ]
        gene_ID = talon.search_for_transcript_suffix(edges, transcript_dict)

        # Make sure that no match got returned
        assert gene_ID == None

    def test_find_match(self):
        """ Example where the toy transcript database contains exactly one 
            suffix match for the transcript.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)
        conn.close()

        edges = [ 11, 12, 13 ]
        gene_ID = talon.search_for_transcript_suffix(edges, transcript_dict)

        # Make sure that correct match got returned
        assert gene_ID == 3

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor

