import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
from helper_fns import *
@pytest.mark.dbunit

class TestSearchForISM(object):

    def test_find_no_match(self):
        """ Example where the toy transcript database contains no matches
            for the edge set.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)
        conn.close()

        edges = ( 5, 4, 3, 2 )
        gene_ID, matches = talon.search_for_ISM(edges, transcript_dict)

        # Make sure that no match got returned
        assert gene_ID == None

    def test_find_match(self):
        """ Example where the toy transcript database contains exactly one  
            ISM match for the transcript.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)

        edges = ( 2, 3 )
        gene_ID, matches = talon.search_for_ISM(edges, transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)

        assert gene_ID == correct_gene_ID
        conn.close()

