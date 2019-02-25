import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
from helper_fns import *
@pytest.mark.dbunit

class TestSearchForTranscript(object):

    def test_find_no_match(self):
        """ Example where the toy transcript database contains no matches
            for the edge set.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)
        conn.close()

        edges = frozenset({ 1, 3, 4, 5 })
        gene_ID, transcript = talon.search_for_transcript(edges, 
                                                             transcript_dict)

        # Make sure that no match got returned
        assert gene_ID == None
        assert transcript == None

    def test_find_match(self):
        """ Example where the toy transcript database contains exactly one match 
            for the transcript.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)

        edges = frozenset({ 12, 13, 14, 15, 16 })
        gene_ID, transcript = talon.search_for_transcript(edges,
                                                             transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        correct_transcript_ID = fetch_correct_ID("TG2-001", "transcript", cursor)

        assert gene_ID == correct_gene_ID
        assert transcript[0] == correct_transcript_ID
        conn.close()

