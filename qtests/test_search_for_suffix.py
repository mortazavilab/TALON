import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
from helper_fns import *
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

        edges = ( 1, 2, 3 )
        gene_ID, transcripts = talon.search_for_transcript_suffix(edges, transcript_dict)

        # Make sure that no match got returned
        assert gene_ID == None

    def test_find_match(self):
        """ Example where the toy transcript database contains exactly one 
            suffix match for the transcript.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)

        edges = ( 14, 15, 16 )
        gene_ID, transcripts = talon.search_for_transcript_suffix(edges, transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        assert gene_ID == correct_gene_ID
        assert transcripts == [(12, 13, 14, 15, 16)]

        conn.close()

    def test_find_match_with_diff_3(self):
        """ Example of an ISM where the final exon is different from the match
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)

        edges = ( 14, 15, 200 )
        gene_ID, transcripts = talon.search_for_transcript_suffix(edges, transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        assert gene_ID  == correct_gene_ID

        conn.close()

