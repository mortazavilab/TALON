import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
from helper_fns import *
@pytest.mark.dbunit

class TestSearchForPrefix(object):

    def test_find_no_match(self):
        """ Example where the toy transcript database contains no matches
            for the edge set, because the edges are actually a suffix.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)
        conn.close()

        edges = ( 14, 15, 16 )
        gene_ID, transcripts = talon.search_for_transcript_prefix(edges, transcript_dict)

        # Make sure that no match got returned
        assert gene_ID == None
        conn.close()

    def test_find_match(self):
        """ Example where the toy transcript database contains exactly one 
            prefix match for the transcript.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)

        edges = ( 12, 13, 14 )
        gene_ID, transcripts = talon.search_for_transcript_prefix(edges, transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        assert gene_ID == correct_gene_ID
        assert transcripts == [(12, 13, 14, 15, 16)]
        conn.close()

    def test_find_match_with_diff_5(self):
        """ Example of an ISM where the start is different from the match
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)

        edges = ( 200, 13, 14 )
        gene_ID, transcripts = talon.search_for_transcript_prefix(edges, transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        assert gene_ID == correct_gene_ID

        conn.close()

    def test_find_monoexon_match(self):
        """ Input is a single exon that matches the start of an existing transcript
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)

        edges = ( 12, )
        gene_ID, matches = talon.search_for_transcript_prefix(edges, 
                                                              transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)

        assert gene_ID == correct_gene_ID
        conn.close()
