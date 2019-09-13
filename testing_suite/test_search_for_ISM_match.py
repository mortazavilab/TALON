import pytest
from talon import talon
from .helper_fns import fetch_correct_ID, get_db_cursor
@pytest.mark.dbunit

class TestSearchForISM(object):

    def test_find_no_match(self):
        """ Example where the toy transcript database contains no matches
            for the edge set.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor, build)
        conn.close()

        edges = ( 100, 200, 300)
        matches = talon.search_for_ISM(edges, transcript_dict)

        # Make sure that no match got returned
        assert matches == None

    def test_find_match(self):
        """ Example where the toy transcript database contains exactly one  
            ISM match for the transcript.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor, build)

        edges = ( 2, 3 )
        matches = talon.search_for_ISM(edges, transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)

        assert matches[0]["gene_ID"] == correct_gene_ID
        conn.close()

    def test_find_monoexon_match(self):
        """ Input is a sinlge exon that matches part of an existing transcript
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor, build)

        edges = ( 14, )
        matches = talon.search_for_ISM(edges, transcript_dict)

        # Make sure that correct match got returned
        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)

        assert matches[0]["gene_ID"] == correct_gene_ID
        conn.close()

