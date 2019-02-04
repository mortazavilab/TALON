import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
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

        edges = [ 9, 11, 12, 13 ]
        gene_ID, transcript_ID = talon.search_for_transcript(edges, 
                                                             transcript_dict)

        # Make sure that no match got returned
        assert gene_ID == None
        assert transcript_ID == None

    def test_find_match(self):
        """ Example where the toy transcript database contains exactly one match 
            for the transcript.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = talon.make_transcript_dict(cursor)
        conn.close()

        edges = [ 9, 10, 11, 12, 13 ]
        gene_ID, transcript_ID = talon.search_for_transcript(edges,
                                                             transcript_dict)

        # Make sure that correct match got returned
        assert gene_ID == 3
        assert transcript_ID == 3

    #    conn, cursor = get_db_cursor()
    #    build = "toy_build"
    #    location_dict = talon.make_location_dict(build, cursor)

    #    chrom = "chr1"
    #    pos = 1
    #    match = talon.search_for_vertex_at_pos(chrom, pos, location_dict)
    #    conn.close()

    #    # Make sure that match is correct and that we can access various 
    #    # attributes using their names
    #    assert match["genome_build"] == "toy_build"
    #    assert match["chromosome"] == "chr1"
    #    assert match["position"] == 1

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor

