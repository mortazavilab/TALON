import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
@pytest.mark.dbunit

class TestSearchForVertex(object):

    def test_find_no_match(self):
        """ Example where the toy transcript database contains no matches
            for the position.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)

        chrom = "chr1"
        pos = 0
        vertex_match = talon.search_for_vertex_at_pos(chrom, pos, location_dict)
        conn.close()

        # Make sure that exactly one match got returned
        assert vertex_match == None

    def test_find_exactly_one_match(self):
        """ Example where the toy transcript database contains exactly one match 
            for the position.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)

        chrom = "chr1"
        pos = 1
        match = talon.search_for_vertex_at_pos(chrom, pos, location_dict)
        conn.close()

        # Make sure that match is correct and that we can access various 
        # attributes using their names
        assert match["genome_build"] == "toy_build"
        assert match["chromosome"] == "chr1"
        assert match["position"] == 1

    #def test_find_two_matches(self):
    #    """ Example where the toy transcript database contains two matches
    #        for the position.
    #    """
    #    conn, cursor = get_db_cursor()
    #    build = "toy_build"
    #    chrom = "chr1"
    #    pos = 1000
    #    vertex_matches = talon.search_for_vertex_at_pos(build, chrom, pos, cursor)
    #    conn.close()
#
 #       # Make sure that exactly two matches got returned
 #       assert len(vertex_matches) == 2

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor

