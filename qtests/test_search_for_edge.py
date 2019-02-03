import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
@pytest.mark.dbunit
@pytest.mark.incremental

class TestSearchForEdge(object):

    def test_find_match(self):
        """ Example where the toy transcript edge dict does not contain the 
            edge being queried.
        """
        conn, cursor = get_db_cursor()
       
        # Create a location dict and then fetch vertices for two psotions
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        edge_dict = talon.make_edge_dict(cursor)
        conn.close()

        chrom = "chr1"
        pos1 = 600
        pos2 = 500
        v1 = talon.search_for_vertex_at_pos(chrom, pos1, location_dict)["location_ID"]
        v2 = talon.search_for_vertex_at_pos(chrom, pos2, location_dict)["location_ID"]        

        assert v1 != None
        assert v2 != None

        # Now look for the edge between them
        edge_match = talon.search_for_edge(v1, v2, "exon", edge_dict)
        assert edge_match == None
 
        # Try them in the opposite order
        edge_match = talon.search_for_edge(v2, v1, "exon", edge_dict)
        assert edge_match["edge_ID"] ==  3

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor

