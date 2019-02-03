import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
import dstruct
@pytest.mark.dbunit

class TestMatchAllEdges(object):

    def test_all_known_edges(self):
        """ Example where the toy transcript database contains matches for all
            vertices.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")
        conn.close()

        chrom = "chr2"
        vertex_IDs = [ 11, 12, 13, 14, 15, 16]
        strand = "+"
        edge_IDs = talon.match_all_transcript_edges(vertex_IDs, strand,
                                                        edge_dict, run_info)

        assert edge_IDs == [ 9, 10, 11, 12, 13 ] 

    def test_antisense(self):
        """ Example where all of the vertices are in the database, but the edges
            are not, because they are antisense to the original transcript """
        
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")
        orig_n_edges = len(edge_dict)
        conn.close()

        chrom = "chr2"
        vertex_IDs = [ 16, 15, 14, 13, 12, 11]
        strand = "-"

        edge_IDs = talon.match_all_transcript_edges(vertex_IDs, strand,
                                                        edge_dict, run_info)
        expected_edges = []
        for i in range(1,6):
            num = orig_n_edges + i
            edge_id = "TALON-%d" % num
            expected_edges.append(edge_id)

        assert edge_IDs == expected_edges

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor

