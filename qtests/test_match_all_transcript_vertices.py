import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
import dstruct
@pytest.mark.dbunit

class TestMatchAllVertices(object):

    def test_all_known_locations(self):
        """ Example where the toy transcript database contains matches for all
            vertices.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")
        orig_vertex_count = run_info['vertex']
        conn.close()

        chrom = "chr1"
        pos = [1, 100, 500, 600, 900, 1000]
        vertex_IDs = talon.match_all_transcript_vertices(chrom, pos, 
                                                        location_dict, run_info)

        # Make sure that no match got returned
        # TODO: 7 should really be 5 once I update the db schema to remove strand
        assert vertex_IDs == [ 2, 3, 4, 7 ] 
        assert run_info['vertex'] == orig_vertex_count

    def test_with_novel_location(self):
        """ Example where the toy transcript database contains a novel position.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")
        orig_vertex_count = run_info['vertex']
        orig_n_locations = len(location_dict)
        conn.close()

        chrom = "chr1"
        pos = [1, 150, 500, 600, 900, 1000]
        vertex_IDs = talon.match_all_transcript_vertices(chrom, pos,
                                                        location_dict, run_info)

        # Make sure that no match got returned
        # TODO: 7 should really be 5 once I update the db schema to remove strand
        new_vertex_count = run_info['vertex']
        assert vertex_IDs == [ "TALON-%d" % new_vertex_count, 3, 4, 7 ]
       
        # Make sure the data structures got updated
        assert new_vertex_count == orig_vertex_count + 1
        assert len(location_dict) == orig_n_locations + 1

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor

