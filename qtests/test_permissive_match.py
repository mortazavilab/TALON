import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
import dstruct
from helper_fns import *
@pytest.mark.dbunit

class TestPermissiveMatch(object):

    def test_edgecase_single_base_exon(self):
        """ Example where the first exon is only one basepair long
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")

        chrom = "chr1"
        pos = [1, 1, 500, 600]
        start = pos[0]
        splice_pos = pos[2]
        cutoff = 500
        strand = "+"
        
        vertex_match, diff = talon.permissive_vertex_search(chrom, start, 
                                                            strand, splice_pos,
                                                            "start", 
                                                            location_dict, 
                                                            cutoff, run_info)
        assert vertex_match == fetch_correct_vertex_ID(chrom, 1, cursor)
        assert diff == 0
        conn.close()

    def test_match_at_cutoff_distance(self):
        """ Example where the correct match is exactly the length of the cutoff
            away from the position """
        
        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")

        chrom = "chr1"
        pos = [1750, 1500, 1000, 900]
        start = pos[0]
        splice_pos = pos[1]
        cutoff = 250
        strand = "-"

        vertex_match, diff = talon.permissive_vertex_search(chrom, start,
                                                            strand, splice_pos,
                                                            "start",
                                                            location_dict,
                                                            cutoff, run_info)

        assert vertex_match == fetch_correct_vertex_ID(chrom, 2000, cursor)
        assert diff == -250
        conn.close()

    def test_beyond_cutoff_distance(self):
        """ Example where the only nearby vertices are beyond the cutoff 
            distance, prompting creation of a new vertex."""

        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")
        orig_vertex_count = run_info['vertex']
        conn.close()

        chrom = "chr1"
        pos = [1700, 1500, 1000, 900]
        start = pos[0]
        splice_pos = pos[1]
        cutoff = 250
        strand = "-"

        vertex_match, diff = talon.permissive_vertex_search(chrom, start,
                                                            strand, splice_pos,
                                                            "start",
                                                            location_dict,
                                                            cutoff, run_info)

        new_vertex_count = run_info['vertex']
        assert vertex_match == "TALON-%d" % new_vertex_count
        assert diff == "NA"

    def test_match_monoexonic(self):
        """ Test the permissive match strategy on a monoexonic transcript """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")

        chrom = "chr2"
        pos = [920, 970]
        start = pos[0]
        splice_pos = pos[1]
        cutoff = 500
        strand = "+"

        start_match, start_diff = talon.permissive_vertex_search(chrom, start,
                                                            strand, splice_pos,
                                                            "start",
                                                            location_dict,
                                                            cutoff, run_info)
        end = pos[1]
        splice_pos = pos[0]
        end_match, end_diff = talon.permissive_vertex_search(chrom, end,
                                                            strand, splice_pos,
                                                            "end",
                                                            location_dict,
                                                            cutoff, run_info)

        assert start_match == fetch_correct_vertex_ID(chrom, 900, cursor)
        assert start_diff == -20
        assert end_match == fetch_correct_vertex_ID(chrom, 1000, cursor)
        assert end_diff == 30
        conn.close()    

