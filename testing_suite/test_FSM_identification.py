import pytest
import sys
import sqlite3
sys.path.append("..")
import talon as talon
from helper_fns import *
@pytest.mark.integration

class TestIdentifyFSM(object):

    def test_FSM_perfect(self):
        """ Example where the transcript is a perfect full splice match.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)        
        transcript_dict = talon.make_transcript_dict(cursor, build)

        chrom = "chr1"
        positions = [1, 100, 500, 600, 900, 1010]
        strand = "+"
        edge_IDs = [2, 3, 4]
        vertex_IDs = [2, 3, 4, 5]
        v_novelty = [0, 0, 0, 0]

        gene_ID, transcript_ID, novelty, start_end_info = talon.process_FSM(chrom, 
                                                            positions, strand, 
                                                            edge_IDs, vertex_IDs, 
                                                            transcript_dict, 
                                                            edge_dict,
                                                            location_dict, run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor) 
        correct_transcript_ID = fetch_correct_ID("TG1-001", "transcript", cursor)
        assert gene_ID == correct_gene_ID
        assert transcript_ID == correct_transcript_ID
        assert novelty == []
        assert start_end_info["start_vertex"] == 1
        assert start_end_info["end_vertex"] == 6
        assert start_end_info["diff_3p"] == 10
        conn.close()

    def test_FSM_end_diff(self):
        """ Example where the transcript is an FSM but has a difference on
            the ends large enough to be novel.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        orig_vertices = run_info['vertex']

        chrom = "chr2"
        positions = [1, 100, 500, 600, 900, 1301] #Last postion is > 300bp away
        strand = "+"
        edge_IDs = [13, 14, 15]
        vertex_IDs = [14, 15, 16, 17] 
        v_novelty = [0, 0, 0, 0]
   
        gene_ID, transcript_ID, novelty, start_end_info = talon.process_FSM(chrom,
                                                            positions, strand,
                                                            edge_IDs, vertex_IDs,
                                                            transcript_dict,
                                                            edge_dict,
                                                            location_dict, run_info)

        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        correct_transcript_ID = fetch_correct_ID("TG2-001", "transcript", cursor)
        assert gene_ID == correct_gene_ID
        assert transcript_ID == correct_transcript_ID
        assert start_end_info["end_vertex"] == orig_vertices + 1
        conn.close()

    def test_FSM_start_diff(self):
        """ Example where the transcript is an FSM but has a difference on
            the start large enough to be novel.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        orig_vertices = run_info['vertex']

        chrom = "chr1"
        positions = [2501, 1500, 1000, 900] #First postion is > 500bp away
        strand = "-"
        edge_IDs = [7]
        vertex_IDs = [7, 6]
        v_novelty = [0, 0]

        gene_ID, transcript_ID, novelty, start_end_info = talon.process_FSM(chrom,
                                                            positions, strand,
                                                            edge_IDs, vertex_IDs,
                                                            transcript_dict,
                                                            edge_dict,
                                                            location_dict, run_info)

        correct_gene_ID = fetch_correct_ID("TG3", "gene", cursor)
        correct_transcript_ID = fetch_correct_ID("TG3-001", "transcript", cursor)
        assert gene_ID == correct_gene_ID
        assert transcript_ID == correct_transcript_ID
        assert start_end_info["start_vertex"] == orig_vertices + 1
        assert start_end_info["end_vertex"] == 5
        conn.close() 

#    def test_no_match(self):
#        """ Example with no FSM or ISM match """

       # conn, cursor = get_db_cursor()
      #  build = "toy_build"
     #   edge_dict = talon.make_edge_dict(cursor)
    #    location_dict = talon.make_location_dict(build, cursor)
   #     run_info = talon.init_run_info(cursor, build)
  #      transcript_dict = talon.make_transcript_dict(cursor, build)
 #       gene_starts, gene_ends = talon.make_gene_start_and_end_dict(cursor)

#        edge_IDs = [ 200 ] 
#        vertex_IDs = [2, 5 ]
#        v_novelty = [ 0, 0 ]

#        gene_ID, transcript_ID, novelty = talon.process_FSM(edge_IDs, vertex_IDs,
#                                                            transcript_dict, gene_starts, gene_ends,
#                                                            run_info)
#        assert gene_ID == transcript_ID == None 
#        conn.close()       

