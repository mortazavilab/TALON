import pytest
import sys
import sqlite3
sys.path.append("..")
import talon as talon
from helper_fns import *
@pytest.mark.integration

class TestIdentifyNNC(object):

    def test_NNC_match(self):
        """ Example where the transcript is an NNC match to an existing one by
            virtue of a new splice donor.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)
        gene_starts, gene_ends = talon.make_gene_start_and_end_dict(cursor, build)

        chrom = "chr1"
        positions = [ 1, 110, 900, 1000]
        edge_IDs = [ run_info.edge + 1 ]
        vertex_IDs = [ run_info.vertex + 1, 5 ]
        strand = "+"
        v_novelty = [0, 0]

        gene_ID, transcript_ID, transcript_novelty, start_end_info = talon.process_NNC(chrom,
                                                            positions,
                                                            strand, edge_IDs,
                                                            vertex_IDs, transcript_dict,
                                                            gene_starts, gene_ends,
                                                            edge_dict, location_dict,
                                                            vertex_2_gene, run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        assert gene_ID == correct_gene_ID
        assert start_end_info["vertex_IDs"] == [1] + vertex_IDs + [6]
        assert transcript_dict[frozenset(start_end_info["edge_IDs"])] != None
        conn.close()
