import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
from helper_fns import *
@pytest.mark.integration

class TestIdentifyNIC(object):

    def test_NIC_match(self):
        """ Example where the transcript is an NIC match to an existing one by 
            virtue of skipping an exon.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build, "TALON")
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex2gene = talon.make_vertex_2_gene_dict(cursor)

        edge_IDs = (1, 200, 5)
        vertex_IDs = (1, 2, 5, 6)
        strand = "+"
        v_novelty = (0, 0, 0, 0)

        gene_ID, transcript_ID, novelty = talon.process_NIC(edge_IDs, vertex_IDs,
                                                      v_novelty, strand, transcript_dict,
                                                      vertex2gene, run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        assert gene_ID == correct_gene_ID
        conn.close()
