import pytest
from talon import talon, init_refs
from .helper_fns import fetch_correct_ID, get_db_cursor
@pytest.mark.integration

class TestIdentifyNNC(object):

    def test_NNC_match(self):
        """ Example where the transcript is an NNC match to an existing one by
            virtue of a new splice donor.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        init_refs.make_temp_novel_gene_table(cursor, build)
        init_refs.make_temp_transcript_table(cursor, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr1"
        positions = [ 1, 110, 900, 1000]
        edge_IDs = [ talon.edge_counter.value() + 1 ]
        vertex_IDs = [ talon.vertex_counter.value() + 1, 5 ]
        strand = "+"
        v_novelty = [0, 0]

        gene_ID, transcript_ID, transcript_novelty, start_end_info, fusion = talon.process_NNC(chrom,
                                                            positions,
                                                            strand, edge_IDs,
                                                            vertex_IDs, transcript_dict,
                                                            gene_starts, gene_ends,
                                                            edge_dict, location_dict,
                                                            vertex_2_gene, run_info,
                                                            cursor, "temp_gene",
                                                            "temp_transcript")

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        assert gene_ID == correct_gene_ID
        assert start_end_info["vertex_IDs"] == [1] + vertex_IDs + [6]
        assert transcript_dict[frozenset(start_end_info["edge_IDs"])] != None
        assert fusion == False
        conn.close()
