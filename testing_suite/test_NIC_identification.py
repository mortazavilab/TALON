import pytest
from talon import talon, init_refs
from .helper_fns import fetch_correct_ID, get_db_cursor
@pytest.mark.integration

class TestIdentifyNIC(object):

    def test_NIC_match(self):
        """ Example where the transcript is an NIC match to an existing one by
            virtue of skipping an exon.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr1"
        positions = [ 1, 100, 900, 1000]
        edge_IDs = [ talon.edge_counter.value() + 1 ]
        vertex_IDs = [ 2, 5 ]
        strand = "+"
        v_novelty = [0, 0]

        gene_ID, transcript_ID, novelty, start_end_info, fusion = talon.process_NIC(chrom,
                                                            positions,
                                                            strand, edge_IDs,
                                                            vertex_IDs, transcript_dict,
                                                            gene_starts, gene_ends,
                                                            edge_dict, location_dict,
                                                            vertex_2_gene, run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        assert gene_ID == correct_gene_ID
        assert start_end_info["vertex_IDs"] == [1,2,5,6]
        assert transcript_dict[frozenset(start_end_info["edge_IDs"])] != None
        assert fusion == False
        conn.close()

    def test_antisense(self):
        """ Example where the vertices are known but there is no same-strand
            match """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        edge_dict = init_refs.make_edge_dict(cursor)
        locations = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        # Construct temp novel gene db
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")

        chrom = "chr1"
        start = 1000
        end = 1
        edge_IDs = [ talon.edge_counter.value() + 1 ]
        positions = [ 1000, 900, 100, 1]
        vertex_IDs = [ 5, 2 ]
        strand = "-"
        anti_strand = "+"
        v_novelty = (0, 0, 0, 0)

        # Find antisense match
        gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info = \
                                      talon.process_spliced_antisense(chrom, positions,
                                                                  strand, edge_IDs,
                                                                  vertex_IDs,
                                                                  transcript_dict,
                                                                  gene_starts,
                                                                  gene_ends,
                                                                  edge_dict, locations,
                                                                  vertex_2_gene, run_info,
                                                                  cursor, "temp_gene")
        #anti_gene_ID = talon.find_gene_match_on_vertex_basis(vertex_IDs,
        #                                                     anti_strand,
        #                                                     vertex_2_gene)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        anti_gene_ID = gene_novelty[-1][-1]
        assert anti_gene_ID == correct_gene_ID
        assert start_end_info["vertex_IDs"] == [6, 5, 2, 1]

        conn.close()
