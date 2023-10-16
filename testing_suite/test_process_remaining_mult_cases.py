import pytest
from talon import talon, init_refs
from .helper_fns import  fetch_correct_ID, get_db_cursor
import logging

logging.basicConfig(level=logging.DEBUG)

@pytest.mark.integration

class TestIdentifyRemaining(object):

    def test_fusion(self):
        """ Example where the transcript is shares splice junctions between
            two different genes
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        db = "scratch/toy.db"
        talon.get_counters(db)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(db, build, create_novel_spliced_genes=True)
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")
        init_refs.make_temp_transcript_table(cursor, "toy_build")
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")
        correct_gene_ID = talon.gene_counter.value() + 1

        chrom = "chr1"
        positions = [100, 500, 600, 900, 1010, 5000, 5550, 6000]
        strand = "+"
        edge_IDs = [2, 3, 4]+[ talon.edge_counter.value() + 1, talon.edge_counter.value() + 2 ]
        vertex_IDs = [2, 3, 4, 5, 9, 10, 11]
        v_novelty = [0, 0, 0, 0, 0, 0]

        # Construct temp novel gene db
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")
        fusion = True

        gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info = \
                             talon.process_remaining_mult_cases(chrom, positions,
                                                                strand, edge_IDs,
                                                                vertex_IDs,
                                                                transcript_dict,
                                                                gene_starts, gene_ends,
                                                                edge_dict, location_dict,
                                                                vertex_2_gene, run_info,
                                                                cursor, "temp_gene",
                                                                "temp_transcript",
                                                                fusion)

        assert gene_ID == correct_gene_ID
        assert transcript_dict[frozenset(start_end_info["edge_IDs"])] != None
        assert gene_novelty[0][-2] == "fusion_novel"
        conn.close()
        conn.close()


    def test_intergenic(self):
        """ Example where the transcript is an NIC match to an existing one by
            virtue of a new splice donor.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")
        init_refs.make_temp_transcript_table(cursor, "toy_build")
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")
        correct_gene_ID = talon.gene_counter.value() + 1

        # Construct temp novel gene db
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")

        chrom = "chrX"
        positions = [ 1, 100, 900, 1000]
        edge_IDs = [ talon.edge_counter.value() + 1, talon.edge_counter.value() + 2 ]
        vertex_IDs = [ talon.vertex_counter.value() + 1, talon.vertex_counter.value() + 2 ]
        strand = "+"
        fusion = False

        gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info = \
                             talon.process_remaining_mult_cases(chrom, positions,
                                                                strand, edge_IDs,
                                                                vertex_IDs,
                                                                transcript_dict,
                                                                gene_starts, gene_ends,
                                                                edge_dict, location_dict,
                                                                vertex_2_gene, run_info,
                                                                cursor, "temp_gene",
                                                                "temp_transcript",
                                                                fusion)

        assert gene_ID == correct_gene_ID
        assert transcript_dict[frozenset(start_end_info["edge_IDs"])] != None
        assert gene_novelty[0][-2] == "intergenic_novel"
        conn.close()

    def test_antisense(self):
        """ Example where the transcript is antisense but contains no known
            splice vertices
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)

        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")
        init_refs.make_temp_transcript_table(cursor, "toy_build")
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")
        correct_gene_ID = talon.gene_counter.value() + 1

        # Construct temp novel gene db
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")

        chrom = "chr2"
        positions = [ 1000, 950, 700, 600]
        edge_IDs = [ talon.edge_counter.value() + 1, talon.edge_counter.value() + 2 ]
        vertex_IDs = [ talon.vertex_counter.value() + 1, talon.vertex_counter.value() + 2 ]
        strand = "-"
        fusion = False

        gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info = \
                             talon.process_remaining_mult_cases(chrom, positions,
                                                                strand, edge_IDs,
                                                                vertex_IDs,
                                                                transcript_dict,
                                                                gene_starts, gene_ends,
                                                                edge_dict, location_dict,
                                                                vertex_2_gene, run_info,
                                                                cursor, "temp_gene",
                                                                "temp_transcript",
                                                                fusion)
        assert gene_ID == correct_gene_ID
        assert transcript_dict[frozenset(start_end_info["edge_IDs"])] != None
        assert gene_novelty[0][-2] == "antisense_gene"
        conn.close()

    def test_genomic(self):
        """ Example where the transcript overlaps a gene but contains no known
            splice vertices
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")
        init_refs.make_temp_transcript_table(cursor, "toy_build")
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        # Construct temp novel gene db
        init_refs.make_temp_novel_gene_table(cursor, "toy_build")

        chrom = "chr1"
        positions = [ 1000, 950, 700, 600]
        edge_IDs = [ talon.edge_counter.value() + 1, talon.edge_counter.value() + 2 ]
        vertex_IDs = [ talon.vertex_counter.value() + 1, talon.vertex_counter.value() + 2 ]
        strand = "-"
        fusion = False

        gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info = \
                             talon.process_remaining_mult_cases(chrom, positions,
                                                                strand, edge_IDs,
                                                                vertex_IDs,
                                                                transcript_dict,
                                                                gene_starts, gene_ends,
                                                                edge_dict, location_dict,
                                                                vertex_2_gene, run_info,
                                                                cursor, "temp_gene",
                                                                "temp_transcript",
                                                                fusion)
        correct_gene_ID = fetch_correct_ID("TG3", "gene", cursor)
        assert gene_ID == correct_gene_ID
        assert transcript_dict[frozenset(start_end_info["edge_IDs"])] != None
        assert gene_novelty == []
        assert transcript_novelty[-1][-2] == "genomic_transcript"
        conn.close()
