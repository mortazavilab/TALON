import pytest
import sys
import sqlite3
sys.path.append("..")
import talon as talon
from helper_fns import *
@pytest.mark.integration

class TestIdentifyMonoexonic(object):

    def test_match(self):
        """ Example where the transcript is a moniexonic match.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        talon.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)
        gene_starts, gene_ends = talon.make_gene_start_and_end_dict(cursor)

        chrom = "chr4"
        strand = "-"
        positions = ( 3900, 1100 )

        annotation = talon.identify_monoexon_transcript(chrom, positions, 
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info)

        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        correct_transcript_ID = fetch_correct_ID("TG6-001", "transcript", cursor)
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['start_delta'] == 100
        assert annotation['end_delta'] == -100

        conn.close()

    def test_partial_match(self):
        """ Example where the transcript is short, so it overlaps the 
            annotated transcript but is not an accepted match.
            the start should get assigned to the annotated start, but the end is
            novel """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        talon.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)
        gene_starts, gene_ends = talon.make_gene_start_and_end_dict(cursor)

        chrom = "chr4"
        strand = "-"
        positions = ( 3900, 2900 )

        annotation = talon.identify_monoexon_transcript(chrom, positions,
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info)

        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['start_delta'] == 100
        assert annotation['end_delta'] == None

        conn.close()

    def test_partial_match_3prime(self):
        """ Example where the transcript is short, so it overlaps the
            annotated transcript but is not an accepted match.
            the end should get assigned to the annotated end, but the end is
            novel """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        talon.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)
        gene_starts, gene_ends = talon.make_gene_start_and_end_dict(cursor)

        chrom = "chr4"
        strand = "-"
        positions = ( 2000, 1100 )

        annotation = talon.identify_monoexon_transcript(chrom, positions,
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info)

        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['start_delta'] == None
        assert annotation['end_delta'] == -100

        conn.close()

    def test_overlap_but_no_vertex_match(self):
        """ Example where the transcript is short, so it overlaps the
            annotated transcript but is not an accepted match.
            the start should get assigned to the annotated end, but the end is
            novel """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        talon.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)
        gene_starts, gene_ends = talon.make_gene_start_and_end_dict(cursor)
        tot_vertices = len(vertex_2_gene)
        query = """ SELECT COUNT(*) FROM temp_monoexon """
        tot_monoexonic = cursor.execute(query).fetchone()[0]

        chrom = "chr4"
        strand = "-"
        positions = ( 2500, 2000 )

        annotation = talon.identify_monoexon_transcript(chrom, positions,
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info)

        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        print(annotation['start_vertex'])
        print(annotation['end_vertex'])
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['start_delta'] == None
        assert annotation['end_delta'] == None

        # Now check if the transcript got added to the right data structures
        assert len(vertex_2_gene) == tot_vertices + 2
        assert cursor.execute(query).fetchone()[0] == tot_monoexonic + 1

        conn.close()

    def test_antisense(self):
        """ Example where the transcript is antisense """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        talon.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)
        gene_starts, gene_ends = talon.make_gene_start_and_end_dict(cursor)

        chrom = "chr4"
        strand = "+"
        positions = ( 1300, 3900 )

        annotation = talon.identify_monoexon_transcript(chrom, positions,
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info)

        anti_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        gene_novelty_types = [ x[-2] for x in annotation['gene_novelty']]
        t_novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_novelty'][0][-1] == "TRUE"
        assert "antisense_gene" in gene_novelty_types
        assert "antisense_transcript" in t_novelty_types

        conn.close() 
