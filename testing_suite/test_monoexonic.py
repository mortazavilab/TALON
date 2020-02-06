import pytest
from talon import talon, init_refs
from .helper_fns import fetch_correct_ID, get_db_cursor
@pytest.mark.integration

class TestIdentifyMonoexonic(object):

    def test_match(self):
        """ Example where the transcript is a monoexonic match.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        init_refs.make_temp_novel_gene_table(cursor, build)
        init_refs.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr4"
        strand = "-"
        positions = ( 3900, 1100 )

        annotation = talon.identify_monoexon_transcript(chrom, positions, 
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info,
                                               'temp_gene', 'temp_monoexon')

        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        correct_transcript_ID = fetch_correct_ID("TG6-001", "transcript", cursor)
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['start_delta'] == 100
        assert annotation['end_delta'] == -100

        conn.close()

    def test_partial_match(self):
        """ Example where the transcript overlaps a single-exon transcript,
            but is shorter. In the past, the start would be assigned to the 
            annotated start, and the end would be novel. This is no longer
            the case- at this time, the transcript will be assigned to
            the annotated match. """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        init_refs.make_temp_novel_gene_table(cursor, build)
        init_refs.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr4"
        strand = "-"
        positions = ( 3900, 2900 )

        annotation = talon.identify_monoexon_transcript(chrom, positions,
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info,
                                               'temp_gene', 'temp_monoexon')

        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        correct_transcript_ID = fetch_correct_ID("TG6-001", "transcript", cursor)
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['transcript_ID'] == correct_transcript_ID
        assert annotation['start_delta'] == 100
        assert annotation['end_delta'] == -1900

        conn.close()

# Commenting out these tests for now because they are redundant. But saving in 
# case they might be useful down the line.

#    def test_partial_match_3prime(self):
#        """ Example where the transcript is short, so it overlaps the
#            annotated transcript but is not an accepted match.
#            the end should get assigned to the annotated end, but the end is
#            novel """
#
#        conn, cursor = get_db_cursor()
#        build = "toy_build"
#        database = "scratch/toy.db"
#        talon.get_counters(database)
#        init_refs.make_temp_novel_gene_table(cursor, build)
#        init_refs.make_temp_monoexonic_transcript_table(cursor, build)
#        edge_dict = init_refs.make_edge_dict(cursor)
#        location_dict = init_refs.make_location_dict(build, cursor)
#        run_info = talon.init_run_info(database, build)
#        transcript_dict = init_refs.make_transcript_dict(cursor, build)
#        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
#        gene_starts, gene_ends = init_refs.make_gene_start_and_end_dict(cursor, build)
#
#        chrom = "chr4"
#        strand = "-"
#        positions = ( 2000, 1100 )
#
#        annotation = talon.identify_monoexon_transcript(chrom, positions,
#                                               strand, cursor,
#                                               location_dict, edge_dict,
#                                               transcript_dict, vertex_2_gene,
#                                               gene_starts, gene_ends, run_info,
#                                               'temp_gene', 'temp_monoexon')
#
#        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
#        assert annotation['gene_ID'] == correct_gene_ID
#        assert annotation['start_delta'] == None
#        assert annotation['end_delta'] == -100
#
#        conn.close()
#
#    def test_overlap_but_no_vertex_match(self):
#        """ Example where the transcript is short, so it overlaps the
#            annotated transcript but is not an accepted match.
#            the start should get assigned to the annotated end, but the end is
#            novel """
#
#        conn, cursor = get_db_cursor()
#        build = "toy_build"
#        database = "scratch/toy.db"
#        talon.get_counters(database)
#        init_refs.make_temp_novel_gene_table(cursor, build)
#        init_refs.make_temp_monoexonic_transcript_table(cursor, build)
#        edge_dict = init_refs.make_edge_dict(cursor)
#        location_dict = init_refs.make_location_dict(build, cursor)
#        run_info = talon.init_run_info(database, build)
#        transcript_dict = init_refs.make_transcript_dict(cursor, build)
#        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
#        gene_starts, gene_ends = init_refs.make_gene_start_and_end_dict(cursor, build)
#        tot_vertices = len(vertex_2_gene)
#        query = """ SELECT COUNT(*) FROM temp_monoexon """
#        tot_monoexonic = cursor.execute(query).fetchone()[0]
#
#        chrom = "chr4"
#        strand = "-"
#        positions = ( 2500, 2000 )
#
#        annotation = talon.identify_monoexon_transcript(chrom, positions,
#                                               strand, cursor,
#                                               location_dict, edge_dict,
#                                               transcript_dict, vertex_2_gene,
#                                               gene_starts, gene_ends, run_info,
#                                               'temp_gene', 'temp_monoexon')
#
#        correct_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
#        print(annotation['start_vertex'])
#        print(annotation['end_vertex'])
#        assert annotation['gene_ID'] == correct_gene_ID
#        assert annotation['start_delta'] == None
#        assert annotation['end_delta'] == None
#
#        # Now check if the transcript got added to the right data structures
#        assert len(vertex_2_gene) == tot_vertices + 2
#        assert cursor.execute(query).fetchone()[0] == tot_monoexonic + 1
#
#        conn.close()
#
    def test_antisense(self):
        """ Example where the transcript is antisense """

        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        init_refs.make_temp_novel_gene_table(cursor, build)
        init_refs.make_temp_monoexonic_transcript_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr4"
        strand = "+"
        positions = ( 1300, 3900 )

        annotation = talon.identify_monoexon_transcript(chrom, positions,
                                               strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info,
                                               'temp_gene', 'temp_monoexon')

        anti_gene_ID = fetch_correct_ID("TG6", "gene", cursor)
        gene_novelty_types = [ x[-2] for x in annotation['gene_novelty']]
        t_novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_novelty'][0][-1] == "TRUE"
        assert "antisense_gene" in gene_novelty_types
        assert "antisense_transcript" in t_novelty_types

        conn.close() 
