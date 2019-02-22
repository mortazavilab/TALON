import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
from helper_fns import *
@pytest.mark.integration

class TestIdentifyFSM(object):

    def test_FSM_perfect(self):
        """ Example where the transcript is a perfect full splice match.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)        
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr1"
        strand = "+"
        positions = ( 1, 100, 500, 600, 900, 1000 )


        annotation = talon.identify_transcript(chrom, positions, strand, cursor, 
                                               location_dict, edge_dict, 
                                               transcript_dict, vertex_2_gene, 
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor) 
        correct_transcript_ID = fetch_correct_ID("TG1-001", "transcript", cursor)
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['transcript_ID'] == correct_transcript_ID
        assert annotation['transcript_novelty'] == []
        conn.close()

    def test_FSM_end_diff(self):
        """ Example where the transcript is an FSM but has a difference on
            the ends large enough to be novel.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr2"
        strand = "+"
        positions = ( 1, 100, 500, 600, 900, 1500 )


        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert novelty_types == ["FSM_transcript", "related_transcript_IDs",
                                 "3p_novel"]
        assert annotation['end_delta'] == None
        conn.close()

    def test_ISM_suffix(self):
        """ Example where the transcript is a suffix ISM
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr1"
        strand = "+"
        positions = ( 550, 600, 900, 1500 )

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert novelty_types == ["ISM_transcript", "ISM-suffix_transcript",
                                 "related_transcript_IDs"]
        assert annotation['start_delta'] == -50
        conn.close()

    def test_ISM_prefix(self):
        """ Example where the transcript is a prefix ISM
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr1"
        strand = "+"
        positions = ( 1, 100, 500, 600 )

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert novelty_types == ["ISM_transcript", "ISM-prefix_transcript",
                                 "related_transcript_IDs"]
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_ISM_internal(self):
        """ Example where the transcript matches an internal exon
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr1"
        strand = "+"
        positions = ( 500, 600 )

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert novelty_types == ["ISM_transcript", "related_transcript_IDs"]
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_NIC(self):
        """ Example where the transcript skips an exon
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr1"
        strand = "+"
        positions = ( 1, 100, 900, 1000 ) 

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert novelty_types == ["NIC_transcript"]
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_NNC(self):
        """ Example where the transcript skips an exon and has a novel splice
            donor
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr1"
        strand = "+"
        positions = ( 1, 50, 900, 1000 )

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert novelty_types == ["NNC_transcript"]
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_spliced_antisense(self):
        """ Example where the transcript matches known vertices but is antisense
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr2"
        strand = "-"
        positions = ( 1000, 900, 600, 500, 100, 1 )

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        anti_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        gene_novelty_types = [ x[-2] for x in annotation['gene_novelty']]
        t_novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_novelty'][0][-1] == "TRUE"
        assert gene_novelty_types == ["antisense_gene", "related_gene_IDs"]
        assert t_novelty_types == ["antisense_transcript"]
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_genomic_unspliced(self):
        """ Monoexonic fragment that overlaps gene 1 """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, build)
        edge_dict = talon.make_edge_dict(cursor)
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)

        chrom = "chr1"
        strand = "+"
        positions = (1, 990)

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               run_info)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert novelty_types == ["genomic_transcript"]
        assert annotation['end_delta'] == 10
        conn.close()        

