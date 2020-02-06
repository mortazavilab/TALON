import pytest
import sqlite3
from talon import talon, init_refs
from .helper_fns import fetch_correct_ID, get_db_cursor
@pytest.mark.integration

class TestIdentifyFSM(object):

    def test_FSM_perfect(self):
        """ Example where the transcript is a perfect full splice match.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)        
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr1"
        strand = "+"
        positions = [ 1, 100, 500, 600, 900, 1000 ]


        annotation = talon.identify_transcript(chrom, positions, strand, cursor, 
                                               location_dict, edge_dict, 
                                               transcript_dict, vertex_2_gene, 
                                               gene_starts, gene_ends, run_info,
                                               "temp_gene")

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
        database = "scratch/toy.db"
        talon.get_counters(database)

        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr2"
        strand = "+"
        positions = [ 1, 100, 500, 600, 900, 1500 ]


        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info,
                                               "temp_gene")

        correct_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert annotation['end_delta'] == None
        conn.close()

    def test_NIC_instead_of_ISM(self):
        """ Test case where the transcript looks like an ISM, but has known 
            starts and ends. In past TALON versions, this was considered NIC,
            but expected behavior is now ISM """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)

        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr3"
        strand = "+"
        positions = ( 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200 )

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info,
                                               "tmp_gene") 

        correct_gene_ID = fetch_correct_ID("TG5", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert "ISM_transcript" in novelty_types
        conn.close()

    def test_ISM_suffix(self):
        """ Example where the transcript is a suffix ISM
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)

        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr1"
        strand = "+"
        positions = [ 550, 600, 900, 1200 ]

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info,
                                               "temp_gene")

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert "ISM_transcript" in novelty_types 
        assert "ISM-suffix_transcript" in novelty_types
        assert annotation['start_delta'] == 50
        conn.close()

    def test_ISM_prefix(self):
        """ Example where the transcript is a prefix ISM
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)

        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr1"
        strand = "+"
        positions = ( 1, 100, 500, 620 )

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info, "temp_gene")

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert "ISM-prefix_transcript" in novelty_types
        assert annotation['start_delta'] == 0
        assert annotation['end_delta'] == 20
        conn.close()

    def test_ISM_internal(self):
        """ Example where the transcript matches an internal exon
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

        chrom = "chr1"
        strand = "+"
        positions = ( 500, 600 )

        annotation = talon.identify_monoexon_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info, "temp_gene", "temp_monoexon")

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert "ISM_transcript" in novelty_types 
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_NIC(self):
        """ Example where the transcript skips an exon
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)

        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr1"
        strand = "+"
        positions = ( 1, 100, 900, 1000 ) 

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info, "temp_gene")

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert "NIC_transcript" in novelty_types
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_NNC(self):
        """ Example where the transcript skips an exon and has a novel splice
            donor
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database) 

        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr1"
        strand = "+"
        positions = [ 1, 50, 900, 1000 ]

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info, "temp_gene")

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert "NNC_transcript" in novelty_types
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_spliced_antisense(self):
        """ Example where the transcript matches known vertices but is antisense
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)

        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr2"
        strand = "-"
        positions = [ 1000, 900, 600, 500, 100, 1 ]

        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info, "temp_gene")

        anti_gene_ID = fetch_correct_ID("TG2", "gene", cursor)
        gene_novelty_types = [ x[-2] for x in annotation['gene_novelty']]
        t_novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_novelty'][0][-1] == "TRUE"
        assert "antisense_gene" in gene_novelty_types
        assert "antisense_transcript" in t_novelty_types
        assert annotation['start_delta'] == annotation['end_delta'] == 0
        conn.close()

    def test_genomic_unspliced(self):
        """ Monoexonic fragment that overlaps gene 1 """
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

        chrom = "chr1"
        strand = "+"
        positions = (1, 990)

        annotation = talon.identify_monoexon_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info, "temp_gene", "temp_monoexon")

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert annotation['gene_ID'] == correct_gene_ID
        assert "genomic_transcript" in novelty_types
        assert annotation['end_delta'] == -10
        conn.close()        

    def test_NIC_with_all_known_edges(self):
        """ Test case derived from a real mouse Map2k4 read. All of edges are
            known (except 3'), yet the read is NIC not FSM/ISM """

        database = "scratch/Map2k4.db"
        talon.get_counters(database)
        conn = sqlite3.connect(database)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()        
        build = "mm10"
        init_refs.make_temp_novel_gene_table(cursor, build)
        edge_dict = init_refs.make_edge_dict(cursor)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)
        gene_starts = init_refs.make_gene_start_or_end_dict(cursor, build, "start")
        gene_ends = init_refs.make_gene_start_or_end_dict(cursor, build, "end")

        chrom = "chr11"
        strand = "-"
        positions = [65788254, 65788136, 65775765, 65775733, 65756371, 65756269, 65735366, 65735192, 65719603, 65719484, 65712297, 65712178, 65709983, 65709932, 65707111, 65706984, 65696365, 65696288, 65693570, 65693422, 65691773, 65691728, 65690804, 65689322]
      
        annotation = talon.identify_transcript(chrom, positions, strand, cursor,
                                               location_dict, edge_dict,
                                               transcript_dict, vertex_2_gene,
                                               gene_starts, gene_ends, run_info, "temp_gene")

        assert annotation['gene_ID'] == 1
        assert annotation['transcript_ID'] == 8
        novelty_types = [ x[-2] for x in annotation['transcript_novelty']]
        assert "NIC_transcript" in novelty_types

        conn.close()

