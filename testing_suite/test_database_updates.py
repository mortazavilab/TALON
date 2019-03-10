import pytest
import sys
import sqlite3
sys.path.append("..")
import talon as talon
import dstruct
from helper_fns import *
@pytest.mark.dbunit

class TestDatabaseUpdates(object):

    def test_abundance(self):
        """ Try to add abundance entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)

        abundance = [ ( 1, "test", 5),
                      ( 2, "test", 1),
                      ( 3, "test", 2)]
        batch_size = 2
        talon.batch_add_abundance(cursor, abundance, batch_size)

        # Test if items are there
        query = "SELECT * FROM abundance"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 3
        conn.close()

    def test_datasets(self):
        """ Try to add dataset metadata to database """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)

        datasets = [ ( 1, "toy", "toy", "toy") ]
        talon.add_datasets(cursor, datasets)

        # Test if items are there
        query = "SELECT * FROM dataset"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 1
        conn.close()

    def test_observed(self):
        """ Try to add observed entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)

        observed = [ ( 1, 1, 1, "read1", "test", 1, 2, 0, 0, 100),
                     ( 2, 1, 1, "read2", "test", 1, 2, 0, 0, 100),
                     ( 3, 1, 1, "read3", "test", 1, 2, 0, 0, 100) ]
        batch_size = 1
        talon.batch_add_observed(cursor, observed, batch_size)

        # Test if items are there
        query = "SELECT * FROM observed"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 3
        conn.close()

    def test_gene_annot(self):
        """ Try to add gene annotation entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)

        annot = [ ( 1, "toy", "TALON", "status", "NOVEL"),
                  ( 2, "toy", "TALON", "status", "NOVEL") ]
        batch_size = 1
        talon.batch_add_annotations(cursor, annot, "gene", batch_size)

        # Test if items are there
        query = "SELECT * FROM gene_annotations WHERE value = 'NOVEL'"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 2
        conn.close()

    def test_transcript_annot(self):
        """ Try to add transcript annotation entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)

        annot = [ ( 1, "toy", "TALON", "status", "NOVEL"),
                  ( 2, "toy", "TALON", "status", "NOVEL") ]
        batch_size = 2
        talon.batch_add_annotations(cursor, annot, "transcript", batch_size)

        # Test if items are there
        query = "SELECT * FROM transcript_annotations WHERE value = 'NOVEL'"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 2
        conn.close()

    def test_exon_annot(self):
        """ Try to add exon annotation entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)

        annot = [ ( 1, "toy", "TALON", "status", "NOVEL"),
                  ( 2, "toy", "TALON", "status", "NOVEL") ]
        batch_size = 3
        talon.batch_add_annotations(cursor, annot, "exon", batch_size)

        # Test if items are there
        query = "SELECT * FROM exon_annotations WHERE value = 'NOVEL'"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 2
        conn.close()

    def test_gene_update(self):
        """ Try to add novel gene entries to database while ignoring duplicates
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)
        talon.make_temp_novel_gene_table(cursor, build)
        talon.create_gene("chr4", 1, 1000, "+", cursor, run_info)

        talon.add_genes(cursor)

        # Test if gene with ID 6 is there, but make sure we didn't add 
        # duplicates of the other genes
        query = "SELECT * FROM genes"
        gene_IDs = [ x['gene_ID'] for x in cursor.execute(query)]
        assert 7 in gene_IDs
        assert len(gene_IDs) == 7
        conn.close()

    def test_transcript_update(self):
        """ Try to add novel transcript entries to database while ignoring 
            duplicates
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor, build)
        talon.create_transcript( 1, (1,2,3), (1,2,3,4), transcript_dict, run_info)

        batch_size = 5
        talon.batch_add_transcripts(cursor, transcript_dict, batch_size)

        # Test if transcript with ID 7 is there, but make sure we didn't add
        # duplicates of the others
        query = "SELECT * FROM transcripts"
        cursor.execute(query)
        transcript_IDs = [ x['transcript_ID'] for x in cursor.fetchall()]
        assert 8 in transcript_IDs
        assert len(transcript_IDs) == 8
        conn.close()

    def test_edge_update(self):
        """ Try to add novel exons and introns. """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)
        orig_n_edges = run_info.edge

        talon.create_edge(2, 1, "exon", "-", edge_dict, run_info)

        batch_size = 10
        talon.batch_add_edges(cursor, edge_dict, batch_size)
        
        # Test if the edge table has the correct number of edges now
        query = "SELECT * FROM edge"
        cursor.execute(query)
        edge_IDs = [ x['edge_ID'] for x in cursor.fetchall()]
        assert orig_n_edges + 1 in edge_IDs
        assert len(edge_IDs) == orig_n_edges + 1
        conn.close()

    def test_location_update(self):
        """ Update locations """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        location_dict = talon.make_location_dict(build, cursor)
        run_info = talon.init_run_info(cursor, build)
        orig_n_pos = run_info.vertex

        talon.create_vertex("chr4", 2000, run_info, location_dict)

        batch_size = 10
        talon.batch_add_locations(cursor, location_dict, batch_size)

        # Test if the table has the correct number of locations now
        query = "SELECT * FROM location"
        cursor.execute(query)
        loc_IDs = [ x['location_ID'] for x in cursor.fetchall()]
        assert orig_n_pos + 1 in loc_IDs
        assert len(loc_IDs) == orig_n_pos + 1
        conn.close()

    def test_vertex2gene_update(self):
        """ Update vertex to gene relationships """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        vertex_2_gene = talon.make_vertex_2_gene_dict(cursor)
        # Pretend that vertex 1 and 2 can now belong to gene 2 as well as 1
        talon.update_vertex_2_gene(2, (1,2), "-", vertex_2_gene) 
        # Add redundant assignments
        talon.update_vertex_2_gene(1, (1,2,3,4,5,6), "+", vertex_2_gene)

        batch_size = 100
        talon.batch_add_vertex2gene(cursor, vertex_2_gene, batch_size)

        # Use queries to check if the insert worked as expected
        query = "SELECT * FROM vertex WHERE vertex_ID = '1'"
        cursor.execute(query)
        gene_IDs = [ x['gene_ID'] for x in cursor.fetchall()]
        assert gene_IDs == [1, 2]

        query = "SELECT * FROM vertex WHERE gene_ID = '1'"
        cursor.execute(query)
        vertex_IDs = [ x['vertex_ID'] for x in cursor.fetchall()]
        assert vertex_IDs == [1, 2, 3, 4, 5, 6]

    def test_counter_update(self):
        """ Update counters """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        run_info = talon.init_run_info(cursor, build)

        # Change the counter values to some arbitrary numbers
        run_info.genes = 10
        run_info.transcripts = 20
        run_info.edge = 2000
        run_info.vertex = 10000
        run_info.dataset = 30
        run_info.observed = 400

        # Now try the update
        talon.update_counter(cursor, run_info)
        run_info = None

        # Check results with queries
        run_info_2 = talon.init_run_info(cursor, build)
        assert run_info_2.genes == 10
        assert run_info_2.transcripts == 20
        assert run_info_2.edge == 2000
        assert run_info_2.vertex == 10000
        assert run_info_2.dataset == 30
        assert run_info_2.observed == 400

