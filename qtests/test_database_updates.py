import pytest
import sys
import sqlite3
sys.path.append("..")
import talonQ as talon
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

        # Test if gene with ID 5 is there, but make sure we didn't add 
        # duplicates of the other genes
        query = "SELECT * FROM genes"
        gene_IDs = [ x['gene_ID'] for x in cursor.execute(query)]
        assert 5 in gene_IDs
        assert len(gene_IDs) == 5
        conn.close()

    def test_transcript_update(self):
        """ Try to add novel transcript entries to database while ignoring 
            duplicates
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        edge_dict = talon.make_edge_dict(cursor)
        run_info = talon.init_run_info(cursor, build)
        transcript_dict = talon.make_transcript_dict(cursor)
        talon.create_transcript( 1, (1,2,3), (1,2,3,4), transcript_dict, run_info)

        batch_size = 5
        talon.batch_add_transcripts(cursor, transcript_dict, batch_size)

        # Test if transcript with ID 5 is there, but make sure we didn't add
        # duplicates of the others
        query = "SELECT * FROM transcripts"
        cursor.execute(query)
        transcript_IDs = [ x['transcript_ID'] for x in cursor.fetchall()]
        assert 5 in transcript_IDs
        assert len(transcript_IDs) == 5
        conn.close()


