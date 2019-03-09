import pytest
import sys
import sqlite3
sys.path.append("..")
import talon as talon
from helper_fns import *
@pytest.mark.unit

class TestIdentifyGeneOnVertexBasis(object):

    def test_perfect_match(self):
        """ Example where the vertices perfectly match a gene.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, "toy_build")
        run_info = talon.init_run_info(cursor, build)
        vertex2gene = talon.make_vertex_2_gene_dict(cursor)

        vertex_IDs = (1, 2, 3, 4, 5, 6)
        strand = "+"

        gene_ID = talon.find_gene_match_on_vertex_basis(vertex_IDs, strand, vertex2gene)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        assert gene_ID == correct_gene_ID
        conn.close()

    def test_NNC_type_match(self):
        """ Example where some vertices match a gene, while others don't.
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, "toy_build")
        run_info = talon.init_run_info(cursor, build)
        vertex2gene = talon.make_vertex_2_gene_dict(cursor)

        vertex_IDs = (1, 200, 3, 4, 5, 6)
        strand = "+"

        gene_ID = talon.find_gene_match_on_vertex_basis(vertex_IDs, strand, vertex2gene)

        correct_gene_ID = fetch_correct_ID("TG1", "gene", cursor)
        assert gene_ID == correct_gene_ID
        conn.close()

    def test_no_match(self):
        """ Example where no match exists """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        talon.make_temp_novel_gene_table(cursor, "toy_build")
        run_info = talon.init_run_info(cursor, build)
        vertex2gene = talon.make_vertex_2_gene_dict(cursor)

        vertex_IDs = (1000, 2000, 3000, 4000)
        strand = "+"

        gene_ID = talon.find_gene_match_on_vertex_basis(vertex_IDs, strand, vertex2gene)

        assert gene_ID == None
        conn.close()
