import pytest
import sqlite3
from talon import query_utils as qutils
@pytest.mark.dbunit

class TestQueries(object):
    """ Make sure that the queries in Query Utils are working correctly """

    def test_count_observed_reads(self):
        """ Count the number of observed reads for the provided dataset. """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Check each dataset in chr11_and_Tcf3 example set
        dataset = "PB65_B017"
        assert qutils.count_observed_reads(cursor, dataset) == 8
        dataset = "PB65_B018"
        assert qutils.count_observed_reads(cursor, dataset) == 3
        dataset = "D12"
        assert qutils.count_observed_reads(cursor, dataset) == 5

        conn.close()
 
    def test_count_known_genes_detected(self):
        """ Count the number of known genes detected in each dataset """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Check each dataset in chr11_and_Tcf3 example set
        dataset = "PB65_B017"
        known_genes = qutils.fetch_all_known_genes_detected(cursor, dataset)
        assert sorted(known_genes) == [ 5, 723, 2987 ]
        assert qutils.count_known_genes_detected(cursor, dataset) == 3

        dataset = "PB65_B018"
        assert qutils.count_known_genes_detected(cursor, dataset) == 3

        dataset = "D12"
        assert qutils.count_known_genes_detected(cursor, dataset) == 1

        conn.close()

    def test_count_novel_genes_detected(self):
        """ Count the number of novel genes detected per dataset """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Check each dataset in chr11_and_Tcf3 example set
        dataset = "PB65_B017"
        novel_genes = qutils.fetch_all_novel_genes_detected(cursor, dataset)
        assert len(novel_genes) == 2

        dataset = "PB65_B018"
        novel_genes = qutils.fetch_all_novel_genes_detected(cursor, dataset)
        assert len(novel_genes) == 0

        dataset = "D12"
        novel_genes = qutils.fetch_all_novel_genes_detected(cursor, dataset)
        assert len(novel_genes) == 1

    def test_count_known_transcripts_detected(self):
        """ Count the number of known transcripts detected in each dataset """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Check each dataset in chr11_and_Tcf3 example set
        dataset = "PB65_B017"
        known_transcripts = qutils.fetch_all_known_transcripts_detected(cursor, dataset)
        assert len(known_transcripts) == 1

        dataset = "PB65_B018"
        known_transcripts = qutils.fetch_all_known_transcripts_detected(cursor, dataset)
        assert len(known_transcripts) == 1

        dataset = "D12"
        known_transcripts = qutils.fetch_all_known_transcripts_detected(cursor, dataset)
        assert len(known_transcripts) == 1

        conn.close() 

    def test_count_novel_transcripts_detected(self):
        """ Count the number of novel transcripts detected in each dataset """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Check each dataset in chr11_and_Tcf3 example set
        dataset = "PB65_B017"
        novel_transcripts = qutils.fetch_novel_transcripts(cursor, dataset)
        assert len(novel_transcripts) == 7

        dataset = "PB65_B018"
        novel_transcripts = qutils.fetch_novel_transcripts(cursor, dataset)
        assert len(novel_transcripts) == 2

        dataset = "D12"
        novel_transcripts = qutils.fetch_novel_transcripts(cursor, dataset)
        assert len(novel_transcripts) == 4

        conn.close()
       
    def test_count_antisense_genes(self):
        """ Count the number of antisense genes in the datasets """
 
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Check each dataset in chr11_and_Tcf3 example set
        dataset = "PB65_B017"
        genes = qutils.fetch_antisense_genes(cursor, dataset)
        assert len(genes) == 1

        dataset = "PB65_B018"
        genes = qutils.fetch_antisense_genes(cursor, dataset)
        assert len(genes) == 0

        dataset = "D12"
        genes = qutils.fetch_antisense_genes(cursor, dataset)
        assert len(genes) == 1

        # Now try all 3 datasets at once
        datasets = ["PB65_B017", "PB65_B018", "D12"]
        genes = qutils.fetch_antisense_genes(cursor, datasets)
        assert len(genes) == 2
        conn.close()

    def test_count_intergenic_genes(self):
        """ Count intergenic novel genes """               
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        genes = qutils.fetch_intergenic_novel_genes(cursor, dataset)
        assert len(genes) == 1
        conn.close()
 
    def test_count_all_ISM_transcripts(self):
        """ """
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        transcripts = qutils.fetch_all_ISM_transcripts(cursor, dataset)
        assert len(transcripts) == 2

        dataset = "D12"
        transcripts = qutils.fetch_all_ISM_transcripts(cursor, dataset)
        assert len(transcripts) == 0
        conn.close()

    def test_count_prefix_ISM_transcripts(self):
        """ """
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B018"
        transcripts = qutils.fetch_prefix_ISM_transcripts(cursor, dataset)
        assert len(transcripts) == 0
        conn.close()

    def test_count_suffix_ISM_transcripts(self):
        """ """
        # TODO: resolve why this test is failing. I believe it has to do with 
        # parallel transcript order switches
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        transcripts = qutils.fetch_suffix_ISM_transcripts(cursor, dataset)
        assert len(transcripts) == 1
        conn.close()

    def test_count_NIC_transcripts(self):
        """ """
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B018"
        transcripts = qutils.fetch_NIC_transcripts(cursor, dataset)
        assert len(transcripts) == 2

        # Now try all 3 datasets at once
        datasets = ["PB65_B017", "PB65_B018", "D12"]
        transcripts = qutils.fetch_NIC_transcripts(cursor, datasets)
        assert len(transcripts) == 3
        conn.close()

    def test_count_NNC_transcripts(self):
        """ """
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        transcripts = qutils.fetch_NNC_transcripts(cursor, dataset)
        assert len(transcripts) == 0

        dataset = "D12"
        transcripts = qutils.fetch_NNC_transcripts(cursor, dataset)
        assert len(transcripts) == 2
        conn.close()

    def test_count_antisense_transcripts(self):
        """ """
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        transcripts = qutils.fetch_antisense_transcripts(cursor, dataset)
        assert len(transcripts) == 1

        dataset = "D12"
        transcripts = qutils.fetch_antisense_transcripts(cursor, dataset)
        assert len(transcripts) == 1
        conn.close()

