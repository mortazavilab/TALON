import pytest
import subprocess
import pandas as pd
#from talon import talon_abudance as abd

@pytest.mark.integration
class TestAbundance(object):
    """ Make sure that the abundance utility is working correctly """

    def test_base_settings(self):
        """ Test abundance utility without a datasets or whitelist file. """
        database =  "scratch/chr11_and_Tcf3.db"
        try:
            subprocess.check_output(
                ["talon_abundance", "--db", database, 
                 "-a", "gencode_vM7",
                 "-b", "mm10",
                 "--o", "scratch/chr11_and_Tcf3_base"])
        except:
            pytest.fail("Talon abundance crashed on basic case") 

        # Now check the correctness of the abundance file
        abd = "scratch/chr11_and_Tcf3_base_talon_abundance.tsv"
        data = pd.read_csv(abd, sep="\t", header = 0)
        print(data)
        assert list(data.columns) == ["gene_ID", "transcript_ID", 
                                      "annot_gene_id", 
                                      "annot_transcript_id", "annot_gene_name", 
                                      "annot_transcript_name", "n_exons", 
                                      "length", "gene_novelty", 
                                      "transcript_novelty", 
                                      "ISM_subtype", 
                                      "PB65_B017", "PB65_B018", "D12"]
        assert data.shape[0] == 15

    def test_with_whitelist(self):
        """ Test abundance utility with a transcript whitelist """
        database =  "scratch/chr11_and_Tcf3.db"
        whitelist = "input_files/chr11_and_Tcf3/testing_whitelist.txt"
        try:
            subprocess.check_output(
                ["talon_abundance", "--db", database,
                 "-a", "gencode_vM7",
                 "-b", "mm10",
                 "--whitelist", whitelist,
                 "--o", "scratch/chr11_and_Tcf3_whitelist"])
        except:
            pytest.fail("Talon abundance crashed on whitelist case")
        
        # Now check the correctness of the abundance file
        abd = "scratch/chr11_and_Tcf3_whitelist_talon_abundance_filtered.tsv"
        data = pd.read_csv(abd, sep="\t", header = 0)
        
        print(data)
        assert list(data.columns) == ["gene_ID", "transcript_ID",
                                      "annot_gene_id",
                                      "annot_transcript_id", "annot_gene_name",
                                      "annot_transcript_name", "n_exons",
                                      "length", "gene_novelty",
                                      "transcript_novelty",
                                      "ISM_subtype",
                                      "PB65_B017", "PB65_B018", "D12"]
        assert set(data.transcript_ID) == set([28,1744])
        assert int(data.loc[data['transcript_ID'] == 28]['PB65_B017']) == 1
        assert int(data.loc[data['transcript_ID'] == 28]['PB65_B018']) == 0
        assert int(data.loc[data['transcript_ID'] == 28]['D12']) == 0

    def test_with_dataset_list(self):
        """ Test abundance utility with a transcript whitelist """
        database =  "scratch/chr11_and_Tcf3.db"
        datasets = "input_files/chr11_and_Tcf3/testing_datasets.txt"
        try:
            subprocess.check_output(
                ["talon_abundance", "--db", database,
                 "-a", "gencode_vM7",
                 "-b", "mm10",
                 "--datasets", datasets,
                 "--o", "scratch/chr11_and_Tcf3_dset"])
        except:
            pytest.fail("Talon abundance crashed on whitelist case")

        # Now check the correctness of the abundance file
        abd = "scratch/chr11_and_Tcf3_dset_talon_abundance.tsv"
        data = pd.read_csv(abd, sep="\t", header = 0)

        print(data)
        assert set(list(data.columns)) == set(["gene_ID", "transcript_ID",
                                               "annot_gene_id",
                                               "annot_transcript_id", "annot_gene_name",
                                               "annot_transcript_name", "n_exons",
                                               "length", "gene_novelty",
                                               "transcript_novelty",
                                               "ISM_subtype",
                                               "PB65_B018", "D12"])
        assert data.shape[0] == 8
        assert set(data.transcript_ID) == set([1744, 8437, 8453, 8456, 8457, 8458, 8459, 8460]) 

