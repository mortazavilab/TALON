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
                                      "annot_transcript_id",	                                             "annot_gene_name", 
                                      "annot_transcript_name", 
                                      "n_exons", 
                                      "length", "gene_novelty", 
                                      "transcript_novelty", 
                                      "ISM_subtype", 
                                      "PB65_B017", "PB65_B018", "D12"]
        assert data.shape[0] == 15
