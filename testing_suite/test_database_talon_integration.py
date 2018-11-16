# Compared to a lot of the other tests in the suite, this one is intended
# to make sure that the major parts of TALON (database init, dataset
# addition over time) are working correctly.

import pytest
import subprocess
@pytest.mark.integration

class TestDatabaseTalonIntegration(object):
    def test_db_initialization(self):
        """ Initializes a TALON database from a very minimal GTF file (KRT17
            example). Outfiles are written to the scratch area of the testing
            suite. 
        """
        try:
            subprocess.check_output(
                     ["python", "../initialize_talon_database.py", 
                         "--f", "input_files/KRT17_test_case/KRT17.gtf",  
                         "--a",  "KRT17_test",  
                         "--g",  "hg38", "--o", "scratch/KRT17"])
        except Exception as e:
            pytest.fail("Database initialization failed on KRT17 test")

    @pytest.mark.incremental
    def test_TALON_simple_run(self):
        """ Once the database has been initialized, try running TALON on a 
            single SAM transcript (that we know matches KRT17-001). First,
            establish that TALON runs without crashing."""

        try:
             subprocess.check_output(
                     ["python", "../talon.py",
                      "--f", "input_files/KRT17_test_case/trial1_config.csv",
                      "-a", "scratch/KRT17.db",
                      "-b", "hg38",
                      "--o", "scratch/trial1"])
        except:
            pytest.fail("TALON trial 1 failed on KRT17 test")

    @pytest.mark.incremental
    def test_TALON_simple_run_correctness(self):
        """ Once we've established that TALON ran on the KRT17 example without
            a crash, we need to check that the correct identity was assigned
            to the transcript. """

        infile = "scratch/trial1_talon.tsv"

        line_num = 0
        gene_id_col_index = None
        transcript_id_col_index = None
        with open(infile, 'r') as f:
            for line in f:
                line = line.strip().split("\t")

                # From the header, figure out which columns contain the assigned 
                # gene and transcript IDs, as well as 3' and 5' end differences
                if line_num == 0:
                    gene_id_col_index = line.index("gene_id")
                    transcript_id_col_index = line.index("transcript_id")

                # Check that the transcript got assigned to gene 1 and 
                # transcript 1
                else:
                    assert line[gene_id_col_index] == "1"
                    assert line[transcript_id_col_index] == "1"

                line_num += 1
                
