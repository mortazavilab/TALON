# Compared to a lot of the other tests in the suite, this one is intended
# to make sure that the major parts of TALON (database init, dataset
# addition over time) are working correctly.

import pytest
import subprocess
import sqlite3
import sys
sys.path.append("..")
import talon as TALON
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
                # gene and transcript IDs
                if line_num == 0:
                    gene_id_col_index = line.index("gene_id")
                    transcript_id_col_index = line.index("transcript_id")

                # Check that the transcript got assigned to gene 1 and 
                # transcript 1
                else:
                    assert line[gene_id_col_index] == "1"
                    assert line[transcript_id_col_index] == "1"

                line_num += 1

    @pytest.mark.incremental
    def test_TALON_second_simple_run(self):
        """ After running TALON on the KRT17 example, run the same transcript
            through under a different dataset name. """ 

        try:
             subprocess.check_output(
                     ["python", "../talon.py",
                      "--f", "input_files/KRT17_test_case/trial2_config.csv",
                      "-a", "scratch/KRT17.db",
                      "-b", "hg38",
                      "--o", "scratch/trial2"])
        except:
            pytest.fail("TALON trial 2 failed on KRT17 test")

    @pytest.mark.incremental
    def test_TALON_second_simple_run_check_counters(self): 
        """ After running TALON twice on the same input sam file (different
            dataset label), check whether the counters match the values we 
            expect. """

        conn = sqlite3.connect("scratch/KRT17.db")
        cursor = conn.cursor()
        counter = TALON.get_counters(cursor) 

        assert counter["genes"] == 1
        assert counter["transcripts"] == 1
        assert counter["vertices"] == 16
        assert counter["edges"] == 15
        assert counter["datasets"] == 2
        assert counter["observed_starts"] == 2
        assert counter["observed_ends"] == 2

 
