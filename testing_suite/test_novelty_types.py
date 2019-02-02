# This set of tests is designed to make sure that novelty designations of 
# transcripts are being handled correctly.

import pytest
import os
import subprocess
import sqlite3
import sys
sys.path.append("..")
import talon as TALON
sys.path.append("../post-TALON_tools")
import filter_talon_transcripts as filt
import create_GTF_from_database as createGTF
@pytest.mark.integration

class TestNoveltyTypes(object):
    def test_db_initialization_toy(self):
        """ Initializes a TALON database from a very minimal GTF file. 
            Outfiles are written to the scratch area of the testing suite. 
        """
        try:
            subprocess.check_output(
                ["python", "../initialize_talon_database.py", 
                 "--f", "input_files/toy_transcript/toy_annot.gtf",  
                 "--a",  "toy_annot",  
                 "--g",  "toy_build", "--o", "scratch/toy"])
        except Exception as e:
            pytest.fail("Database initialization failed on toy artificial transcript")


    @pytest.mark.incremental
    def test_TALON_toy_ISMs(self):
        """ Once the database has been initialized, try running TALON on toy 
            incomplete splice match transcripts (that we know match the two
            final exons of the annotated toy transcript). First,
            establish that TALON runs without crashing."""

        try:
             subprocess.check_output(
                     ["python", "../talon.py",
                      "--f", "input_files/toy_transcript/config.csv",
                      "-a", "scratch/toy.db",
                      "-b", "toy_build",
                      "-l", "0",
                      "--o", "scratch/toy"])
        except:
            pytest.fail("TALON run failed on toy data")

    @pytest.mark.incremental
    def test_assign_novelty_type(self):
        """ For each toy read that we processed, compute the novelty type and
            compare it to the expected assignment. 
        """
        pass

def is_suffix_ISM(transcript_ID, database):
    """ For now, this function checks whether a transcript is an incomplete 
        splice match (ISM) and a suffix of a known transcript."""
    
    # Connect to the database
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    # Formulate query
    query = """ SELECT * FROM observed WHERE transcript_ID = %s"""
    print query
    cursor.execute(query % transcript_ID)
    print cursor.fetchall()


