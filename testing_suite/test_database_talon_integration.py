# Compared to a lot of the other tests in the suite, this one is intended
# to make sure that the major parts of TALON (database init, dataset
# addition over time) are working correctly.

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

class TestDatabaseTalonIntegration(object):
    def test_db_initialization(self):
        """ Initializes a TALON database from a very minimal GTF file. 
            Outfiles are written to the scratch area of the testing suite. 
        """
        try:
            subprocess.check_output(
                ["python", "../initialize_talon_database.py", 
                 "--f", "input_files/KRT17_test_case/KRT17.gtf",  
                 "--a",  "KRT17_test",  
                 "--g",  "hg38", "--o", "scratch/KRT17"])
        except Exception as e:
            pytest.fail("Database initialization failed on KRT17 test")

        try:
            subprocess.check_output(
                ["python", "../initialize_talon_database.py",
                 "--f", "input_files/known_and_novel_test_case/known_and_novel_test_case.gtf",
                 "--a",  "test",
                 "--g",  "hg38", "--o", "scratch/known_and_novel_test_case"])
        except Exception as e:
            pytest.fail("Database initialization failed on known_and_novel test case")

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

    @pytest.mark.incremental
    def test_TALON_known_novel_run(self):
        """ This is a more complex test of TALON functionality. Four transcripts
            must be annotated. Two of them are known, and two are the same novel
            transcript (of a known gene). This function is mostly looking for
            whether there is a crash or not."""

        try:
             subprocess.check_output(
                     ["python", "../talon.py",
                      "--f", "input_files/known_and_novel_test_case/config.csv",
                      "-a", "scratch/known_and_novel_test_case.db",
                      "-b", "hg38",
                      "--o", "scratch/known_and_novel"])
        except:
            pytest.fail("TALON failed on known_novel test case") 

    @pytest.mark.incremental
    def test_TALON_known_novel_correctness(self):
        """ Once we've established that TALON ran on the known-novel example without
            a crash, we need to check that the correct identities were assigned
            to the transcript. """

        infile = "scratch/known_and_novel_talon.tsv"

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
                    if line_num == 1:
                        assert line[gene_id_col_index] == "1"
                        assert line[transcript_id_col_index] == "18"
                    if line_num == 2 or line_num == 3:
                        assert line[gene_id_col_index] == "2"
                        assert line[transcript_id_col_index] == "41"
                    if line_num == 4:
                        assert line[gene_id_col_index] == "3"
                        assert line[transcript_id_col_index] == "13"
                line_num += 1

    @pytest.mark.incremental
    def test_TALON_second_known_novel_run(self):
        """ After running TALON on the known_novel  example, run two more novel
            transcripts through under a different dataset name. """

        try:
             subprocess.check_output(
                     ["python", "../talon.py",
                      "--f", "input_files/known_and_novel_test_case/config2.csv",
                      "-a", "scratch/known_and_novel_test_case.db",
                      "-b", "hg38",
                      "--o", "scratch/known_and_novel_2"])
        except:
            pytest.fail("TALON failed on phase 2 of known_novel test")

    @pytest.mark.incremental
    def test_post_TALON_filter(self):
        """ Take the known_novel database and attempt to filter the transcripts.
            Known transcripts should be accepted, as should novel transcripts
            appearing in more than one dataset. Reject otherwise. """

        database = "scratch/known_and_novel_test_case.db"
        annot = "test"
        expected_transcripts = [ "18", "41", "13"]
        transcript_tuples = filt.filter_talon_transcripts(database, annot)
        transcripts = [ x[1] for x in transcript_tuples ]

        assert sorted(transcripts) == sorted(expected_transcripts)

    #--------- GTF testing section, which requires a TALON database ------

    @pytest.mark.incremental
    def test_gtf_get_gene_2_transcripts(self):
        """ Make sure that all genes and transcripts get included in this
            data structure that maps each TALON gene ID to its transcript
            tuples (from database). Also ensure that transcripts are
            sorted based on start position (ascending). """

        database = "scratch/known_and_novel_test_case.db"
        genome_build = "hg38"
        whitelist = [ str(x) for x in range(1,43) ]

        gene_2_transcripts = createGTF.get_gene_2_transcripts(database, 
                                                              genome_build, 
                                                              whitelist)

        # Check if the keys (genes) match what's expected
        assert sorted(gene_2_transcripts.keys()) == [str(x) for x in range(1,5)] 
        
        # For each gene, check that the transcripts expected are present and 
        # in the correct order
        for gene in sorted(gene_2_transcripts.keys()):
            transcript_tuples = gene_2_transcripts[gene]
            transcript_ids = [ x[1] for x in transcript_tuples ]
            if gene == "1":
                expected_transcripts = [15,26,36,5,27,33,29,12,24,38,34,14,18,7,23,35,16,11,39,9]
                assert transcript_ids == expected_transcripts
            if gene == "2":
                expected_transcripts = [25,6,8,40,17,32,19,28,41,37,21,20,4,22,31,1]
                assert transcript_ids == expected_transcripts
            if gene == "3":
                assert transcript_ids == [13,30,10,2,3]
            if gene == "4":
                assert transcript_ids == [42]

    @pytest.mark.incremental
    def test_validate_GTF(self):
        """ Run code to create a GTF from the test example, then validate
            the GTF file by running ..."""

        # Create a new database from the GTF
        try:
            subprocess.check_output(
                ["python", "../initialize_talon_database.py",
                 "--f", "input_files/gtf_database_test/test.gtf",
                 "--a",  "test", "--l", "0",
                 "--g",  "hg38", "--o", "scratch/gtftest"])
        except:
            pytest.fail("TALON database init failed on known_novel test case")

        # Create the GTF
        try:
            subprocess.check_output(
                 ["python", "../post-TALON_tools/create_GTF_from_database.py",
                  "--db", "scratch/gtftest.db",
                  "-b", "hg38", "-a", "test", "--o", "scratch/gtftest"])
        except:
            pytest.fail("GTF creation code crashed")

        # Check whether Bedtools sort runs on the GTF without crashing
        try:
            subprocess.check_output(
                ["bedtools", "sort", "-i", "scratch/gtftest_talon.gtf"])
        except:
            pytest.fail("Bedtools crashed while trying to sort GTF- likely a formatting problem")

        # Check to see that the GTF we generated matches the original, pre-TALON
        # GTF that was used to build the database originally (known transcripts
        # only).
        # ---------------------------------------------------------------------
        # Check the first 8 fields first for transcripts and exons. These have
        # very fixed start and endpoints.
        subset_gtf(["transcript", "exon"], range(1,9), 
                  "input_files/gtf_database_test/test.gtf",
                  "scratch/original_coords.txt")    
        subset_gtf(["transcript", "exon"], range(1,9),
                    "scratch/gtftest_talon.gtf", "scratch/talon_coords.txt")

        diff = subprocess.check_output(["diff", "scratch/original_coords.txt", 
                                        "scratch/talon_coords.txt"])
        assert diff == ""

        # Check the first 8 fields (minus column 5) for genes. The ends
        # will differ because I compute them based on all of the transcripts
        # of each gene 
        subset_gtf(["gene"], [1,2,3,4,6,7,8], "input_files/gtf_database_test/test.gtf",
                  "scratch/original_gene_coords.txt")
        subset_gtf(["gene"], [1,2,3,4,6,7,8], "scratch/gtftest_talon.gtf", 
                   "scratch/talon_gene_coords.txt")
     
        gdiff = subprocess.check_output(["diff", "scratch/original_gene_coords.txt", 
                                        "scratch/talon_gene_coords.txt"])   
        assert gdiff == ""

        # Check that the attributes before and after match (ignoring attribute
        # fields with 'talon' in the name, of course)
        subset_gtf(["gene", "transcript", "exon"], range(1,10),
                  "input_files/gtf_database_test/test.gtf",
                  "scratch/sorted_original_gtf.txt")
        subset_gtf(["gene", "transcript", "exon"], range(1,10),
                    "scratch/gtftest_talon.gtf", "scratch/sorted_talon_gtf.txt")

        orig = open("scratch/sorted_original_gtf.txt", 'r')
        talon = open("scratch/sorted_talon_gtf.txt", 'r')
        for line_o,line_t in zip(orig, talon):        
            attributes_before = create_attribute_dict(line_o)
            attributes_after = create_attribute_dict(line_t)
            assert attributes_before == attributes_after
        orig.close()
        talon.close()

def subset_gtf(features, cols, infile, outfile):
    """ Selects lines of GTF that match features, and trims to include only
        the columns provided. Writes clipped lines in sorted order to outfile"""

    feat_str = 'if($3 == ' + "|| $3 ==".join(['"' + x + '"' for x in features]) + ")"
    col_str = "$" + ",$".join([str(x) for x in cols])

    cmd = " ".join(["awk -F\t '{", feat_str, "print", col_str, "}'", infile, 
                    "| sort >", outfile])

    os.system(cmd)
    return

def create_attribute_dict(gtf_line):
    """ Creates a dictionary from the key-value GTF attribute pairs for the
        provided GTF line. Ignores attributes with the string 'talon' in them. 
    """

    attributes = {}
    att_col = gtf_line.strip().split("\t")[8]
    atts = att_col.split(";").strip()
    for att_pair in atts:
        key, value = att_pair.split(" ")
        if "talon" in key:
            continue

        attributes[key] = value

    return attributes
            
