import pytest
import sys
import sqlite3
sys.path.append("..")
import talon as talon
from helper_fns import *
@pytest.mark.integration

class TestAssignments(object):
    """ The objective here is to make sure that each transcript in the 
        chr11_and_Tcf3 example set was assigned the expected identity. """

    def test_ISM_of_Canx(self):
        """ m54284_180814_204240/72352410/ccs is an ISM transcript of Canx.
            Comes from BC017 data. """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        read_ID = "m54284_180814_204240/72352410/ccs"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("Canx", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID
        assert assignment['transcript_ID'] == 8451
        assert assignment['start_delta'] == -4
        assert assignment['end_delta'] == -40

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])

        assert annot_dict["ISM_transcript"] == "TRUE"
        assert annot_dict["transcript_status"] == "NOVEL"
        assert annot_dict["ISM_to_IDs"] == "1743,1744"
        conn.close()

    def test_prefix_ISM_of_Canx(self):
        """ m54284_180814_002203/18677911/ccs is an ISM transcript of Canx.
            Comes from BC017 data. """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        read_ID = "m54284_180814_002203/18677911/ccs"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("Canx", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID
        assert assignment['transcript_ID'] == 8452
        assert assignment['start_delta'] == 64
        assert assignment['end_delta'] == -290 #1

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["ISM_transcript"] == "TRUE"
        assert annot_dict["ISM-prefix_transcript"] == "TRUE"
        assert annot_dict["transcript_status"] == "NOVEL"
        assert annot_dict["ISM_to_IDs"] == "1744"
        assert annot_dict["ISM-prefix_to_IDs"] == "1744"
        conn.close()

    def test_genomic_of_Tcf3(self):
        """ m54284_180814_002203/19268005/ccs is a genomic transcript of Tcf3
            from the BC017 data """
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        read_ID = "m54284_180814_002203/19268005/ccs"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("Tcf3", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["genomic_transcript"] == "TRUE"
        assert annot_dict["transcript_status"] == "NOVEL"
        conn.close()

    def test_suffix_ISM_of_Tcf3(self):
        """ m54284_180814_002203/18809472/ccs is an ISM suffix transcript of Tcf3.
            Comes from BC017 data. """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        read_ID = "m54284_180814_002203/18809472/ccs"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("Tcf3", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID
        assert assignment['transcript_ID'] == 8454
        assert assignment['start_delta'] == 77
        assert assignment['end_delta'] == 8

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["ISM_transcript"] == "TRUE"
        assert annot_dict["ISM-suffix_transcript"] == "TRUE"
        assert annot_dict["transcript_status"] == "NOVEL"
        assert annot_dict["ISM_to_IDs"] == "8437,8438,8440,8443,8445,8446"
        assert annot_dict["ISM-suffix_to_IDs"] == "8437,8438,8440,8443,8445,8446"
        conn.close()

    def test_NIC_of_Drg1(self):
        """ For this example, the same read was planted in two different 
            datasets (m54284_180814_002203/49414590/ccs) """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset_1 = "PB65_B017"
        dataset_2 = "PB65_B018"
        read_ID = "m54284_180814_002203/49414590/ccs"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset IN 
                   ('PB65_B017', 'PB65_B018') AND read_name = ?"""
        cursor.execute(query, [read_ID])
        correct_gene_ID = fetch_correct_ID("Drg1", "gene", cursor)
        for assignment in cursor.fetchall():
            assert assignment['gene_ID'] == correct_gene_ID
            assert assignment['transcript_ID'] == 8455
            assert assignment['start_delta'] == 2
            assert assignment['end_delta'] == -15

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, 8455)
        assert annot_dict["NIC_transcript"] == "TRUE"
        assert annot_dict["transcript_status"] == "NOVEL"
        conn.close()

    def test_FSM_of_Drg1(self):
        """ Read m54284_180814_002203/40042763/ccs is an FSM of the Drg1 gene
            (BC017) """
        
        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        read_ID = "m54284_180814_002203/40042763/ccs"
        
        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("Drg1", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID
        assert assignment['transcript_ID'] == 28
        assert assignment['start_vertex'] == 153
        assert assignment['end_vertex'] == 139

        # Now make sure that the novel transcript was annotated correctly
        # TODO: in the future, I would like this transcript to prioritize
        # the known transcript ID 28. Just a thought
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        #assert annot_dict["FSM_transcript"] == "TRUE"
        #assert annot_dict["FSM_to_IDs"] == "28"
        assert annot_dict["transcript_status"] == "KNOWN"
        conn.close()

    def antisense_to_Grb10(self):
        """ Read m54284_180814_002203/8126905/ccs from BC017 """

        conn = sqlite3.connect("scratch/chr11_and_Tcf3.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "PB65_B017"
        read_ID = "m54284_180814_002203/8126905/ccs"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_antisense_gene_ID = fetch_correct_ID("Grb10", "gene", cursor)
        assert assignment['gene_ID'] == 2988
        assert assignment['transcript_ID'] == 8457

        # Now make sure that the novel gene was annotated correctly
        annot_dict_gene = make_annot_dict(cursor, 8457)
        assert annot_dict_gene["antisense_gene"] == "TRUE"
        assert annot_dict_gene["gene_status"] == "NOVEL"
        assert annot_dict_gene["gene_antisense_to_IDs"] == str(correct_antisense_gene_ID)

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["antisense_transcript"] == "TRUE"
        assert annot_dict["transcript_status"] == "NOVEL"
        conn.close()
        
        

def make_annot_dict_gene(cursor, gene_ID):
    """ Extracts all gene annotations for the transcript ID and puts
        them in a dict """
    query = """SELECT * from gene_annotations WHERE ID = ?"""
    annotations = cursor.execute(query, [gene_ID]).fetchall()
    annot_dict = {}
    for annot in annotations:
        attribute = annot["attribute"]
        value = annot["value"]
        annot_dict[attribute] = value
    return annot_dict

def make_annot_dict(cursor, transcript_ID):
    """ Extracts all transcript annotations for the transcript ID and puts
        them in a dict """
    query = """SELECT * from transcript_annotations WHERE ID = ?"""
    annotations = cursor.execute(query, [transcript_ID]).fetchall()
    annot_dict = {}
    for annot in annotations:
        attribute = annot["attribute"]
        value = annot["value"]
        annot_dict[attribute] = value
    return annot_dict
