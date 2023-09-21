import pytest
import sqlite3
from .helper_fns import fetch_correct_ID

@pytest.mark.integration

# All data comes from hl60_1_1 from the ENCODE data

class TestAssignments(object):
    """ The objective here is to make sure that each transcript in the
        readthrough example set was assigned the expected identity. """

    def test_FSM_of_annot_rt(self):
        """ cenps_cort_fsm is a FSM to the annotated readthrough locus of
            CENPS-CORT. Comes from ENCODE hl60 data"""

        conn = sqlite3.connect("scratch/readthrough.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "hl60_1_1"
        read_ID = "cenps_cort_fsm"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("CENPS-CORT", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID

        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["transcript_status"] == "KNOWN"
        conn.close()

    def test_ISM_of_annot_rt(self):
       """ cenps_cort_ism is an ISM of readthrough locus of CENPS-CORT"""

       conn = sqlite3.connect("scratch/readthrough.db")
       conn.row_factory = sqlite3.Row
       cursor = conn.cursor()

       dataset = "hl60_1_1"
       read_ID = "cenps_cort_ism"

       # Fetch observed entry from table
       query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
       assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

       correct_gene_ID = fetch_correct_ID("CENPS-CORT", "gene", cursor)
       assert assignment['gene_ID'] == correct_gene_ID

       # Now make sure that the novel transcript was annotated correctly
       annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
       assert annot_dict["ISM_transcript"] == "TRUE"
       conn.close()

    def test_NNC_of_annot_rt(self):
       """ cenps_cort_nnc shares most sjs with CENPS-CORT redthrough locus """

       conn = sqlite3.connect("scratch/readthrough.db")
       conn.row_factory = sqlite3.Row
       cursor = conn.cursor()

       dataset = "hl60_1_1"
       read_ID = "cenps_cort_nnc"

       # Fetch observed entry from table
       query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
       assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]
       correct_gene_ID = fetch_correct_ID("CENPS-CORT", "gene", cursor)
       assert assignment['gene_ID'] == correct_gene_ID

       # Now make sure that the novel transcript was annotated correctly
       annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
       assert annot_dict["NNC_transcript"] == "TRUE"
       assert annot_dict["transcript_status"] == "NOVEL"
       conn.close()

    def test_NIC_of_annot_rt(self):
        """ cenps_cort_nnc shares all ss, but has a novel sj with CENPS-CORT redthrough locus """

        conn = sqlite3.connect("scratch/readthrough.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "hl60_1_1"
        read_ID = "cenps_cort_nic"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]
        correct_gene_ID = fetch_correct_ID("CENPS-CORT", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["NIC_transcript"] == "TRUE"
        assert annot_dict["transcript_status"] == "NOVEL"
        conn.close()

    def test_FSM_of_novel_rt_1(self):
        """ eloa_rpl11_fsm_1 is FSM to 2 different genes, ELOA and RPL11"""
        conn = sqlite3.connect("scratch/readthrough.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "hl60_1_1"
        read_ID = "eloa_rpl11_fsm_1"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        # we had 5 annotated genes (CENPS, CORT, CENPS-CORT, ELOA, and RPL11)
        # so new gene should be 6
        correct_gene_ID = 6
        assert assignment['gene_ID'] == correct_gene_ID

        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict['transcript_status'] == 'NOVEL'
        assert annot_dict["fusion_transcript"] == "TRUE"

        annot_dict = make_annot_dict_gene(cursor, assignment['gene_ID'])
        assert annot_dict['gene_status'] == 'NOVEL'
        assert annot_dict['fusion_novel'] == 'TRUE'

        conn.close()

    def test_FSM_of_novel_rt_2(self):
        """ eloa_rpl11_fsm_2 is FSM to 2 different genes, ELOA and RPL11"""
        conn = sqlite3.connect("scratch/readthrough.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "hl60_1_1"
        read_ID = "eloa_rpl11_fsm_2"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        # we had 5 annotated genes (CENPS, CORT, CENPS-CORT, ELOA, and RPL11)
        # so new gene should be 6
        correct_gene_ID = 6
        assert assignment['gene_ID'] == correct_gene_ID

        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict['transcript_status'] == 'NOVEL'
        assert annot_dict["ISM_transcript"] == "TRUE"
        assert annot_dict['ISM-suffix_transcript'] == 'TRUE'

        annot_dict = make_annot_dict_gene(cursor, assignment['gene_ID'])
        assert annot_dict['gene_status'] == 'NOVEL'
        assert annot_dict['fusion_novel'] == 'TRUE'

        conn.close()

    def test_FSM_of_overlapping_single_gene(self):
        """ rpl11_fsm is an FSM an annotated gene that is subsumed by the
            RPL11-ELOA readthrough loci. However it should just
            be annotated to RPL11"""

        conn = sqlite3.connect("scratch/readthrough.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "hl60_1_1"
        read_ID = "rpl11_fsm"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("RPL11", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID

        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["transcript_status"] == "KNOWN"
        conn.close()

    def test_ISM_of_overlapping_single_gene(self):
        """ rpl11_ism is an ISM an annotated gene that is subsumed by the
            RPL11-ELOA readthrough loci. However it should just
            be annotated to RPL11"""

        conn = sqlite3.connect("scratch/readthrough.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "hl60_1_1"
        read_ID = "rpl11_ism"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]

        correct_gene_ID = fetch_correct_ID("RPL11", "gene", cursor)
        assert assignment['gene_ID'] == correct_gene_ID

        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict['transcript_status'] == 'NOVEL'
        assert annot_dict["ISM_transcript"] == "TRUE"
        assert annot_dict['ISM-suffix_transcript'] == 'TRUE'
        conn.close()

    def test_NNC_of_annot_rt(self):
       """ eloa_rpl11_nnc shares most sjs with novel ELOA-RPL11 rt locus """

       conn = sqlite3.connect("scratch/readthrough.db")
       conn.row_factory = sqlite3.Row
       cursor = conn.cursor()

       dataset = "hl60_1_1"
       read_ID = "eloa_rpl11_nnc"

       # Fetch observed entry from table
       query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
       assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]
       correct_gene_ID = 6
       assert assignment['gene_ID'] == correct_gene_ID

       # Now make sure that the novel transcript was annotated correctly
       annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
       assert annot_dict["NNC_transcript"] == "TRUE"
       assert annot_dict["transcript_status"] == "NOVEL"
       conn.close()

    def test_NNC_of_annot_rt(self):
        """ eloa_rpl11_nic shares all sjs but one new sj w/ novel ELOA-RPL11 rt locus """

        conn = sqlite3.connect("scratch/readthrough.db")
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        dataset = "hl60_1_1"
        read_ID = "eloa_rpl11_nic"

        # Fetch observed entry from table
        query = """SELECT * from observed WHERE dataset = ? AND read_name = ?"""
        assignment = cursor.execute(query, [dataset, read_ID]).fetchall()[0]
        correct_gene_ID = 6
        assert assignment['gene_ID'] == correct_gene_ID

        # Now make sure that the novel transcript was annotated correctly
        annot_dict = make_annot_dict(cursor, assignment['transcript_ID'])
        assert annot_dict["NIC_transcript"] == "TRUE"
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
