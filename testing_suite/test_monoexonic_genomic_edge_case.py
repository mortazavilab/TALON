import pytest
import sqlite3
#from talon import talon, init_refs
#from .helper_fns import fetch_correct_ID
@pytest.mark.integration

class TestMonoexonicGenomicEdgeCase(object):
    """ The aim here is to make sure that TALON does not apply the same type of
        edge constraint to monoexonic transcripts of monoexonic genes that it 
        does to spliced alignments. This can result in spurious genomic transcripts.
    """
    def test_within_bounds(self):
        """ Example where the monoexonic start and end positions are within the
            allowable cutoff distances. Expectation is for the transcript to
            be called as known. TALON run is done as part of the build_dataabases
            step, so we simply check the results."""

        with sqlite3.connect("scratch/monoexon.db") as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()

            cursor.execute("""SELECT * FROM observed 
                              WHERE read_name = "read_1" """)
            result = cursor.fetchone()
            print([x for x in result])
            assert result['gene_ID'] == 1
            assert result['transcript_ID'] == 1
            assert result['start_delta'] == 10
            assert result['end_delta'] == -91
         
    def test_out_of_bounds(self):
        """ Example where the distance at the ends is > the cutoff of 99 bp.
        """
    
        with sqlite3.connect("scratch/monoexon.db") as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()

            cursor.execute("""SELECT * FROM observed
                              WHERE read_name = "read_2" """)
            result = cursor.fetchone()
            print([x for x in result])
            assert result['gene_ID'] == 1
            assert result['transcript_ID'] == 1
            assert result['start_delta'] == 100
            assert result['end_delta'] == -401
