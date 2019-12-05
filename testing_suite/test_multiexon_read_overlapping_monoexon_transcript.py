import pytest
from talon import talon, init_refs
import sqlite3
import pysam
@pytest.mark.integration

class TestMultExonReadOverlapsMonoGene(object):

    def test_transcript_assigned_intergenic(self):
        """ This test covers a case reported by a user where a read overlaps
            the ~600bp mono-exonic pseudogene HMGB1P1. The read itself has
            2 exons, the second of which contains the small pseudogene inside.
            earlier versions of TALON classified the read as intergenic,
            when it was actually supposed to be genomic """

        # Set up references
        database = "scratch/multiexon_read_overlapping_monoexon_transcript/talon.db" 
        conn = sqlite3.connect(database)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        build = "hg38"
        talon.get_counters(database)
        run_info = talon.init_run_info(database, build)
        struct_collection = talon.prepare_data_structures(cursor, run_info)

        # Use pysam to get the read from the SAM file
        sam_file = "input_files/multiexon_read_overlapping_monoexon_transcript/read.sam"
        with pysam.AlignmentFile(sam_file) as sam:
            for entry in sam:
                sam_record = entry
                break

        # Get read attributes
        chrom = sam_record.reference_name
        strand = "-" if sam_record.is_reverse else "+"
        sam_start = sam_record.reference_start 
        sam_end = sam_record.reference_end

        # Do we get any overlap with the reference gene?
        best_gene, match_strand = talon.search_for_overlap_with_gene(chrom, min(sam_start, sam_end),
                                                                     max(sam_start, sam_end), strand, 
                                                                     cursor, run_info, 
                                                                     struct_collection.tmp_gene)
        assert best_gene == 1
        assert match_strand == "-"

        annotation_info = talon.annotate_read(sam_record, cursor, run_info, 
                                              struct_collection, mode = 0)
        
        assert annotation_info['gene_ID'] == 1
        assert annotation_info['transcript_ID'] == 2
        assert 'genomic_transcript' in annotation_info['transcript_novelty'][0]

