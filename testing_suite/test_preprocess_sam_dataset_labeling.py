import pytest
import pysam
from talon import process_sams as procsam
@pytest.mark.unit

class TestMergeReads(object):
    def test_read_labels(self):
        """ Given two sam files and two corresponding dataset labels, ensure 
            that the merged BAM file preserves the RB tags assigned by 
            pysam.merge """

        sams = ["input_files/preprocess_sam/read1.sam", 
               "input_files/preprocess_sam/read2.sam" ]
        datasets = ["dataset1", "dataset2"]
        tmp_dir = "scratch/test_read_labels/"

        merged_bam = procsam.preprocess_sam(sams, datasets, tmp_dir = tmp_dir)

        with pysam.AlignmentFile(merged_bam) as bam:
            for entry in bam:
                if entry.query_name == "read1":
                    assert entry.get_tag("RG") == "dataset1"
                elif entry.query_name == "read2":
                    assert entry.get_tag("RG") == "dataset2"
                else:
                    pytest.fail("Unexpected read encountered")   
            
