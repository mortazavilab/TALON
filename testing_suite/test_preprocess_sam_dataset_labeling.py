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
        tmp_dir = "scratch/test_read_labels/test1/"

        merged_bam = procsam.preprocess_sam(sams, datasets, tmp_dir = tmp_dir, use_cb_tag = False)

        with pysam.AlignmentFile(merged_bam) as bam:
            for entry in bam:
                if entry.query_name == "read_1":
                    assert entry.get_tag("RG") == "dataset1"
                elif entry.query_name == "read_2":
                    assert entry.get_tag("RG") == "dataset1"
                elif entry.query_name == "read_3":
                    assert entry.get_tag("RG") == "dataset2"
                else:
                    pytest.fail("Unexpected read encountered")

    def test_unsorted_sam_file(self):
        """ Make sure that the function can handle a SAM input that has not
            been sorted """

        sams = ["input_files/chr11_and_Tcf3/BC017.sam"]
        datasets = ["BC017"]
        tmp_dir = "scratch/test_read_labels/test2/"

        merged_bam = procsam.preprocess_sam(sams, datasets, tmp_dir = tmp_dir, use_cb_tag = False)
        with pysam.AlignmentFile(merged_bam) as bam:
            # Check the first read
            for entry in bam:
                assert entry.query_name == "m54284_180814_002203/19268005/ccs"
                break
