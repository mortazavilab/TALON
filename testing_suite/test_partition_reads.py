import pytest
import pysam
from talon import process_sams as procsam
@pytest.mark.unit

class TestPartitionReads(object):
    def test_partition_reads(self):
        """ Given a set of toy transcripts, find the minimum number of
            non-overlapping intervals to contain them. 
            Expected read groups:
            - read_1 and read_2 (they are in interval chr1:1-1000)
            - read_3 (interval chr2:1-100)
        """

        sams = ["input_files/preprocess_sam/read1.sam",
                "input_files/preprocess_sam/read2.sam"]
        datasets = ["dataset1", "dataset2"]
        tmp_dir = "scratch/test_read_partition/"

        read_groups, intervals, merged_bam = procsam.partition_reads(sams, datasets, tmp_dir = tmp_dir)

        # Check length and membership of read groups
        assert len(read_groups) == 2 

        assert [ entry.query_name for entry in read_groups[0] ] == ["read_1", "read_2"]
        assert [ entry.query_name for entry in read_groups[1] ] == ["read_3"]
        
        # Check the intervals 
        assert intervals[0] == ("chr1", 1, 1004)
        assert intervals[1] == ("chr2", 1, 100)
