import pysam
from talon import talon_label_reads as tlr
import pytest
import os

def test_invalid_input():
    """ Should crash when given a file that is not SAM/BAM """

    infile = "talon_label_reads/test_inputs/toy_genome.fa" 
    with pytest.raises(ValueError): 
        tlr.split_reads_by_chrom(infile, tmp_dir = "scratch/tlr/raw")

def test_split_reads_by_chrom_sam_input():
    """ When a SAM file is provided, the function should convert it to BAM,
        sort and index it, and then split it into chromosome-wise files.
        Chromosomes with no reads should be omitted. """

    infile = "talon_label_reads/test_inputs/test_split_by_chrom/sample_reads.sam" 
    tmp_dir = "scratch/tlr/sam"
    tlr.split_reads_by_chrom(infile, tmp_dir = tmp_dir)

    assert os.path.isfile(tmp_dir + "/raw/all_reads.sorted.bam")
    assert os.path.isfile(tmp_dir + "/raw/all_reads.sorted.bam" + ".bai")

    chrom_dir = tmp_dir + "/raw/chroms"
    # Make sure chr3 was left out- it had no reads
    assert os.path.isfile(chrom_dir + "/chr3.sam") == False

    # Check results for chroms 1 and 2
    check_result(chrom_dir + "/chr1.sam", "chr1", 2)
    check_result(chrom_dir + "/chr2.sam", "chr2", 1)

def test_split_reads_by_chrom_bam_input():
    """ When a BAM file is provided, the function should index it if necessary
        and then split it into chromosome-wise files.
        Chromosomes with no reads should be omitted. The content in this BAM 
        input file is the same as for the SAM test. """

    infile = "talon_label_reads/test_inputs/test_split_by_chrom/sample_reads.bam"
    tmp_dir = "scratch/tlr/bam"
    tlr.split_reads_by_chrom(infile, tmp_dir = tmp_dir)

    chrom_dir = tmp_dir + "/raw/chroms"
    # Make sure chr3 was left out- it had no reads
    assert os.path.isfile(chrom_dir + "/chr3.sam") == False

    # Check results for chroms 1 and 2
    check_result(chrom_dir + "/chr1.sam", "chr1", 2)
    check_result(chrom_dir + "/chr2.sam", "chr2", 1)

def check_result(sam_file, chrom, expected_reads):
    """ Make sure that the reads in the provided SAM file are on the right
        chromosome in the expected quantity."""
    count = 0
    with pysam.AlignmentFile(sam_file, "r") as sam:
        for record in sam:
            assert record.reference_name == chrom 
            count += 1
    assert count == expected_reads
