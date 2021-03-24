# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Functions related to processing the input SAM files and partitioning them
# for processing in parallel

import pybedtools
import pysam
import os
import time

def convert_to_bam(sam, bam):
    """ Convert provided sam file to bam file (provided name).  """
    
    try:
        infile = pysam.AlignmentFile(sam, "r")
        outfile = pysam.AlignmentFile(bam, "wb", template=infile)
        for s in infile:
            outfile.write(s)

    except Exception as e:
        print(e)
        raise RuntimeError("Problem converting sam file '%s' to bam." % (sam))
             

def preprocess_sam(sam_files, datasets, tmp_dir = "/dev/shm/talon/", n_threads = 0):
    """ Copy and rename the provided SAM/BAM file(s), merge them, and index.
        This is necessary in order to use Pybedtools commands on the reads.
        The renaming is necessary in order to label the reads according to
        their dataset."""

    # Create the tmp dir
    os.system("mkdir -p %s " % (tmp_dir))

    # Copy and rename SAM files with dataset names to ensure correct RG tags
    renamed_sams = []
    for sam, dataset in zip(sam_files, datasets):
        suffix = "." + sam.split(".")[-1]
        if suffix == ".sam":
            bam_copy = tmp_dir + dataset + "_unsorted.bam"
            convert_to_bam(sam, bam_copy)
            sam = bam_copy
        sorted_bam = tmp_dir + dataset + ".bam"
        pysam.sort("-@", str(n_threads), "-o", sorted_bam, sam)
        renamed_sams.append(sorted_bam)

    merged_bam = tmp_dir + "merged.bam"
    merge_args = [merged_bam] + renamed_sams + ["-f", "-r", "-@", str(n_threads)]
    # index_args = [merged_bam, "-@", str(n_threads)]

    # Merge datasets and use -r option to include a read group tag
    try:
        pysam.merge(*merge_args)
        pysam.index(merged_bam)
        ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        print("[ %s ] Merged input SAM/BAM files" % (ts))
    except:
        raise RuntimeError(("Problem merging and indexing SAM/BAM files. "
                            "Check your file paths and make sure that all "
                            "files have headers."))
    return merged_bam

def partition_reads(sam_files, datasets, tmp_dir = "/dev/shm/talon/", n_threads = 0):
    """ Use bedtools merge to create non-overlapping intervals from all of the
        transcripts in a series of SAM/BAM files. Then, iterate over the intervals
        to extract all reads inside of them from the pysam object.
       
        Returns:
            - List of lists: sublists contain pysam reads from a given interval
            - List of tuple intervals
            - filename of merged bam file (to keep track of the header)
           """

    merged_bam = preprocess_sam(sam_files, datasets, tmp_dir = tmp_dir, 
                                n_threads = n_threads)

    try:
        all_reads = pybedtools.BedTool(merged_bam).bam_to_bed()
    except Exception as e:
        print(e)
        raise RuntimeError("Problem opening sam file %s" % (merged_bam))

    # Must sort the Bedtool object
    sorted_reads = all_reads.sort()
    intervals = sorted_reads.merge(d = 100000000)

    # Now open each sam file using pysam and extract the reads
    coords = []
    read_groups = []
    with pysam.AlignmentFile(merged_bam) as bam:  # type: pysam.AlignmentFile
        for interval in intervals:
            reads = get_reads_in_interval(bam, interval.chrom,
                                          interval.start, interval.end)
            read_groups.append(reads)
            coords.append((interval.chrom, interval.start + 1, interval.end))

    return read_groups, coords, merged_bam

def write_reads_to_file(read_groups, intervals, header_template, tmp_dir = "/dev/shm/talon/"):
    """ For each read group, iterate over the reads and write them to a file
        named for the interval they belong to. This step is necessary because
        Pysam objects cannot be pickled. """
 
    tmp_dir = tmp_dir + "interval_files/"
    os.system("mkdir -p %s " % (tmp_dir))
 
    outbam_names = []
    with pysam.AlignmentFile(header_template, "rb") as template:
        for group, interval in zip(read_groups, intervals):
            fname = tmp_dir + "_".join([str(x) for x in interval]) + ".bam"
            outbam_names.append(fname)
            with pysam.AlignmentFile(fname, "wb", template = template) as f:
                for read in group:
                    f.write(read)
        
    return outbam_names


def get_reads_in_interval(sam, chrom, start, end):
    """ Given an open pysam.AlignmentFile, return only the reads that overlap
        the provided interval. Note that this means there may be reads that
        extend beyond the bounds of the interval. """
    iterator = sam.fetch(chrom, start, end)
    reads = [ x for x in iterator ]
    return reads

