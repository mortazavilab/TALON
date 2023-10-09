# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Functions related to processing the input SAM files and partitioning them
# for processing in parallel

import pyranges as pr
import pysam
import os
import time
import logging

save = pysam.set_verbosity(0)
# pysam.set_verbosity(save)

def convert_to_bam(sam, bam, threads):
    """ Convert provided sam file to bam file (provided name).  """

    try:
        infile = pysam.AlignmentFile(sam, "r", threads=threads)
        outfile = pysam.AlignmentFile(
            bam, "wb", template=infile, threads=threads)
        for s in infile:
            outfile.write(s)

    except Exception as e:
        logging.error(e)
        msg = f'Problem converting SAM file {sam} to BAM'
        logging.error(msg)
        raise RuntimeError(msg)
        # raise RuntimeError("Problem converting sam file '%s' to bam." % (sam))


def preprocess_sam(sam_files, datasets, use_cb_tag,
                   tmp_dir="talon_tmp/", n_threads=0):
    """ Copy and rename the provided SAM/BAM file(s), merge them, and index.
        This is necessary in order to use following commands on the reads.
        The renaming is necessary in order to label the reads according to
        their dataset."""

    # Create the tmp dir
    os.system("mkdir -p %s " % (tmp_dir))

    renamed_sams = []

    # merge sam files if we're not pulling out dataset IDs from RG tag
    if not use_cb_tag:
        # Copy and rename SAM files with dataset names to ensure correct RG tags
        for sam, dataset in zip(sam_files, datasets):
            suffix = "." + sam.split(".")[-1]
            if suffix == ".sam":
                bam_copy = tmp_dir + dataset + "_unsorted.bam"
                convert_to_bam(sam, bam_copy, n_threads)
                sam = bam_copy
            sorted_bam = tmp_dir + dataset + ".bam"
            pysam.sort("-@", str(n_threads), "-o", sorted_bam, sam)
            renamed_sams.append(sorted_bam)

        merged_bam = tmp_dir + "merged.bam"
        merge_args = [merged_bam] + renamed_sams + \
            ["-f", "-r", "-@", str(n_threads)]
        # index_args = [merged_bam, "-@", str(n_threads)]

        # # Merge datasets and use -r option to include a read group tag
        # try:
        #     pysam.merge(*merge_args)
        #     pysam.index(merged_bam)
        #     ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        #     print("[ %s ] Merged input SAM/BAM files" % (ts))
        # except:
        #     raise RuntimeError(("Problem merging and indexing SAM/BAM files. "
        #                         "Check your file paths and make sure that all "
        #                         "files have headers."))

    # otherwise, merge without adding the RG tag and convert to bam
    elif use_cb_tag:
        for i, sam in enumerate(sam_files):
            fname_split = sam.split(".")
            suffix = "."+fname_split[-1]
            if suffix == ".sam":
                bam_copy = '{}{}_unsorted.bam'.format(tmp_dir, i)
                convert_to_bam(sam, bam_copy, n_threads)
                sam = bam_copy
            sorted_bam = '{}{}.bam'.format(tmp_dir, i)
            pysam.sort("-@", str(n_threads), "-o", sorted_bam, sam)
            renamed_sams.append(sorted_bam)

        merged_bam = tmp_dir + "merged.bam"
        merge_args = [merged_bam] + renamed_sams + ["-f", "-@", str(n_threads)]

    # Merge datasets and use -r option to include a read group tag
    try:
        pysam.merge(*merge_args)
        sorted_bam = tmp_dir + "merged_sorted.bam"
        pysam.sort("-@", str(n_threads), "-o", sorted_bam, merged_bam)
        pysam.index(sorted_bam)
        # ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        # print("[ %s ] Merged input SAM/BAM files" % (ts))
        logging.info('Merged input SAM/BAM files')
    except:
        # raise RuntimeError(("Problem merging and indexing SAM/BAM files. "
        #                     "Check your file paths and make sure that all "
        #                     "files have headers."))
        msg = "Problem merging and indexing SAM/BAM files. "+\
                            "Check your file paths and make sure that all "+\
                            "files have headers."
        logging.error(msg)
        raise RuntimeError(msg)
    return sorted_bam


def partition_reads(sam_files, datasets, use_cb_tag,
                    tmp_dir="talon_tmp/", n_threads=0):
    """ Use bedtools merge to create non-overlapping intervals from all of the
        transcripts in a series of SAM/BAM files. Then, iterate over the intervals
        to extract all reads inside of them from the pysam object.

        Returns:
            - List of lists: sublists contain pysam reads from a given interval
            - List of tuple intervals
            - filename of merged bam file (to keep track of the header)
           """
    merged_bam = preprocess_sam(sam_files, datasets, use_cb_tag,
                                tmp_dir=tmp_dir, n_threads=n_threads)

    try:
        gr = pr.read_bam(merged_bam)
    except Exception as e:
        # print(e)
        logging.error(e)
        msg = f'Problem opening SAM file {merged_bam}'
        logging.error(msg)
        raise RuntimeError(msg)

    gr = gr.merge(slack=100000000, strand=False)

    # Now open each sam file using pysam and extract the reads
    coords = []
    read_groups = []
    with pysam.AlignmentFile(merged_bam) as bam:  # type: pysam.AlignmentFile
        for _, interval in gr.df.iterrows():
            reads = get_reads_in_interval(bam, interval.Chromosome,
                                          interval.Start, interval.End)
            read_groups.append(reads)
            coords.append((interval.Chromosome,
                           interval.Start + 1, interval.End))

    return read_groups, coords, merged_bam


def write_reads_to_file(read_groups, intervals,
                        header_template, tmp_dir="talon_tmp/"):
    """ For each read group, iterate over the reads and write them to a file
        named for the interval they belong to. This step is necessary because
        Pysam objects cannot be pickled. """

    tmp_dir = tmp_dir + "interval_files/"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    outbam_names = []
    with pysam.AlignmentFile(header_template, "rb") as template:
        for group, interval in zip(read_groups, intervals):
            fname = tmp_dir + "_".join([str(x) for x in interval]) + ".bam"
            outbam_names.append(fname)
            with pysam.AlignmentFile(fname, "wb", template=template) as f:
                for read in group:
                    f.write(read)

    return outbam_names


def get_reads_in_interval(sam, chrom, start, end):
    """ Given an open pysam.AlignmentFile, return only the reads that overlap
        the provided interval. Note that this means there may be reads that
        extend beyond the bounds of the interval. """
    iterator = sam.fetch(chrom, start, end)
    reads = [x for x in iterator]
    return reads
