# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Useful functions for processing transcripts in the SAM format

import itertools
import re
import pysam

def check_read_quality(sam_record: pysam.AlignedSegment, run_info):
    """ Process an individual sam read and return quality attributes. """
    read_ID = sam_record.query_name
    flag = sam_record.flag
    cigar = sam_record.cigarstring
    seq = sam_record.query
    read_length = sam_record.query_length
    dataset = sam_record.get_tag('RG')
 
    # Only use uniquely mapped transcripts
    if flag not in [0, 16]:
        return [dataset, read_ID, 0, 0, read_length, "NA", "NA"]

    # Only use reads that are greater than or equal to length threshold
    if read_length < run_info.min_length:
        return [dataset, read_ID, 0, 1, read_length, "NA", "NA"]

    # Locate the MD field of the sam transcript
    try:
        md_tag = sam_record.get_tag('MD')
    except KeyError:
        raise ValueError("SAM transcript %s lacks an MD tag" % read_ID)

    # Only use reads where alignment coverage and identity exceed
    # cutoffs
    coverage = compute_alignment_coverage(cigar)
    identity = compute_alignment_identity(md_tag, seq)

    if coverage < run_info.min_coverage or \
       identity < run_info.min_identity:
        return [dataset, read_ID, 0, 1, read_length, coverage, identity]

    # At this point, the read has passed the quality control
    return [dataset, read_ID, 1, 1, read_length, coverage, identity]

def compute_alignment_coverage(CIGAR):
    """ This function computes what fraction of the read is actually aligned to
        the genome by excluding hard or soft-clipped bases."""

    total_bases = 0.0
    unaligned_bases = 0.0
    ops, counts = split_cigar(CIGAR)
    for op,ct in zip(ops, counts):
        if op == "N":
            continue
        if op == "H" or op == "S":
            unaligned_bases += ct
        total_bases += ct

    return (total_bases - unaligned_bases)/total_bases

def compute_alignment_identity(MD_tag, SEQ):
    """ This function computes what fraction of the read matches the reference
        genome."""

    total_bases = len(SEQ)
    matches = 0.0
    ops, counts = splitMD(MD_tag)
    for op,ct in zip(ops, counts):
        if op == "M":
            matches += ct
        if op == "D":
            total_bases += ct

    return matches/total_bases

def splitMD(MD):
        """ Takes MD tag and splits into two lists:
            one with capital letters (match operators), and one with
            the number of bases that each operation applies to. """

        operations = []

        # Split MD string where type changes.
        # Digits are separated from base changes.
        # Deletions (with ^) are captured together.
        counts = ["".join(x) for _, x in itertools.groupby(MD, key=str.isdigit)]

        # Get operations
        for i in range(0,len(counts)):
            curr = counts[i]
            try:
                counts[i] = int(curr)
                operations.append("M")
            except ValueError:
                # Handle deletion
                if curr.startswith("^"):
                    operations.append("D")
                    counts[i] = len(counts[i]) - 1
                else:
                    operations.append("X")
                    counts[i] = len(counts[i])

        return operations, counts

def split_cigar(cigar):
    """ Takes CIGAR string from SAM and splits it into two lists:
        one with capital letters (match operators), and one with
        the number of bases that each operation applies to. """

    alignTypes = re.sub('[0-9]', " ", cigar).split()
    counts = re.sub('[=A-Z]', " ", cigar).split()
    counts = [int(i) for i in counts]

    return alignTypes, counts

def compute_transcript_end(start, cigar):
    """ Given the start position and CIGAR string of a mapped SAM transcript,
        compute the end position in the reference genome.
        Args:
            start: The start position of the transcript with respect to the
            forward strand

            cigar: SAM CIGAR string describing match operations to the reference
            genome

        Returns:
            end position of the transcript.
    """
    end = start

    ops, counts = split_cigar(cigar)
    for op,ct in zip(ops, counts):
        if op in ["=", "H", "M", "N", "D"]:
            end += ct

    return end - 1

def compute_jI(start, cigar):
    """ If the input sam file doesn't have the custom STARlong-derived jI tag,
        we need to compute it. This is done by stepping through the CIGAR 
        string, where introns are represented by the N operation.

        start: The start position of the transcript with respect to the
               forward strand
        cigar: SAM CIGAR string describing match operations to the reference
               genome
        Returns: jI string representation of intron start and end positions.

        Example jI strings:
            no introns: jI:B:i,-1
            two introns: jI:B:i,167936516,167951806,167951862,167966628
    """

    operations, counts = split_cigar(cigar)
    jI = ["jI:B:i"]
    genomePos = start

    # Iterate over cigar operations
    for op,ct in zip(operations, counts):
        if op == "N":
            # This is an intron
            intronStart = genomePos
            intronEnd = genomePos + ct - 1

            jI.append(str(intronStart))
            jI.append(str(intronEnd))

        if op not in ["S", "I"]:
            genomePos += ct

    # If the transcript has no introns, add -1 to the tag
    if len(jI) == 1:
        jI.append("-1")

    jIstr = ",".join(jI)
    return jIstr


def get_introns(sam_record: pysam.AlignedSegment, start, cigar):
    """ Locates the jI field in a list of SAM fields or computes
        it from the CIGAR string and start position if it isn't found.
        Note that positions refer to start and endpoints of introns, not exons,
        so adjustments are needed to avoid an off-by-one error if you want exons.

        Example jI strings:
            no introns: jI:B:i,-1
            two introns: jI:B:i,167936516,167951806,167951862,167966628
        Args:
            sam_record: a pysam AlignedSegment
            start: The start position of the transcript with respect to the
            forward strand
            cigar: SAM CIGAR string describing match operations to the reference
            genome
        Returns:
            intron_list: intron starts and ends in a list (sorted order)
    """
    try:
        intron_list = sam_record.get_tag("jI").tolist()
    except KeyError:
        jI = compute_jI(start, cigar)
        intron_list = [int(x) for x in jI.split(",")[1:]]

    if intron_list[0] == -1:
        return []
    else:
        return intron_list
