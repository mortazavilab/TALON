# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Useful functions for processing transcripts in the SAM format

import itertools
import re

def compute_alignment_coverage(CIGAR):
    """ This function computes what fraction of the read is actually aligned to
        the genome by excluding hard or soft-clipped bases."""

    total_bases = 0.0
    unaligned_bases = 0.0
    ops, counts = split_cigar(CIGAR)
    for op,ct in zip(ops, counts):
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

def splitMD(MD_tag):
        """ Takes MD tag and splits into two lists:
            one with capital letters (match operators), and one with
            the number of bases that each operation applies to. """

        MD = str(MD_tag).split(":")[2]
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
    counts = re.sub('[A-Z]', " ", cigar).split()
    counts = [int(i) for i in counts]

    return alignTypes, counts
