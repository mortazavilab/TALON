# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program reads in SAM-formatted long read alignments and adds a custom 
# tag to reflect the fraction of As in the sequence immediately following the
# alignment. This can help indicate the likelihood of an internal priming 
# artifact. 

import pyfaidx
import pysam
from optparse import OptionParser

def get_options():
    """ Read input args """

    parser = OptionParser(description=(""))
    parser.add_option("--f", dest = "sam_file",
                      help = "SAM file of transcripts")
    parser.add_option("--g", dest = "genome_file",
                      help = "Reference genome fasta file")
    parser.add_option("--ar", dest = "fracA_range_size", type = int,
                      help = ("Size of post-transcript interval to compute "
                              "fraction As on. Default = 10"), default = 10)
    
    parser.add_option("--o", dest = "outprefix",
                      help = "Prefix for outfiles")

    (opts, args) = parser.parse_args()
    return opts

def fetch_seq(chrom=str, start=int, stop=int, strand=str, genome=pyfaidx.Fasta,
              indexing = 0):
    """ Given a genomic interval, return the sequence with respect to the
        strand supplied.
        If 1-based indexing is specified, then 1 will be subtracted from the
        position to convert to the Python indexing. """

    if start > stop:
        raise ValueError("Start must be less than or equal to stop")

    if indexing != 0:
        if indexing == 1:
            start -= 1
        else:
            raise ValueError("Valid indexing modes include: 1 or 0")

    seq = genome[chrom][start:stop]

    if strand == "-":
        seq = seq.reverse.complement

    return str(seq)

def compute_frac_As(seq=str):
    """ Compute fraction of sequence made up of As """

    a = seq.count('A')
    n = len(seq)
    if n == 0:
        return 0
    else:
        return float(a)/n

def fetch_range_after_transcript(transcript_end=int, strand=str, length=int):
    """ Given the 1-based stop position of a transcript and its strand,
        return a 1-based genomic range of the specified length that starts with
        the base just after the end position. The smaller position is always
        reported first.
        Example:
              fetch_range_after_transcript(4, '+', 2) would yield (5, 6)
              fetch_range_after_transcript(4, '-', 2) would yield (2, 3)
    """
    if length < 1:
        raise ValueError("Length must be greater than or equal to 1")

    if strand == '+':
        range_start = transcript_end + 1
        range_end = range_start + length - 1
    elif strand == '-':
        range_start = transcript_end - 1
        range_end = range_start - length + 1
    else:
        raise ValueError("Strand must be + or -")

    return (min(range_start, range_end), max(range_start, range_end))

def compute_transcript_end(transcript=pysam.AlignedSegment):
    """ Compute the position of the final transcript base relative to the genome,
        taking strand into account. Position is 1-based. """

    strand = "-" if transcript.is_reverse else "+"
    if strand == '+':
        return transcript.reference_end
    if strand == '-':
        return transcript.reference_start + 1 # (make 1-based)

def compute_frac_as_after_transcript(chrom=str, transcript_end=int, strand=str,
                                     range_size=int, genome=pyfaidx.Fasta):
    """ Given a transcript end, strand, range size, and genome object,
        compute the fraction of sequence in the range immediately after
        the transcript end that is made up of As."""

    # Get sequence of range immediately after transcript
    range_start, range_end = fetch_range_after_transcript(transcript_end,
                                                          strand, range_size)
    range_seq = fetch_seq(chrom, range_start, range_end, strand, genome,
                          indexing = 1)

    # Get fraction As in sequence
    return compute_frac_As(range_seq)

def main(options=None):
    if options == None:
        options = get_options()

    genome = pyfaidx.Fasta(options.genome_file, sequence_always_upper=True,
                           one_based_attributes=False)

    frac_A_outfile = options.outprefix + "_fraction_As_%dbp_after_transcript.tsv" \
                     % (options.fracA_range_size)
    o_afrac = open(frac_A_outfile, 'w')
    o_afrac.write("\t".join(["read_name", "fraction_As"]) + '\n')
    with pysam.AlignmentFile(options.sam_file) as sam:
        for record in sam:  # type: pysam.AlignedSegment
            if record.is_secondary == True or record.is_unmapped == True:
                continue
            read_id = record.query_name
            chrom = record.reference_name
            strand = "-" if record.is_reverse else "+"
            transcript_end = compute_transcript_end(record)
            frac_As = compute_frac_as_after_transcript(chrom, transcript_end, strand,
                                                       options.fracA_range_size,
                                                       genome)
            o_afrac.write("\t".join([read_id, str(frac_As)]) + '\n')

    o_afrac.close()

if __name__ == '__main__':
    options = get_options()
    main(options) 
