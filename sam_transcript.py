# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from transcript import *
import re

class SamTranscript(Transcript):
    """ Stores information about a gene transcript that comes from a sam file, 
        including its location and constitutive exons. Inherits from the 
        Transcript class.
  
        Attributes:
           sam_id: Unique identifier for the transcript from the SAM file

           chromosome: Chromosome that the transcript is located on
           (format "chr1")

           start: The start position of the transcript with respect to the
           forward strand

           end: The end position of the transcript with respect to the
           forward strand

           strand: "+" if the transcript is on the forward strand, and "-" if
           it is on the reverse strand

           exons: Data structure containing at least one exon

       Optional Attributes:
           gene:

           transcript_id: Accession ID of transcript, i.e. and Ensembl ID

           transcript_name: Human-readable name of the transcript

    """

    def __init__(self, sam_id, identifier, name, chromosome, start, end, \
                 strand, samFields):
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.exons = []
        self.samFields = samFields

        self.sam_id = sam_id
        self.identifier = identifier
        self.name = name

def get_sam_transcript(samFields):
    """ Creates a SamTranscript object from a SAM entry.
        
        Args:
            samFields: List containing fields from a sam entry.

        Returns:
            A SamTranscript object
    """
    
    sam_id = samFields[0]
    flag = samFields[1]
    chromosome = samFields[2]
    start = int(samFields[3])
    cigar = samFields[5]
    seq = samFields[9]

    end = compute_transcript_end(start, cigar)
    print end
    exit()

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
        if op in ["M", "N", "D"]:
            end += ct

    return end - 1

def split_cigar(cigar):
        """ Takes CIGAR string from SAM and splits it into two lists:
            one with capital letters (match operators), and one with
            the number of bases that each operation applies to. """

        alignTypes = re.sub('[0-9]', " ", cigar).split()
        counts = re.sub('[A-Z]', " ", cigar).split()
        counts = [int(i) for i in counts]

        return alignTypes, counts

