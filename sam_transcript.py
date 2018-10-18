# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

import edge as Edge
import genetree as GeneTree
import re
from transcript import *

class SamTranscript(Transcript):
    """ Stores information about a gene transcript that comes from a sam file, 
        including its location and constitutive exons. Inherits from the 
        Transcript class.
  
        Attributes:
           identifier: Unique identifier for the transcript from the SAM file
           chromosome: Chromosome that the transcript is located on
           (format "chr1")
           start: The start position of the transcript with respect to the
           forward strand
           end: The end position of the transcript with respect to the
           forward strand
           strand: "+" if the transcript is on the forward strand, and "-" if
           it is on the reverse strand
           exons: Data structure containing at least one exon
    """

    def __init__(self, sam_id, chromosome, start, end, strand, introns,
                 samFields, dataset):
        self.identifier = str(sam_id)
        self.chromosome = str(chromosome)
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.intron_coords = introns

        self.samFields = samFields
        self.dataset = dataset

        self.name = None
        self.exons = []
        self.introns = []
        self.n_exons = 0
 
        self.create_sam_exons()
        self.create_sam_introns()

        # Fields that will be set during TALON run
        self.gene_ID = None
        self.transcript_ID = None
        self.diff_5 = None
        self.diff_3 = None
        novel = None

    def create_sam_exons(self):
        """ Uses intron coordinates and transcript start/end to create exon
            objects for the transcript """
        introns = self.intron_coords
        starts = [self.start]
        ends = []
     
        # First, compute exon coordinates and put starts and ends in lists
        for i in range(0,len(introns)):
            if i % 2 == 0:
                # This is an intron start, i.e. an exon end. Subtract 1 to get 
                # the exon base
                ends.append(introns[i] - 1)
            else:
                # This is an intron end, i.e. an exon start. Add 1 to get the
                # exon base
                starts.append(introns[i] + 1)
        ends.append(self.end)

        # Now iterate over start and end pairs and create exon objects
        ct = 1
        for curr_start,curr_end in zip(starts,ends):
            exon_id = None
            if curr_start == curr_end:
                raise ValueError("Transcript " + self.identifier + \
                                 "has exon of length zero.")
            exon = Edge.Edge(exon_id, self.chromosome, curr_start, curr_end, self.strand, 
                        None, None, {})
            self.add_exon(exon)
            ct += 1
        return
     
    def create_sam_introns(self):
        introns = []
        i = 0
        while i < len(self.intron_coords):
            s = self.intron_coords[i] - 1 # adjust to vertex-based system
            e = self.intron_coords[i+1] + 1 # adjust to vertex-based system
            intron = Edge.Edge(None, self.chromosome, s, e, self.strand,
                        None, None, {})
            introns.append(intron)
            i += 2
        self.introns = introns
        return
             

def get_sam_transcript(samFields, dataset):
    """ Creates a SamTranscript object from a SAM entry.
        
        Args:
            samFields: List containing fields from a sam entry.
            dataset: Name of the dataset that the sam entry comes from
        Returns:
            A SamTranscript object
    """
    sam_id = samFields[0]
    flag = int(samFields[1])
    chromosome = samFields[2]
    start = int(samFields[3])

    cigar = samFields[5]
    seq = samFields[9]
    otherFields = samFields[11:len(samFields)]

    end = compute_transcript_end(start, cigar)
    introns = get_introns(otherFields, start, cigar)

    if flag in [16, 272]:
        strand = "-"
    else:
        strand = "+" 
    sam = SamTranscript(sam_id, chromosome, start, end, strand, introns, 
                        samFields, dataset)
    return sam


def get_introns(fields, start, cigar):
    """ Locates the jI field in a list of SAM fields or computes
        it from the CIGAR string and start position if it isn't found. 
   
        Example jI strings:
            no introns: jI:B:i,-1
            two introns: jI:B:i,167936516,167951806,167951862,167966628
        Args:
            fields: List containing fields from a sam entry.
            start: The start position of the transcript with respect to the
            forward strand
            cigar: SAM CIGAR string describing match operations to the reference
            genome
        Returns:
            intron_list: intron starts and ends in a list (sorted order)
    """
    indices = [i for i, s in enumerate(fields) if 'jI:B:i' in s]

    if len(indices) == 1:
        jI = fields[indices[0]]
    elif len(indices) == 0:
        jI = compute_jI(start, cigar)
    else:
        raise ValueError('SAM entry contains more than one jI:B:i field')

    intron_list = [int(x) for x in jI.split(",")[1:]]
    if intron_list[0] == -1:
        return []
    else:
        return intron_list

def compute_jI(start, cigar):
    """ If the input sam file doesn't have the custom STARlong-derived jI tag, 
        we need to compute it. This is done by stepping
            through the CIGAR string, where introns are represented by the N
            operation.
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
        if op in ["H", "M", "N", "D"]:
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
