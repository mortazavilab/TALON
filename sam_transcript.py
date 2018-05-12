# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from exon import *
from genetree import *
import re
from transcript import *


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

    def __init__(self, sam_id, chromosome, start, end, strand, introns,
                 samFields):
        self.sam_id = sam_id
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.introns = introns
        self.samFields = samFields

        self.identifier = None
        self.gene_id = None
        self.name = None
        self.exons = self.create_sam_exons()

    def create_sam_exons(self):
        """ Uses intron coordinates and transcript start/end to create exon
            objects for the transcript """
        introns = self.introns
        starts = [self.start]
        ends = []

        # First, compute exon coordinates and put starts and ends in lists
        for i in range(0,len(introns)):
            if i % 2 == 0:
                # This is an intron start, i.e. an exon end. Subtract 1 to get 
                # the exon base
                starts.append(introns[i] - 1)
            else:
                # This is an intron end, i.e. an exon start. Add 1 to get the
                # exon base
                ends.append(introns[i] + 1)
        ends.append(self.end)

        # Now iterate over start and end pairs and create exon objects
        exons = []
        ct = 1
        for s,e in zip(starts,ends):
            exon_id = self.sam_id + "_" + str(ct)
            exon = Exon(exon_id, self.chromosome, s, e, self.strand, 
                        self.gene_id, self.identifier)
            exons.append(exon)
            ct += 1
        return exons
      

def get_sam_transcript(samFields):
    """ Creates a SamTranscript object from a SAM entry.
        
        Args:
            samFields: List containing fields from a sam entry.

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
                        samFields)
    return sam


#def get_exons(start, end, introns):
#    """ Transforms intron coordinates and adds in start and end to create
#        a list representing the exon starts and ends in the transcript."""
    
#    exons = [start]
#    for i in range(0,len(introns)):
#        if i % 2 == 0:
            # This is an intron start, i.e. an exon end. Subtract 1 to get the
            # exon base
#            exons.append(introns[i] - 1)
#        else:
            # This is an intron end, i.e. an exon start. Add 1 to get the
            # exon base
#            exons.append(introns[i] + 1)
#    exons.append(end)
#    return exons
    

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

def get_loose_exon_matches(chromosome, start, end, strand, exons, cutoff_5, \
                           cutoff_3):
    """ Finds and returns all annotation exons that overlap the query location 
        and which start and end within 10 basepairs of the query.

        Args:
            chromosome: Chromosome that the query is located on
            (format "chr1")

            start: The start position of the query with respect to the
            forward strand

            end: The end position of the query with respect to the
            forward strand

            strand: "+" if the query is on the forward strand, and "-" if
            it is on the reverse strand

            exons: An ExonTree object with genes organized by chromosome in an
            interval tree data structure

            cutoff_5: Permitted difference at 5' end
     
            cutoff_3: Permitted difference at 3' end
    """
    loose_matches = []
    differences = []

    # Get all exons that overlap the query location
    exon_matches = exons.get_exons_in_range(chromosome, start, end, strand)
    
    # Compute the 5' and 3' differences
    for match in exon_matches:
        query_range = [start, end]
        match_range = [match.start, match.end]
         
        diff_5, diff_3 = get_difference(query_range, match_range, strand)
        
        # Enforce similarity cutoff
        if abs(diff_5) <= cutoff_5 and abs(diff_3) <= cutoff_3:
            loose_matches.append(match)
            differences.append([diff_5, diff_3])
    return loose_matches, differences

def get_difference(a, b, strand):
    """ Computes the 5' and 3' difference between two exon intervals. 

        Example: a = [ 0, 10]  b = [ 2, 8 ] on + strand
            5' difference =  -2
            3' difference =  +2

        Args:
            a: First interval, formattted as a list
            b: Second interval, formatted as a list
    """
    if strand == "+":
        diff_5 = a[0] - b[0] 
        diff_3 = a[1] - b[1]

    elif strand == "-":    
        diff_5 = b[1] - a[1]
        diff_3 = b[0] - a[0]

    return diff_5, diff_3

def get_overlap(a, b):
    """ Computes the amount of overlap between two intervals.
        Returns 0 if there is no overlap. The function treats the start and 
        ends of each interval as inclusive, meaning that if a = b = [10, 20], 
        the overlap reported would be 11, not 10.

        Args: 
            a: First interval, formattted as a list
            b: Second interval, formatted as a list
    """
    overlap = max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
    return overlap


