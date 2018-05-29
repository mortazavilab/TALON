# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from intervaltree import Interval, IntervalTree
from exon import *

class ExonTree(object):
    """ Stores locations of exons as intervals with the goal of querying them 
        for overlap.

        The structure is a dictionary of chromosome names, each mapped to an
        interval tree. Each interval tree is made up of intervals that map to 
        corresponding exon ids.

        Attributes:
            chromosomes: A dictionary mapping chromosome name to an interval 
            tree containing exons that are located on that chromosome.

            exons: A dictionary mapping exon accession IDs to actual exon 
            objects

            novel: A counter keeping track of the number of novel exons added.
            Used to generate next novel ID
 
    """
    def __init__(self):
        self.chromosomes = {}
        self.exons = {}
        self.novel = 0
        

    def add_chromosome(self, chrom_name):
        """ Add a chromosome (with empty interval tree) to the GeneTree

            Args:
                chrom_name: Name of the chromosome. This name will be used as 
                a key to access the interval tree belonging to this chromosome.
        """
        self.chromosomes[chrom_name] = IntervalTree()
        return

    def add_exon(self, exon):
        """ Adds an Exon object to the ExonTree. The exon's 
            start-end interval is added to the chromosome's interval tree, and 
            is used as a key to retrieve the exon. 

            All positions are 1-based.

            Args:
                exon: Exon object to be added
                exon_id: Accession ID to assign to the exon
                
        """
        chromosome = exon.chromosome
        start = exon.start
        end = exon.end
        exon_id = exon.identifier

        if chromosome not in self.chromosomes:
            self.add_chromosome(chromosome)
        if exon.start == exon.end:
            return

        # Check for collisions. There are two ways a collision can happen.
        # By far the most common case is an exon that is in multiple 
        # transcripts. For this case, merge the exon transcript sets. A more
        # rare special case is when two exons have the same ID but belong to 
        # different genes
        if exon_id in self.exons:
             oldExon = self.exons[exon_id]
             exon.transcript_ids = oldExon.transcript_ids | exon.transcript_ids
        else: 
            self.chromosomes[chromosome][start:end] = exon_id
    
        self.exons[exon_id] = exon
        return

    def add_novel_exon(self, chromosome, start, end, strand, gene_id, transcript_id):
        """ Creates an exon from the provided information and adds it to the
            exon tree. It is assigned an ID based on the number of novel exons
            that are in the tree already.

            Args:
                chromosome: Name of the chromosome of the exon

                start: The start position of the exon with respect to the
                forward strand (int). Should always be less than or equal to
                end.

                end: The end position of the exon with respect to the
                forward strand (int). Should always be greater than or equal
                to start.

                strand: "+" if the query interval is on the forward strand,
                "-" if it is on the reverse strand

        """
        self.novel += 1
        new_id = "novel_exon." + str(self.novel)
        new_exon = Exon(new_id, chromosome, start, end, strand, gene_id, transcript_id)
        self.add_exon(new_exon, new_id)
        return

    def get_exons_in_range(self, chromosome, start, end, strand):
        """ Finds exons that overlap with the provided start-end interval.
            
            Args:
                chromosome: Name of the chromosome of the query interval
                
                start: The start position of the interval with respect to the
                forward strand (int). Should always be less than or equal to
                end.

                end: The end position of the interval with respect to the 
                forward strand (int). Should always be greater than or equal 
                to start.

                strand: "+" if the query interval is on the forward strand, 
                "-" if it is on the reverse strand

            Returns:
                List of exons (objects of class Exon) that overlap the query 
                interval. 
        """
        if start > end:
            raise ValueError('Exon start must be less than or equal to end.')

        if chromosome not in self.chromosomes:
            raise ValueError("Query chromosome not found: " + chromosome)

        # Add 1 to the query end because in the interval tree data structure,
        # ranges are inclusive of the lower limit, but non-inclusive of the 
        # upper limit
        end += 1
 
        overlapping_exon_intervals  = self.chromosomes[chromosome][start:end]

        # Only report exons on the same strand as the query
        overlapping_exons = []
        for interval in overlapping_exon_intervals:
            exon_id = interval.data
            exon = self.exons[exon_id]
            if exon.strand == strand:
                overlapping_exons.append(exon)

        return overlapping_exons

    def print_tree(self):
        """ Print a rudimentary visual representation of the ExonTree. """
        # TODO: It would be nice if it printed the chromosomes in order 
        for chrom in self.chromosomes:
            print chrom + ":"
            for exon_interval in self.chromosomes[chrom]:
                print "\t" + str(exon_interval)
        return               


