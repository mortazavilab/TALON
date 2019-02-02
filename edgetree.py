# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from intervaltree import Interval, IntervalTree
from edge import *
import pdb

class EdgeTree(object):
    """ Stores locations of edges as intervals with the goal of querying them 
        for overlap.
        The structure is a dictionary of chromosome names, each mapped to an
        interval tree. Each interval tree is made up of intervals that map to 
        corresponding edge ids.
        Attributes:
            chromosomes: A dictionary mapping chromosome name to an interval 
            tree containing edges that are located on that chromosome.
            edges: A dictionary mapping edge accession IDs to actual edge 
            objects
            novel: A counter keeping track of the number of novel edges added.
            Used to generate next novel ID
 
    """
    def __init__(self):
        self.chromosomes = {}
        self.edges = {}
        self.novel = 0
        

    def add_chromosome(self, chrom_name):
        """ Add a chromosome (with empty interval tree) to the GeneTree
            Args:
                chrom_name: Name of the chromosome. This name will be used as 
                a key to access the interval tree belonging to this chromosome.
        """
        self.chromosomes[chrom_name] = IntervalTree()
        return

    def add_edge(self, edge):
        """ Adds an Edge object to the EdgeTree. The edge's 
            start-end interval is added to the chromosome's interval tree, and 
            is used as a key to retrieve the edge. 
            All positions are 1-based.
            Args:
                edge: Edge object to be added
                
        """
        chromosome = edge.chromosome
        start = edge.start
        end = edge.end
        edge_id = str(edge.identifier)

        if chromosome not in self.chromosomes:
            self.add_chromosome(chromosome)
        if edge.start == edge.end:
            print("Ignoring edge with ID " + edge_id + " because its length is zero.")
            return

        # Check for collisions. There are two ways a collision can happen.
        # By far the most common case is an edge that is in multiple 
        # transcripts. For this case, merge the edge transcript sets. 
        # IntervalTrees allow redundancy for intervals with the exact same
        # positions, so it is OK if two exons from different genes occupy 
        # identical intervals
        if edge_id in self.edges:
             oldEdge = self.edges[edge_id]
             edge.transcript_ids = oldEdge.transcript_ids | edge.transcript_ids
        else: 
            self.chromosomes[chromosome][start:end] = edge_id
    
        self.edges[edge_id] = edge
        return

    #def add_novel_edge(self, chromosome, start, end, strand, gene_id, transcript_id):
    #    """ Creates an edge from the provided information and adds it to the
    #        edge tree. It is assigned an ID based on the number of novel edges
    #        that are in the tree already.

    #        Args:
    #            chromosome: Name of the chromosome of the edge

    #            start: The start position of the edge with respect to the
    #            forward strand (int). Should always be less than or equal to
    #            end.

    #            end: The end position of the edge with respect to the
    #            forward strand (int). Should always be greater than or equal
    #            to start.

    #            strand: "+" if the query interval is on the forward strand,
    #            "-" if it is on the reverse strand

    #    """
    #    self.novel += 1
    #    new_id = "novel_edge." + str(self.novel)
    #    new_edge = Edge(new_id, chromosome, start, end, strand, gene_id, transcript_id)
    #    self.add_edge(new_edge, new_id)
    #    return

    def get_edges_in_range(self, chromosome, start, end, strand):
        """ Finds edges that overlap with the provided start-end interval.
            
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
                List of edges (objects of class Edge) that overlap the query 
                interval. 
        """
        if start > end:
            raise ValueError('Edge start must be less than or equal to end.')

        if chromosome not in self.chromosomes:
            return []
            #raise ValueError("Query chromosome not found: " + chromosome)

        # Add 1 to the query end because in the interval tree data structure,
        # ranges are inclusive of the lower limit, but non-inclusive of the 
        # upper limit
        end += 1
 
        overlapping_edge_intervals  = self.chromosomes[chromosome][start:end]

        # Only report edges on the same strand as the query
        overlapping_edges = []
        for interval in overlapping_edge_intervals:
            edge_id = interval.data
            edge = self.edges[edge_id]
            if edge.strand == strand:
                overlapping_edges.append(edge)

        return overlapping_edges

    def print_tree(self):
        """ Print a rudimentary visual representation of the EdgeTree. """
        # TODO: It would be nice if it printed the chromosomes in order 
        for chrom in self.chromosomes:
            print(chrom + ":")
            for edge_interval in self.chromosomes[chrom]:
                print("\t" + str(edge_interval))
        return               

