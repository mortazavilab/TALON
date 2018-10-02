# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from intervaltree import Interval, IntervalTree
from gene import Gene

class GeneTree(object):
    """ Stores locations of genes as intervals with the goal of querying them 
        for overlap.
        The structure is a dictionary of chromosome names, each mapped to an
        interval tree. Each interval tree is made up of intervals that map to 
        corresponding gene names.
        Attributes:
            chromosomes: A dictionary mapping chromosome name to an interval 
            tree containing genes that are located on that chromosome.
 
            gene_ids: A dictionary mapping the ID of every gene to its gene
            object. Also useful for detectingcollisions.
    """
    def __init__(self):
        self.chromosomes = {}
        self.genes = {}

    def add_chromosome(self, chrom_name):
        """ Add a chromosome (with empty interval tree) to the GeneTree
            Args:
                chrom_name: Name of the chromosome. This name will be used as 
                a key to access the interval tree belonging to this chromosome.
        """
        self.chromosomes[chrom_name] = IntervalTree()
        return

    def add_gene(self, gene):
        """ Adds a Gene object to the GeneTree. The gene's 
            start-end interval is added to the chromosome's interval tree, and 
            is used as a key to retrieve the gene. 
            All positions are 1-based.
            Args:
                gene: Gene object to be added
        """
        gene_id = gene.identifier
        chromosome = gene.chromosome
        start = gene.start
        end = gene.end

        if chromosome not in self.chromosomes:
            self.add_chromosome(chromosome)

        if gene_id in self.genes:
            raise KeyError('Gene IDs must be unique. ' + gene_id + \
                           " is duplicated.")       
 
        self.chromosomes[chromosome][start:end] = gene_id
        self.genes[gene_id] = gene

        return

    def get_genes_in_range(self, chromosome, start, end, strand):
        """ Finds genes that overlap with the provided start-end interval.
            
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
                List of genes (objects of class Gene) that overlap the query 
                interval. 
        """
        if start > end:
            raise ValueError('Gene start must be less than or equal to end.')

        if chromosome not in self.chromosomes:
            raise KeyError('Chromosome not found among genes: ' + chromosome)

        # Add 1 to the query end because in the interval tree data structure,
        # ranges are inclusive of the lower limit, but non-inclusive of the 
        # upper limit
        end += 1
 
        overlapping_gene_intervals  = self.chromosomes[chromosome][start:end]

        # Only report genes on the same strand as the query
        overlapping_genes = []
        for interval in overlapping_gene_intervals:
            gene_id = interval.data
            gene = self.genes[gene_id]
            if gene.strand == strand:
                overlapping_genes.append(gene)

        return overlapping_genes

    def print_tree(self):
        """ Print a rudimentary visual representation of the GeneTree. """
        for chrom in self.chromosomes:
            print chrom + ":"
            for gene_interval in self.chromosomes[chrom]:
                print "\t" + str(gene_interval)
        return 
