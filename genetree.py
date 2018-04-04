# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from intervaltree import Interval, IntervalTree

class GeneTree(object):
    """ Stores locations of genes as intervals with the goal of querying them 
        for overlap.

        The structure is a dictionary of chromosome names, each mapped to an
        interval tree. Each interval tree is made up of intervals that map to 
        corresponding gene names.

        Attributes:
            chromosomes: A dictionary mapping chromosome name to an interval 
            tree containing genes that are located on that chromosome.
    """
    def __init__(self):
        self.chromosomes = {}

    def add_chromosome(self, chrom_name):
        """ Add a chromosome (with empty interval tree) to the GeneTree

            Args:
                chrom_name: Name of the chromosome. This name will be used as 
                a key to access the interval tree belonging to this chromosome.
        """
        self.chromosomes[chrom_name] = IntervalTree()
        return

    def add_gene(self, gene, chromosome, start, end, strand):
        # TODO: gene object instead of gene name
        """ Add a gene to the GeneTree. The gene's start-end interval is added 
            to the chromosome's interval tree, and is used as a key to retrieve 
            the gene name. 

            Args:
                gene: Name of the gene (string)

                chromosome: Name of the chromosome that the gene is on 
                (string).

                start: The start position of the gene with respect to the 
                forward strand (int). Should always be less than or equal to 
                end.
            
                end: The end position of the gene with respect to the forward
                strand (int). Should always be greater than or equal to start.

                strand: "+" if the gene is on the forward strand, "-" if it is 
                on the reverse strand
        """
        if start > end:
            raise ValueError('Gene start must be less than or equal to end.')

        if chromosome not in self.chromosomes:
            self.add_chromosome(chromosome)
        
        self.chromosomes[chromosome][start:end] = gene
        return  

    def print_tree(self):
       """ Print a rudimentary visual representation of the GeneTree. """
       # TODO: It would be nice if it printed the chromosomes in order 
       for chrom in self.chromosomes:
           print chrom + ":"
           for gene_interval in self.chromosomes[chrom]:
               print "\t" + str(gene_interval)
       return               


