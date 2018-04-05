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
 
            gene_ids: A dictionary containing the ID of every gene in the
            structure. Used to detect collisions.
    """
    def __init__(self):
        self.chromosomes = {}
        self.gene_ids = {}

    def add_chromosome(self, chrom_name):
        """ Add a chromosome (with empty interval tree) to the GeneTree

            Args:
                chrom_name: Name of the chromosome. This name will be used as 
                a key to access the interval tree belonging to this chromosome.
        """
        self.chromosomes[chrom_name] = IntervalTree()
        return

    def add_gene(self, gene_id, gene_name, chromosome, start, end, strand):
        """ Creates a gene object and adds it to the GeneTree. The gene's 
            start-end interval is added to the chromosome's interval tree, and 
            is used as a key to retrieve the gene. 

            All positions are 1-based.

            Args:
                gene_id: Accession ID of the gene. Must be unique.

                gene_name: Human-readable name of the gene. Optional, and does
                not have to be unique (although that is generally preferable)

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

        if gene_id in self.gene_ids:
            raise KeyError('Gene IDs must be unique. ' + gene_id + \
                           " is duplicated.")       
 
        gene = Gene(gene_id, chromosome, start, end, strand)

        # Set gene name if available
        if gene_name != None:
            gene.set_name(gene_name)

        self.chromosomes[chromosome][start:end] = gene
        self.gene_ids[gene_id] = 1

        return  

    def add_gene_from_gtf(self, gene_info):
        """ Adds gene from a GTF file to the GeneTree

            Args:
                gene_info: A list containing fields from a GTF file gene entry.
                Example:
                ['chr1', 'HAVANA', 'gene', '11869', '14409', '.', '+', '.',
                'gene_id "ENSG00000223972.5";
                gene_type "transcribed_unprocessed_pseudogene";
                gene_status "KNOWN"; gene_name "DDX11L1"; level 2;
                havana_gene "OTTHUMG00000000961.2";']

        """
        gene_name = None
        chromosome = gene_info[0]
        description = gene_info[-1]
        if "gene_id" not in description:
            raise ValueError('GTF entry lacks a gene_id field')
        gene_id = (description.split("gene_id ")[1]).split('"')[1]

        if "gene_name" in description:
            gene_name = (description.split("gene_name ")[1]).split('"')[1]
        start = int(gene_info[3])
        end = int(gene_info[4])
        strand = gene_info[6]

        self.add_gene(gene_id, gene_name, chromosome, start, end, strand)
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
            gene = interval.data
            if gene.strand == strand:
                overlapping_genes.append(gene)

        return overlapping_genes

    def print_tree(self):
        """ Print a rudimentary visual representation of the GeneTree. """
        # TODO: It would be nice if it printed the chromosomes in order 
        for chrom in self.chromosomes:
            print chrom + ":"
            for gene_interval in self.chromosomes[chrom]:
                print "\t" + str(gene_interval)
        return               


