# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

class Gene(object):
    """ Contains high-level information about a gene, such as its identifiers, 
        genomic location, and transcripts. Does not contain exon information.

        Attributes:

            identifier: Accession ID of gene, i.e. an Ensembl ID. Required.

            name: Human-readable name of the gene. This attribute can be left 
            empty if the gene does not have an assigned name.

            chromosome: Chromosome that the gene is located on (format "chr1")

            start: The start position of the gene with respect to the forward 
            strand (int). Should always be less than or equal to end.

            end: The end position of the gene with respect to the forward strand 
            (int). Should always be greater than or equal to start.

            strand: "+" if the gene is on the forward strand, "-" if it is on 
            the reverse strand

            transcripts:

    """
    # TODO: Update transcript description in comment

    def __init__(self, identifier, chromosome, start, end, strand):
        self.name = ""
        self.identifier = identifier
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.transcripts = {}

    def set_name(self, name):
        """ Sets the name attribute of the Gene to the provided value.
        """
        self.name = name
        return

    def print_gene(self):
        """ Print a string representation of the Gene. Good for debugging. """

        if self.name != "":
            # Include identifier in output if there is one
            print self.identifier + " (" + self.name + "):"
        else:
            print self.identifier + ":"

        print "\tLocation: " + self.chromosome + ":" + str(self.start) + "-" + \
              str(self.end) + "(" + self.strand + ")"
        
        # TODO: Print transcripts too, at least in shorthand 
        return
