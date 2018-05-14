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

            transcripts: A dictionary that contains transcript IDs mapped to 
            Transcript objects. These objects are the transcripts that come 
            from this gene.

    """

    def __init__(self, identifier, name, chromosome, start, end, strand):
        start = int(start)
        end = int(end)

        self.name = name
        self.identifier = identifier
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.transcripts = {}

        if start > end:
            raise ValueError('Gene start must be less than or equal to end.')

    def set_name(self, name):
        """ Sets the name attribute of the Gene to the provided value.
        """
        self.name = name
        return

    def add_transcript(self, transcript):
        """ Adds a key-value pair (transcript identifier -> Transcript oject)
            to the gene's transcript dictionary

            Args:
                transcript: object of type Transcript. Must overlap with the 
                location of the gene.
        """
        if transcript.start >= self.end or transcript.end <= self.start:
            raise ValueError('Transcript must overlap the gene it is assigned to')
 
        if transcript.gene_id == self.identifier:
            # In order to belong to a gene, the transcript gene_id must
            # match
            transcript_id = transcript.identifier
            self.transcripts[transcript_id] = transcript
        else:
            raise ValueError('Gene ID of transcript must match gene ' + \
                  'in order for assignment to be made.')
        return             


    def print_gene(self):
        """ Print a string representation of the Gene. Good for debugging. """

        if self.name != "":
            # Include name in output if there is one
            print self.identifier + " (" + self.name + "):"
        else:
            print self.identifier + ":"

        print "\tLocation: " + self.chromosome + ":" + str(self.start) + "-" + \
              str(self.end) + "(" + self.strand + ")"
        
        # Print transcripts in shorthand 
        for transcript in self.transcripts:
            print "\t Transcript: " + transcript

        return

def get_gene_from_db(gene_row):
    """ Uses information from a database gene entry to create a
    Gene object.

        Args:
            gene_row: Tuple-formatted row from 'genes' table of a
            TALON database
    """
    gene_id = gene_row['identifier']
    name = gene_row['name']
    chromosome = gene_row['chromosome']
    start = gene_row['start']
    end = gene_row['end']
    strand = gene_row['strand']

    #transcripts = {} #gene_row['transcript_ids'].split(",")

    gene = Gene(gene_id, name, chromosome, start, end, strand)
    return gene

def get_gene_from_gtf(gene_info):
    """ Creates a Gene object from a GTF file entry
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

    gene = Gene(gene_id, gene_name, chromosome, start, end, strand)
    return gene
