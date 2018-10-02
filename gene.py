# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

class Gene(object):
    """ Contains high-level information about a gene, such as its identifiers, 
        genomic location, and transcripts. Does not contain exon information.
        Attributes:
            - identifier: Accession ID of gene, i.e. an Ensembl ID. Required.
            - name: Human-readable name of the gene. This attribute can be left 
              empty if the gene does not have an assigned name.
            - chromosome: Chromosome that the gene is located on (format "chr1")
            - start: The start position of the gene with respect to the forward 
              strand (int). Should always be less than or equal to end.
            - end: The end position of the gene with respect to the forward strand 
              (int). Should always be greater than or equal to start.
            - strand: "+" if the gene is on the forward strand, "-" if it is on 
              the reverse strand
            - annotations: a dictionary of miscellaneous annotation categories
              extracted from a GTF
            
    """

    def __init__(self, identifier, chromosome, start, end, strand, annotations):
        start = int(start)
        end = int(end)

        self.identifier = str(identifier)
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.transcripts = {}
        self.length = end - start + 1
        self.annotations = annotations

        if start > end:
            raise ValueError("""Plus strand gene start must be less than or 
                             equal to end.""")

    def set_name(self, name):
        """ Sets the name attribute of the Gene to the provided value.
        """
        self.annotations['name'] = name
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
            # In order to belong to a gene, the transcript gene_id must match
            transcript_id = transcript.identifier
            self.transcripts[transcript_id] = transcript
        else:
            raise ValueError('Gene ID of transcript must match gene ' + \
                  'in order for assignment to be made.')
        return             


    def print_gene(self):
        """ Print a string representation of the Gene. Good for debugging. """

        if "name" in self.annotations != "":
            # Include name in output if there is one
            print self.identifier + " (" + self.annotations['name']  + "):"
        else:
            print self.identifier + ":"

        print "\tLocation: " + self.chromosome + ":" + str(self.start) + "-" + \
              str(self.end) + "(" + self.strand + ")"
        
        # Print transcripts in shorthand 
        for transcript in self.transcripts:
            print "\t Transcript: " + transcript

        return

def get_gene_from_db(gene_start_row, gene_end_row):
    """ Uses information from a database gene entry to create a
    Gene object.
        Args:
            gene_row: Tuple-formatted row from 'genes' table of a
            TALON database
    """
    if gene_start_row['gene_id'] != gene_end_row['gene_id']:
            raise ValueError("get_gene_from_db: provided start and stop " + \
                             "come from different genes")
    gene_id = gene_start_row['gene_ID']
    chromosome = gene_start_row['chromosome']
    start = gene_start_row[2]
    end = gene_end_row[2]
    strand = gene_start_row['strand']

    #transcripts = {} #gene_row['transcript_ids'].split(",")

    gene = Gene(gene_id, chromosome, start, end, strand, {})
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
    chromosome = gene_info[0]
    start = int(gene_info[3])
    end = int(gene_info[4])
    strand = gene_info[6]
    annotations = extract_gene_annotations_from_GTF(gene_info)
    if "gene_id" not in gene_info[-1]:
            raise ValueError('GTF entry lacks a gene_id field')
    gene_id = annotations['gene_id']

    gene = Gene(gene_id, chromosome, start, end, strand, annotations)
    return gene

def extract_gene_annotations_from_GTF(tab_fields):
    """Parses the description field of a gene GTF in order to organize the 
       information therein into a dictionary.
    """
    attributes = {}

    description = tab_fields[-1].strip()
    # Parse description
    for pair in [x.strip() for x in description.split(";")]:
        if pair == "": continue

        pair = pair.replace('"', '')
        key, val = pair.split()
        attributes[key] = val

    attributes["source"] = tab_fields[1]
    return attributes  


def get_gene_from_exon(exon, gene_id):
    """ In rare cases, GTF exons are listed with gene and transcript IDs that
        do not have corresponding entries. In this case, we create a gene
        for this exon for bookkeeping purposes."""

    gene_name = gene_id
    chromosome = exon.chromosome
    start = exon.start
    end = exon.end
    strand = exon.strand
    gene = Gene(gene_id, gene_name, None, chromosome, start, end, strand)
    return gene

def create_novel_gene(chromosome, start, end, strand, counter):
    """ Creates a novel gene with a unique identifier (obtained using
        counter). Returns the gene object as well as the updated counter.
    """
    gene_id = str(counter["genes"] + 1)
    counter["genes"] += 1
    gene = Gene(gene_id, chromosome, start, end, strand, None)
    return gene
