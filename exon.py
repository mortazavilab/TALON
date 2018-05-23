# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

class Exon(object):
    """Stores information about an exon, including its location
       and the gene/transcript(s) it belongs to.

       Attributes:
           identifier: Accession ID of the exon

           gene: Accession ID of the gene that the exon belongs to

           transcript_ids: Set of transcript accession IDs that the exon 
           belongs to

           chromosome: Chromosome that the transcript is located on 
           (format "chr1")

           start: The start position of the exon with respect to the
           forward strand 

           end: The end position of the exon with respect to the
           forward strand

           strand: "+" if the exon is on the forward strand, and "-" if
           it is on the reverse strand
 
           length: The length of the exon

    """

    def __init__(self, identifier, chromosome, \
                 start, end, strand, gene_id, transcript_id):
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.length = self.end - self.start + 1
        self.name = None

        self.identifier = identifier
        self.gene_id = gene_id
        self.transcript_ids = set()
        if transcript_id != None:
            self.transcript_ids.add(transcript_id)

    def print_exon(self):
        """ Prints a string representation of the exon"""
        print self.identifier + ": " + self.chromosome + ":" + \
              str(self.start) + "-" + str(self.end)
        print self.transcript_ids
        return

def create_exon_from_gtf(exon_info):
    """ Creates an exon object using information from a GTF entry

            Args:
               exon_info: A list containing fields from a GTF file exon entry.
               Example:   
               ['chr1', 'HAVANA', 'exon', '11869', '12227', '.', '+', '.', 
                'gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; 
                gene_type "transcribed_unprocessed_pseudogene"; 
                gene_status "KNOWN"; gene_name "DDX11L1"; 
                transcript_type "processed_transcript"; 
                transcript_status "KNOWN"; transcript_name "DDX11L1-002"; 
                exon_number 1; exon_id "ENSE00002234944.1"; level 2; 
                tag "basic"; transcript_support_level "1"; 
                havana_gene "OTTHUMG00000000961.2"; 
                havana_transcript "OTTHUMT00000362751.1";'] 
    """
    description = exon_info[-1]
    start = int(exon_info[3])
    end = int(exon_info[4])
    chromosome = exon_info[0]
    strand = exon_info[6]

    if "exon_id" not in description:
        raise ValueError('GTF exon entry lacks an exon_id field')
    exon_id = "_".join([chromosome, str(start), str(end), strand])
    exon_name = (description.split("exon_id ")[1]).split('"')[1]
    gene_id = (description.split("gene_id ")[1]).split('"')[1]
    transcript_id = (description.split("transcript_id ")[1]).split('"')[1]
     
    exon = Exon(exon_id, chromosome, start, end, strand, gene_id, \
                transcript_id)
    exon.name = exon_name
    return exon

def get_exon_from_db(exon_row):
    """ Uses information from a database exon entry to create a
    Gene object.

        Args:
            gene_row: Tuple-formatted row from 'exons' table of a
            TALON database
    """
    exon_id = exon_row['identifier']
    exon_name = exon_row['name']
    chromosome = exon_row['chromosome']
    start = exon_row['start']
    end = exon_row['end']
    strand = exon_row['strand']
    gene_id = exon_row['gene_id']

    exon = Exon(exon_id, chromosome, start, end, strand, gene_id, None)
    exon.name = exon_name
    for transcript_id in exon_row['transcript_ids'].split(","):
        exon.transcript_ids.add(transcript_id)
    return exon
