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

            novel: A counter that keeps track of how many novel transcripts have
            been added to this gene (so that an appropriate ID can be made for
            the next one)

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
        self.novel = 0        

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
        if transcript.gene_id == self.identifier:
            # In order to belong to a gene, the transcript gene_id must
            # match
            transcript_id = transcript.identifier
            self.transcripts[transcript_id] = transcript
        else:
            raise ValueError('Gene ID of transcript must match gene ' + \
                  'in order for assignment to be made.')
        return             

    def lookup_transcript_strict(self, query_transcript):
        """ Checks whether the gene contains a transcript with an exon string
            that matches that of the query_transcript (including exact 3' and 
            5' ends). If yes, the matching transcript is returned. If not, the 
            function returns None.

            Args:
                query_transcript: Transcript object
        """
        query_str = "_".join([str(x) for x in query_transcript.exons])
        for transcript in self.transcripts:
            transcript_str = "_".join([str(x) for x in transcript.exons])
            if query_str == transcript_str:
                return transcript
        return None

    def lookup_transcript_permissive_both(self, query_transcript, verbose):
        """ Checks whether the gene contains a transcript with an exon string
            that matches that of the query_transcript, but allow differences at
            the 5' and 3' end. If yes, the matching transcript is returned. 
            If not, the function returns None.

            Args:
                query_transcript: Transcript object

                verbose: If true, print strings along the way
        """
        query_exons = query_transcript.exons[:]
        query_strand = query_transcript.strand

        query_exons.pop(0)
        query_exons.pop(-1)
        query_str = "_".join([str(x) for x in query_exons])

        if verbose == True:
            print "========================================"
            print query_str
            print "----------------------------------------"

        for transcript_id in self.transcripts:
            transcript = self.transcripts[transcript_id]
            transcript_exons = transcript.exons[:]
            transcript_exons.pop(0)
            transcript_exons.pop(-1)
            transcript_str = "_".join([str(x) for x in transcript_exons])

            if verbose == True:
                print transcript.name
                print transcript_str
                print "................."
            if query_str == transcript_str:
                return transcript
        return None

    def lookup_transcript_permissive5(self, query_transcript):
        """ Checks whether the gene contains a transcript with an exon string
            that matches that of the query_transcript, but allow differences at
            the 5' end. If yes, the matching transcript is returned. If not, the
            function returns None.

            Args:
                query_transcript: Transcript object
        """
        query_exons = query_transcript.exons[:]
        query_strand = query_transcript.strand

        # Find index of 5' end, dependent on strand
        if query_strand == "+":
            index_5prime = 0
        elif query_strand == "-":
            index_5prime = -1
        
        query_exons.pop(index_5prime)
        query_str = "_".join([str(x) for x in query_exons])
        for transcript_id in self.transcripts:
            transcript = self.transcripts[transcript_id]
            transcript_exons = transcript.exons[:]
            transcript_exons.pop(index_5prime)
            transcript_str = "_".join([str(x) for x in transcript_exons])
            if query_str == transcript_str:
                return transcript
        return None

    def lookup_transcript_permissive3(self, query_transcript):
        """ Checks whether the gene contains a transcript with an exon string
            that matches that of the query_transcript, but allow differences at
            the 3' end. If yes, the matching transcript is returned. If not, the
            function returns None.

            Args:
                query_transcript: Transcript object
        """
        query_exons = query_transcript.exons[:]
        query_strand = query_transcript.strand
 
        # Find index of 3'end, dependent on strand
        if query_strand == "+":
            index_3prime = -1
        elif query_strand == "-":
            index_3prime = 0
        query_exons.pop(index_3prime)
        query_str = "_".join([str(x) for x in query_exons])
            
        for transcript_id in self.transcripts:
            transcript = self.transcripts[transcript_id]
            transcript_exons = transcript.exons[:]
            transcript_exons.pop(index_3prime)
            transcript_str = "_".join([str(x) for x in transcript_exons])
            if query_str == transcript_str:
                return transcript
        return None

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
