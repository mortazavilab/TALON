# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

class Transcript(object):
    """Stores information about a gene transcript, including its location
       and constitutive exons.

       Attributes:
           identifier: Accession ID of transcript, i.e. an Ensembl ID. Must
           be unique.

           name: Human-readable name of the transcript. Does not have to be 
           unique

           chromosome: Chromosome that the transcript is located on 
           (format "chr1")

           start: The start position of the transcript with respect to the
           forward strand 

           end: The end position of the transcript with respect to the
           forward strand

           strand: "+" if the transcript is on the forward strand, and "-" if
           it is on the reverse strand

           gene_id: unique ID of the gene that this transcript belongs to

           exons: List of exon objects belonging to this transcript, in sorted
           order.

    """

    def __init__(self, identifier, name, transcript_type, chromosome, start, 
                 end, strand, gene_id):

        self.identifier = identifier
        self.name = name
        self.transcript_type = transcript_type
        self.gene_id = gene_id

        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.n_exons = 0
        self.exons = []


    def get_length(self):
        """ Computes the length of the transcript by summing the lengths of
            its exons """

        if len(self.exons) == 0:
            raise ValueError('Cannot compute length: Transcript does not ' + \
                             'have any exons')
        
        transcript_length = 0
        for exon in self.exons:
            transcript_length += exon.length
        return transcript_length

    def get_exon_coords(self):
        """ Returns a list of the exon coordinates in order """
        exon_coords = []
        for exon in self.exons:
            exon_coords.append(exon.start)
            exon_coords.append(exon.end)
        return exon_coords

    def add_exon(self, exon):
        """Adds an exon object to the transcript."""

        if exon.start > exon.end:
            raise ValueError('Exon start (' + str(exon_start) + ')' + \
                'is supposed to be before the exon end (' + str(exon_end) + ')')

        # Check where in the list the exon should be added
        for i in range(0,len(self.exons)):
            existing_exon = self.exons[i]
            if exon.end < existing_exon.start:
                self.exons = self.exons[0:i] + [exon] + self.exons[i:]
                self.check_exon_validity()
                self.n_exons += 1
                return
        self.exons.append(exon)
        self.check_exon_validity()
        self.n_exons += 1
        return
                    
    def check_exon_validity(self):
        """ The transcript's exons are valid if:
            1) Exons are in sorted order (ascending)
            2) Exon bounds do not exceed transcript start and end
            3) Exons are all on the appropriate chromosome
            If these conditions are violated, this function raises an error.
        """
        prev = 0
        for exon in self.exons:
            if exon.chromosome != self.chromosome:
                raise ValueError('Invalid exon in transcript ' + \
                      self.identifier + ': wrong chromosome')
            if exon.start < self.start or exon.end > self.end:
                raise ValueError('Invalid exon in transcript ' + \
                      self.identifier + ': (' + str(exon.start) + "-" + \
                      str(exon.end) + \
                      ') is located beyond start or end of transcript')
            if exon.start <= prev:
                # This error would indicate a TALON bug rather than user error,
                # so we shouldn't see it. 
                raise ValueError('Exons of transcript ' + \
                      self.identifier + ' are not stored in ascending order.')
            prev = exon.end
        return


    def print_transcript(self):
        """ Print a string representation of the Transcript. Good for debugging
        """
        transcript_id = self.identifier
        if transcript_id == None:
            transcript_id = "Transcript"
        if self.name != None:
            # Include name in output if there is one
            print transcript_id + " (" + self.name + "):"
        else:
            print transcript_id + ":"

        print "\tLocation: " + self.chromosome + ":" + str(self.start) + "-" + \
              str(self.end) + "(" + self.strand + ")"

        # Print exons
        print "\tExons: " + "\n".join([str(x.start) + "-" + str(x.end) for x in self.exons])
        return 

def get_transcript_from_db(transcript_row):
    """ Uses information from a database transcript entry to create a
    Transcript object.

        Args:
            transcript_row: Tuple-formatted row from transcripts table of a 
            TALON database
    """
    transcript_id = transcript_row['identifier']
    name = transcript_row['name']
    chromosome = transcript_row['chromosome']
    start = transcript_row['start']
    end = transcript_row['end']
    strand = transcript_row['strand']
    gene_id = transcript_row['gene_id']

    transcript = Transcript(transcript_id, name, chromosome, start, end, \
                            strand, gene_id)
    return transcript
    

def get_transcript_from_gtf(transcript_info):
    """ Uses information from a GTF-formatted transcript entry to create a
    Transcript object.

        Args:
            transcript_info: A list containing fields from a GTF file gene 
            entry. Example:
          
            chr1	HAVANA	transcript	12010	13670	.	+
            .	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; 
            gene_type "transcribed_unprocessed_pseudogene"; 
            gene_status "KNOWN"; gene_name "DDX11L1"; 
            transcript_type "transcribed_unprocessed_pseudogene"; 
            transcript_status "KNOWN"; transcript_name "DDX11L1-001"; 
            level 2; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; 
            transcript_support_level "NA"; havana_gene "OTTHUMG00000000961.2"; 
            havana_transcript "OTTHUMT00000002844.2";
    """
    chromosome = transcript_info[0]
    description = transcript_info[-1]
    start = int(transcript_info[3])
    end = int(transcript_info[4])
    strand = transcript_info[6]
    transcript_type = None

    name = None
    gene_id = None
    if "transcript_id" not in description:
            raise ValueError('GTF entry lacks a transcript_id field')
    transcript_id = (description.split("transcript_id ")[1]).split('"')[1]

    if "transcript_name" in description:
        name = (description.split("transcript_name ")[1]).split('"')[1]

    if "gene_id" in description:
        gene_id = (description.split("gene_id ")[1]).split('"')[1]

    if "transcript_type" in description:
        transcript_type=(description.split("transcript_type ")[1]).split('"')[1]

    transcript = Transcript(transcript_id, name, transcript_type, chromosome, 
                            start, end, strand, gene_id)
    return transcript

def get_transcript_from_exon(exon, transcript_id):
    """ In rare cases, GTF exons are listed with gene and transcript IDs that
        do not have corresponding entries. In this case, we create a transcript
        for this exon for bookkeeping purposes."""

    gene_id = exon.gene_id
    transcript_id = transcript_id
    name = transcript_id
    chromosome = exon.chromosome
    start = exon.start
    end = exon.end
    strand = exon.strand
    transcript = Transcript(transcript_id, name, None, chromosome, start, end,
                            strand, gene_id)
    return transcript

