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

    def __init__(self, identifier, chromosome, start, end, strand, gene_id, 
                 annotations):

        self.identifier = str(identifier)
        self.gene_id = str(gene_id)

        self.chromosome = str(chromosome)
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.n_exons = 0
        self.exons = []
        self.introns = []
        self.annotations = annotations


    def get_edge_path(self):
        edges = self.get_all_edges()
        if len(edges) == 0:
            return None
        path = [ x.identifier for x in edges]
        return ",".join(path)

    def get_all_edges(self):
        all_edges = []
        for i in range(0,self.n_exons):
            all_edges.append(self.exons[i])
            try:
                all_edges.append(self.introns[i])
            except:
                pass
            
        return all_edges

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
            exon_coords.append(int(exon.start))
            exon_coords.append(int(exon.end))
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

    def add_intron(self, intron):
        """Adds an edge object to the transcript."""

        if intron.start > intron.end:
            raise ValueError('Intron start (' + str(intron.start) + ')' + \
                'is supposed to be before the intron end (' + str(intron.end) + ')')

        # Check where in the list the intron should be added
        for i in range(0,len(self.introns)):
            existing_intron = self.introns[i]
            if intron.end < existing_intron.start:
                self.introns = self.introns[0:i] + [intron] + self.introns[i:]
                return
        self.introns.append(intron)
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

    def get_introns(self):
        """
        Computes introns based on the exon list
        """
        exon_coords = self.get_exon_coords()
        intron_list = []

        i = 1
        while (i < len(exon_coords) - 1):
            j = i + 1

            intron_list.append(exon_coords[i] + 1)
            intron_list.append(exon_coords[j] - 1)
            i += 2

        return intron_list


    def print_transcript(self):
        """ Print a string representation of the Transcript. Good for debugging
        """
        transcript_id = self.identifier
        if transcript_id == None:
            transcript_id = "Transcript"

        print "\tLocation: " + self.chromosome + ":" + str(self.start) + "-" + \
              str(self.end) + "(" + self.strand + ")"

        # Print exons
        print "\tExons: " + "\n".join([str(x.start) + "-" + str(x.end) for x in self.exons])
        return 

def get_transcript_from_db(transcript_row, exon_tree, intron_tree):
    """ Uses information from a database transcript entry to create a
    Transcript object.
        Args:
            transcript_row: Tuple-formatted row from transcripts table of a 
            TALON database
    """
    transcript_id = str(transcript_row['transcript_id'])
    gene_id = str(transcript_row['gene_id'])

    edges = transcript_row['path'].split(",")

    # Check strand
    sample_edge = str(edges[0])
    strand = (exon_tree.edges[sample_edge]).strand

    # Reverse the edge list if the transcript is on the - strand
    if strand == "-":
        edges = edges[::-1]    

    # Get start and end of transcript
    chromosome = (exon_tree.edges[edges[0]]).chromosome
    start = (exon_tree.edges[edges[0]]).start
    end = (exon_tree.edges[edges[-1]]).end

    transcript = Transcript(transcript_id, chromosome, start, end, strand, 
                            gene_id,{})

    # Make sure that all of the exons and introns in this transcript have a 
    # non-zero length. Otherwise, return None
    for i in range(0,len(edges)):
        # Even indices are exons
        if i % 2 == 0:
            curr_exon_id = str(edges[i])
            if curr_exon_id not in exon_tree.edges:
                print "Warning: Ignoring transcript with ID " + transcript_id +\
                " because exon " + curr_exon_id + " not found in exon tree."
                return None
        else:
            curr_intron_id = str(edges[i])
            if curr_intron_id not in intron_tree.edges:
                print "Warning: Ignoring transcript with ID " + transcript_id +\
                " because intron " + curr_intron_id + " not found in intron tree."
                return None

    for i in range(0,len(edges)):
        # Even indices are exons
        if i % 2 == 0:
            curr_exon_id = str(edges[i])
            curr_exon = exon_tree.edges[curr_exon_id]
            transcript.add_exon(curr_exon)
            (curr_exon.transcript_ids).add(transcript_id)
        else:
            curr_intron_id = str(edges[i])
            curr_intron = intron_tree.edges[curr_intron_id]
            transcript.add_intron(curr_intron)
            (curr_intron.transcript_ids).add(transcript_id)

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
    start = int(transcript_info[3])
    end = int(transcript_info[4])
    strand = transcript_info[6]

    if "transcript_id" not in transcript_info[-1]:
            raise ValueError('GTF entry lacks a transcript_id field')
    annotations = extract_transcript_annotations_from_GTF(transcript_info)
    gene_id = annotations['gene_id']
    transcript_id = annotations['transcript_id']

    transcript = Transcript(transcript_id, chromosome, start, end, strand, 
                            gene_id, annotations)

    return transcript

def extract_transcript_annotations_from_GTF(tab_fields):
    """ Extracts key-value annotations from the GTF description field
    """
    attributes = {}

    description = tab_fields[-1].strip()
    # Parse description
    for pair in [x.strip() for x in description.split(";")]:
        if pair == "": continue

        pair = pair.replace('"', '')
        key, val = pair.split()
        attributes[key] = val

    # Put in placeholders for important attributes (such as gene_id) if they
    # are absent
    if "gene_id" not in attributes:
        attributes["gene_id"] = "NULL"

    attributes["source"] = tab_fields[1]

    return attributes    


def get_transcript_from_exon(exon, gene_id, transcript_id):
    """ In rare cases, GTF exons are listed with gene and transcript IDs that
        do not have corresponding entries. In this case, we create a transcript
        for this exon for bookkeeping purposes."""

    name = transcript_id
    chromosome = exon.chromosome
    start = exon.start
    end = exon.end
    strand = exon.strand
    transcript = Transcript(transcript_id, name, None, chromosome, start, end,
                            strand, gene_id)
    return transcript

def create_novel_transcript(chromosome, start, end, strand, gene_id, counter,
                             exons, introns):
    """ Creates a novel transcript with a unique identifier (obtained using
        counter). Returns the transcript object as well as the updated counter.
    """
    counter["transcripts"] += 1
    transcript_id = str(counter["transcripts"])
    
    transcript = Transcript(transcript_id, chromosome, start, end, strand, 
                            gene_id, None)

    for exon in exons:
        transcript.add_exon(exon)
    for intron in introns:
        transcript.add_intron(intron)

    return transcript
