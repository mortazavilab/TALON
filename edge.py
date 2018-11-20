# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

class Edge(object):
    """Stores information about an edge, including its location
       and the gene/transcript(s) it belongs to.
       Attributes:
           identifier: Accession ID of the edge
           gene: Accession ID of the gene that the edge belongs to
           transcript_ids: Set of transcript accession IDs that the edge 
           belongs to
           chromosome: Chromosome that the transcript is located on 
           (format "chr1")
           start: The start position of the edge with respect to the
           forward strand 
           end: The end position of the edge with respect to the
           forward strand
           strand: "+" if the edge is on the forward strand, and "-" if
           it is on the reverse strand
 
           length: The length of the edge
    """

    def __init__(self, identifier, chromosome, start, end, strand, gene_id,
                 transcript_id, annotations):
        self.chromosome = str(chromosome)
        self.gene_id = gene_id
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.length = abs(self.end - self.start + 1)
        self.annotations = annotations

        self.identifier = str(identifier)
        self.transcript_ids = set()
        if transcript_id != None:
            self.transcript_ids.add(transcript_id)
        self.v1 = None
        self.v2 = None

    def print_edge(self):
        """ Prints a string representation of the edge"""
        print self.identifier + ": " + self.chromosome + ":" + \
              str(self.start) + "-" + str(self.end)
        print self.transcript_ids
        return

def create_edge_from_gtf(edge_info):
    """ Creates an edge object using information from a GTF entry
            Args:
               edge_info: A list containing fields from a GTF file edge entry.
               Example:   
               ['chr1', 'HAVANA', 'exon', '11869', '12227', '.', '+', '.', 
                'gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; 
                gene_type "transcribed_unprocessed_pseudogene"; 
                gene_status "KNOWN"; gene_name "DDX11L1"; 
                transcript_type "processed_transcript"; 
                transcript_status "KNOWN"; transcript_name "DDX11L1-002"; 
                edge_number 1; edge_id "ENSE00002234944.1"; level 2; 
                tag "basic"; transcript_support_level "1"; 
                havana_gene "OTTHUMG00000000961.2"; 
                havana_transcript "OTTHUMT00000362751.1";'] 
    """
    description = edge_info[-1]
    start = int(edge_info[3])
    end = int(edge_info[4])
    chromosome = edge_info[0]
    strand = edge_info[6]

    annotations = extract_edge_annotations_from_GTF(edge_info)
    if "exon_id" not in annotations:
        annotations["exon_id"] = "_".join([chromosome, str(start), str(end), strand])
    gene_id = annotations['gene_id']
    transcript_id = annotations['transcript_id']
    edge_id = "_".join([chromosome, str(start), str(end), strand])

    if "gene_id" in description:
        gene_id = (description.split("gene_id ")[1]).split('"')[1]
    if "transcript_id" in description:
        transcript_id = (description.split("transcript_id ")[1]).split('"')[1]
    
    edge = Edge(edge_id, chromosome, start, end, strand, gene_id, transcript_id,
                annotations)
    return edge

def extract_edge_annotations_from_GTF(tab_fields):
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
    if "transcript_id" not in attributes:
        attributes["transcript_id"] = "NULL"

    attributes["source"] = tab_fields[1]

    return attributes

def get_edge_from_db(vertex_info_1, vertex_info_2):
    """ Uses information from a database edge entry to create an edge object.
    """
    if vertex_info_1["edge_id"] != vertex_info_2["edge_id"]:
        raise ValueError('Tried to create edge from endpoints with different IDs')
    edge_id = vertex_info_1["edge_id"]
    chromosome = vertex_info_1['chromosome']
    start = min(vertex_info_1['position'], vertex_info_2['position'])
    end = max(vertex_info_1['position'], vertex_info_2['position']) 
    strand = vertex_info_1['strand']
    gene_id = vertex_info_1['gene_id']

    edge = Edge(edge_id, chromosome, start, end, strand, gene_id, None, None)
    edge.v1 = str(vertex_info_1["vertex_ID"])
    edge.v2 = str(vertex_info_2["vertex_ID"])
    return edge

def create_novel_edge(chromosome, start, end, strand, gene_id, transcript_id, counter):
    """ Creates a novel edge with a unique identifier (obtained using
        counter). Returns the edge object as well as the updated counter.
    """
    counter["edges"] += 1
    curr_novel = counter["edges"]
    edge = Edge(curr_novel, chromosome, start, end, strand, gene_id, transcript_id,
                None)
    return edge
