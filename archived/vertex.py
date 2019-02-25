# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

import pdb

class Vertex(object):
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

    def __init__(self, identifier, chromosome, pos, strand, gene_id):
        self.identifier = str(identifier)
        self.chromosome = str(chromosome)
        self.pos = int(pos)
        self.strand = strand
        self.gene_id = str(gene_id)

def fetch_vertex(known_vertices, chromosome, pos, gene_id):
    """ Look for a vertex matching the input criteria. Throw error if not found"""
  
    try:
        vertex_matches = known_vertices[chromosome][pos]
        for v in vertex_matches:
            if v.gene_id == gene_id:
                return v

    except Exception as e:
        raise ValueError('Vertex at ' + chromosome + ':' + str(pos) + 
                         ' does not exist!')
                

def search_for_gene(query_transcript, vertices):
    """ Given a query transcript, the function attempts to find vertices
        that correspond to its exon coordinates. If successful, the ID of the
        gene matching the most vertices is returned. If no vertices are found,
        the function returns None. """

    chromosome = query_transcript.chromosome
    strand = query_transcript.strand
    exon_coords = query_transcript.get_exon_coords()
    genes_seen = []
    for pos in exon_coords:
        try:
            matches = vertices[chromosome][pos]
            for vertex in matches:
                if vertex.strand == strand:
                    genes_seen.append(vertex.gene_id)

        except:
            pass

    if len(genes_seen) == 0:
        return None

    return max(set(genes_seen), key=genes_seen.count)
    


def try_vertex_update(edge, known_vertices, novel_ids, counter):
    """ Given a novel edge, this function determines whether novel vertice(s) are
        needed on one or both ends. If so, new objects are created and added
        to the known_vertices data structure """

    # TODO: Clean up code using fetch_vertex as a helper

    chromosome = edge.chromosome
    gene_id = edge.gene_id
    strand = edge.strand

    if strand == "+":
        v1_pos = edge.start
        v2_pos = edge.end
    if strand == "-":
        v2_pos = edge.start
        v1_pos = edge.end
     
    # First vertex
    found = False

    if chromosome in known_vertices:
        if v1_pos in known_vertices[chromosome]: 
            vertex_matches = known_vertices[chromosome][v1_pos]
            for v in vertex_matches:
                if v.gene_id == gene_id:
                    edge.v1 = v.identifier
                    found = True
    if not found:
        counter["vertices"] += 1
        curr_novel = str(counter["vertices"])
        new_vertex = Vertex(curr_novel, chromosome, v1_pos, strand, gene_id)
        edge.v1 = new_vertex.identifier
        vertex_tuple = (new_vertex.identifier, gene_id, chromosome, v1_pos, 
                        strand)
        novel_ids["vertices"][curr_novel] = vertex_tuple

        if chromosome in known_vertices:
            try:
                known_vertices[chromosome][v1_pos].append(new_vertex)
            except:
                known_vertices[chromosome][v1_pos] = [new_vertex]      
        else:
            known_vertices[chromosome] = {}
            known_vertices[chromosome][v1_pos] = [new_vertex]

    # Second vertex
    found = False 
    if chromosome in known_vertices:
        if v2_pos in known_vertices[chromosome]:
            vertex_matches = known_vertices[chromosome][v2_pos]
            for v in vertex_matches:
                if v.gene_id == gene_id:
                    edge.v2 = v.identifier
                    found = True
    if not found:
        counter["vertices"] += 1
        curr_novel = str(counter["vertices"])
        new_vertex = Vertex(curr_novel, chromosome, v2_pos, strand, gene_id) 
        edge.v2 = new_vertex.identifier
        vertex_tuple = (new_vertex.identifier, gene_id, chromosome, v2_pos,
                        strand)
        novel_ids["vertices"][curr_novel] = vertex_tuple
        if chromosome in known_vertices:
            try:
                known_vertices[chromosome][v2_pos].append(new_vertex)
            except:
                known_vertices[chromosome][v2_pos] = [new_vertex]
        else:
            known_vertices[chromosome] = {}
            known_vertices[chromosome][v2_pos] = [new_vertex]

    return 
