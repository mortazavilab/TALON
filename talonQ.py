# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

import argparse
from functools import reduce
import sqlite3
import dstruct
import operator
import warnings
import transcript_utils as tutils

def get_args():
    """ Fetches the arguments for the program """

    program_desc = """TALON takes transcripts from one or more long read
                      datasets (SAM format) and assigns them transcript and gene 
                      identifiers based on a database-bound annotation. 
                      Novel events are assigned new identifiers."""
    parser = argparse.ArgumentParser(description=program_desc)

    parser.add_argument("--f", dest = "config_file",
        help = "Dataset config file: dataset name, sample description, " + \
               "platform, sam file (comma-delimited)", type = str)
    parser.add_argument('--db', dest = 'database', metavar='FILE,', type = str,
        help='TALON database. Created using build_talon_annotation.py')
    parser.add_argument('--build', dest = 'build', metavar='STRING,', type = str,
        help='Genome build (i.e. hg38) to use. Must be in the database.')
    parser.add_argument("--cov", "-c", dest = "min_coverage",
        help = "Minimum alignment coverage in order to use a SAM entry. Default = 0.9",
        type = str, default = 0.9)
    parser.add_argument("--identity", "-i", dest = "min_identity",
        help = "Minimum alignment identity in order to use a SAM entry. Default = 0",
        type = str, default = 0)
    parser.add_argument("--o", dest = "outprefix", help = "Prefix for output files",
        type = str)

    args = parser.parse_args()
    return args

def str_wrap_double(s):
    """ Adds double quotes around the input string """
    s = str(s)
    return '"' + s + '"'

def make_transcript_dict(cursor):
    """ Format of dict:
            Key: tuple consisting of edges in transcript path
            Value: SQLite3 row from transcript table
         This query returns the novelty designation of the transcript in 
         addition to basic attributes (KNOWN if transcript is known, None
         otherwise).
    """
    transcript_dict = {}
    query = """SELECT * FROM transcripts """
    #query = """ SELECT t.*, 
    #                   ta.value AS annot_status 
    #	            FROM transcripts AS t
    #	            LEFT JOIN transcript_annotations AS ta 
    #                    ON t.transcript_ID = ta.ID
    #	            AND ta.attribute = 'transcript_status' 
    #                AND ta.value = 'KNOWN'
    #        """
    cursor.execute(query)
    for transcript in cursor.fetchall():
        transcript_path = transcript["path"].split(",")
        transcript_path = tuple([ int(x) for x in transcript_path ])
        transcript_dict[transcript_path] = transcript

    return transcript_dict

def make_location_dict(genome_build, cursor):
    """ Format of dict:
            Key: chromosome, pos
            Value: SQLite3 row from location table
    """
    location_dict = {}
    query = """SELECT * FROM location WHERE genome_build = ? """
    cursor.execute(query, [genome_build])
    for location in cursor.fetchall():
        chromosome = location["chromosome"]
        position = location["position"]
        key = (chromosome, position)
        location_dict[key] = location

    return location_dict
    
def make_edge_dict(cursor):
    """ Format of dict:
            Key: vertex1_vertex2_type
            Value: SQLite3 row from edge table
    """
    edge_dict = {}
    query = """SELECT * FROM edge"""
    cursor.execute(query)
    for edge in cursor.fetchall():
        vertex_1 = edge["v1"]
        vertex_2 = edge["v2"]
        edge_type = edge["edge_type"]
        key = (vertex_1, vertex_2, edge_type)
        edge_dict[key] = edge

    return edge_dict

def make_vertex_2_gene_dict(cursor):
    """ Create a dictionary that maps vertices to the genes that they belong to.
    """
    vertex_2_gene = {}
    query = """SELECT * FROM vertex LEFT JOIN genes ON vertex.gene_ID = genes.gene_ID"""
    cursor.execute(query)
    for vertex_line in cursor.fetchall():
        vertex = vertex_line["vertex_ID"]
        gene = vertex_line["gene_ID"]
        strand = vertex_line["strand"]

        if vertex in vertex_2_gene:
            vertex_2_gene[vertex].add((gene, strand))
        else:
            vertex_2_gene[vertex] = set()
            vertex_2_gene[vertex].add((gene, strand))

    return vertex_2_gene

def make_temp_novel_gene_table(cursor, build):
    """ Attaches a temporary database with a table that has the following fields:
            - gene_ID
            - chromosome
            - start
            - end
            - strand
        The purpose is to track novel genes from this run in order to match
        transcripts to them when other forms of gene assignment have failed.
    """
    command = """ CREATE TEMPORARY TABLE IF NOT EXISTS temp_gene AS 
                  SELECT gene_ID,
                    chromosome,
                    start,
                    end,
                    strand
                   FROM (SELECT g.gene_ID,
                             loc.chromosome,
                             MIN(loc.position) as start,
                             MAX(loc.position) as end,
                             g.strand
                       FROM genes as g
                       LEFT JOIN vertex as v ON g.gene_ID = v.gene_ID
                       LEFT JOIN location as loc ON loc.location_ID = v.vertex_ID
                       WHERE loc.genome_build = '%s'
                       GROUP BY g.gene_ID); """

    cursor.execute(command % (build))
    return

def search_for_vertex_at_pos(chromosome, position, location_dict):
    """ Given a chromosome and a position (1-based), this function queries the 
        location dict to determine whether a vertex 
        fitting those criteria exists. Returns the row if yes, and __ if no.
    """
    query_key = (chromosome, position)
    try:
        return location_dict[query_key]
    except:
        return None

def search_for_edge(vertex_1, vertex_2, edge_type, edge_dict):
    """ Search the edge dict for an edge linking vertex_1 and vertex_2"""
    query_key = (vertex_1, vertex_2, edge_type)
    try:
        return edge_dict[query_key]
    except:
        return None

def match_all_transcript_vertices(chromosome, positions, strand, location_dict, 
                                  run_info):
    """ Given a chromosome and a list of positions from the transcript in 5' to
        3' end order, this function looks for a matching vertex for each 
        position. Also returns a list where each index indicates whether that
        vertex is novel to the data structure (0 for known, 1 for novel) """ 
 
    # Returned by function
    vertex_matches = []
    novelty = []
    diff_5p = None
    diff_3p = None

    # helpers
    start = 0
    end = len(positions) - 1

    # Iterate over positions
    for curr_index in range(0,len(positions)):
        position = positions[curr_index]

        # Start and ends require permissive matching approach
        if curr_index == start:
            sj_pos = positions[curr_index + 1]
            pos_type = "start"
            vertex_match, diff_5p = permissive_vertex_search(chromosome, position, 
                                                    strand, sj_pos, pos_type,
                                                    location_dict, run_info)
        elif curr_index == end:
            sj_pos = positions[curr_index - 1]
            pos_type = "end"
            vertex_match, diff_3p = permissive_vertex_search(chromosome, position,
                                                    strand, sj_pos, pos_type,
                                                    location_dict, run_info)

        # Remaining mid-transcript positions go through strict matching process
        else:
            vertex_match = search_for_vertex_at_pos(chromosome, position, 
                                                     location_dict)
        if vertex_match == None:
            # If no vertex matches the position, one is created.
            vertex_match = create_vertex(chromosome, position, run_info, 
                                         location_dict)
            novelty.append(1)
        else:
            novelty.append(0)

        # Add to running list of matches
        vertex_matches.append(vertex_match['location_ID'])

    return tuple(vertex_matches), tuple(novelty), diff_5p, diff_3p

def permissive_vertex_search(chromosome, position, strand, sj_pos, pos_type,
                             locations, run_info):
    """ Given a position, this function tries to find a vertex match within the
        cutoff distance that also comes before the splice junction begins.
        If no vertex is found, the function returns None. """

    # Try a strict match first
    if (chromosome, position) in locations:
        match = locations[(chromosome, position)]
        dist = 0
        return match, dist

    if pos_type != "start" and pos_type != "end":
        raise ValueError("Please set pos_type to either 'start' or 'end'.") 
    if strand != "+" and strand != "-":
        raise ValueError("Invalid strand specified: %s" % strand)

    if pos_type == "start": 
        cutoff = run_info.cutoff_5p
    else:
        cutoff = run_info.cutoff_3p

    # If there is no strict match, look for vertices that are
    #     (1) On the correct chromosome
    #     (2) Are not located at or past the adjoining splice junction
    #     (3) Are within the permitted cutoff distance

    if (strand == "+" and pos_type == "start") or \
       (strand == "-" and pos_type == "end"):
        location_keys = list(filter(lambda t: t[0]==chromosome and
                             abs(t[1] - position) <= cutoff and
                             t[1] < sj_pos,
                             list(locations.keys())))
    else:
        location_keys = list(filter(lambda t: t[0]==chromosome and
                             abs(t[1] - position) <= cutoff and
                             t[1] > sj_pos,
                             list(locations.keys())))

    min_dist = cutoff*2
    closest = None
    for candidate_match in location_keys:
        if strand == "+":
            dist = candidate_match[1] - position
        else:
            dist = position - candidate_match[1]

        if abs(dist) <= abs(min_dist):
            min_dist = dist
            closest = locations[candidate_match]

    if closest == None:
        min_dist = None

    return closest, min_dist

            
def create_vertex(chromosome, position, run_info, location_dict):
    """ Creates a novel vertex and adds it to the location data structure. """
    run_info.vertex += 1
    new_ID = run_info.vertex
    new_vertex = {'location_ID': new_ID,
                  'genome_build': run_info.build,
                  'chromosome': chromosome,
                  'position': position}

    location_dict[(chromosome, position)] = new_vertex

    return new_vertex

def create_edge(vertex_1, vertex_2, edge_type, strand, edge_dict, run_info):
    """ Creates a novel edge and adds it to the edge data structure. """
    run_info.edge += 1
    new_ID = run_info.edge
    new_edge = {'edge_ID': new_ID,
                'v1': vertex_1,
                'v2': vertex_2,
                'edge_type': edge_type,
                'strand': strand }
    
    edge_dict[(vertex_1, vertex_2, edge_type)] = new_edge

    return new_edge

def create_gene(chromosome, start, end, strand, memory_cursor, run_info):
    """ Create a novel gene and add it to the temporary table.
    """
    run_info.genes += 1
    new_ID = run_info.genes

    new_gene = ( new_ID, chromosome, min(start, end), max(start, end), strand )
    cols = " (" + ", ".join([str_wrap_double(x) for x in ["gene_ID", 
           "chromosome", "start", "end", "strand"]]) + ") "
    command = 'INSERT INTO temp_gene ' + cols + ' VALUES ' + '(?,?,?,?,?)'
    memory_cursor.execute(command, new_gene)
    return new_ID

def create_transcript(gene_ID, edge_IDs, vertex_IDs, transcript_dict, run_info):
    """Creates a novel transcript and adds it to the transcript data structure.
    """
    run_info.transcripts += 1
    new_ID = run_info.transcripts
    new_transcript = {'transcript_ID': new_ID,
                      'gene_ID': gene_ID,
                      'path': ",".join(map(str, edge_IDs)),
                      'start_vertex': vertex_IDs[0],
                      'end_vertex': vertex_IDs[-1]}

    transcript_dict[tuple(edge_IDs)] = new_transcript

    return new_transcript

def check_all_exons_known(novelty):
    """ Given a list in which each element represents the novelty (1) or
        known-ness of a transcript edge (0), determine whether all of the
        exons are known or not. Return True if all are known, and False
        otherwise """

    if len(novelty) == 1:
        return novelty[0] == 0

    exons = novelty[::2]

    if sum(exons) > 0:
        return False
    else:
        return True

def check_all_SJs_known(novelty):
    """ Given a list in which each element represents the novelty (1) or 
        known-ness of a transcript edge (0), determine whether all of the
        introns are known or not. Return True if all are known, and False 
        otherwise """

    if len(novelty) == 1:
        return None

    introns = novelty[1::2]
        
    if sum(introns) > 0:
        return False
    else:
        return True

def match_all_transcript_edges(vertices, strand, edge_dict, run_info):
    """ Given a list of vertex IDs from the transcript in 5' to
        3' end order, this function looks for a matching edge ID for each
        position. If none exists, it creates one. """

    edge_matches = []
    novelty = []
    edge_type = "exon"

    for index_1 in range(0, len(vertices) - 1):
        index_2 = index_1 + 1

        if index_1 % 2 != 0:
            edge_type = "intron"
        else:
            edge_type = "exon"
       
        vertex_1 = vertices[index_1]
        vertex_2 = vertices[index_2]
        
        edge_match = search_for_edge(vertex_1, vertex_2, edge_type, edge_dict)
                                                
        if edge_match == None:
            # If no edge matches the position, one is created.
            edge_match = create_edge(vertex_1, vertex_2, edge_type, strand, 
                                     edge_dict, run_info)
            novelty.append(1)
        else:
            novelty.append(0)

        # Add to running list of matches
        edge_matches.append(edge_match['edge_ID'])

    return tuple(edge_matches), tuple(novelty)

def search_for_transcript_suffix(edge_IDs, transcript_dict):
    """ Given a list of edges in a query transcript, determine whether it is
        a suffix for any transcript in the dict. It is OK for the final exon ID
        to be different (and the first one in case there is a 5' end difference), 
        but the splice junctions must match.
    """  

    if len(edge_IDs) > 1:
        suffix_matches = list(filter(lambda t: edge_IDs[1:-1] == t[-len(edge_IDs) + 1:-1],
                                     list(transcript_dict.keys())))
    else:
        suffix_matches = list(filter(lambda t: edge_IDs[0] == t[-1],
                                     list(transcript_dict.keys())))

    if len(suffix_matches) == 0:
        return None, None

    gene_ID = transcript_dict[suffix_matches[0]]["gene_ID"]
    return gene_ID, suffix_matches

def search_for_transcript_prefix(edge_IDs, transcript_dict):
    """ Given a list of edges in a query transcript, determine whether it is
        a prefix for any transcript in the dict. It is OK for the first and 
        last exon IDs to be different, but the splice junctions must match.
    """

    if len(edge_IDs) > 1:
        prefix_matches = list(filter(lambda t: edge_IDs[1:-1] == t[1:len(edge_IDs)-1],
                                     list(transcript_dict.keys())))
    else: 
        prefix_matches = list(filter(lambda t: edge_IDs[0] == t[0],
                                     list(transcript_dict.keys())))

    if len(prefix_matches) == 0:
        return None, None

    gene_ID = transcript_dict[prefix_matches[0]]["gene_ID"]
    return gene_ID, prefix_matches


def search_without_transcript_ends(edge_IDs, transcript_dict):
    """ Search for the body of the query transcript (i.e. leave out the 3' and 
        5' exons). Number of edges in the query and match must be the same. 
    """
  
    try:
        matches = list(filter(lambda t: edge_IDs[1:-1] == t[1:-1],
                                     list(transcript_dict.keys())))
        transcripts = [ transcript_dict[match] for match in matches ]
        gene_ID = transcripts[0]["gene_ID"]
        return gene_ID, transcripts

    except:
        return None,None
    

def search_for_ISM(edge_IDs, transcript_dict):
    """ Given a list of edges in a query transcript, determine whether it is an
        incomplete splice match (ISM) of any transcript in the dict. Similarly
        to the suffix case, we're looking for the gene ID here rather than 
        worrying about exactly which transcript it came from. """               


    if len(edge_IDs) > 1:
        edge_IDs = edge_IDs[1:-1]
        ISM_matches = list(filter(lambda t: set(edge_IDs) <= set(t),
                              list(transcript_dict.keys())))
    else:
        ISM_matches = list(filter(lambda t: set(edge_IDs) <= set(t),
                              list(transcript_dict.keys())))
     

    if len(ISM_matches) == 0:
        return None, None

    else:
        match = ISM_matches[0]
        gene_ID = transcript_dict[match]["gene_ID"]
        return gene_ID, ISM_matches

 
def search_for_overlap_with_gene(chromosome, start, end, strand, 
                                 cursor, run_info):
    """ Given a start and an end value for an interval, query the database to
        determine whether the interval overlaps with any genes. If it there is
        more than one match, prioritize same-strand first and foremost. 
        If there is more than one same-strand option, prioritize amount of
        overlap. Antisense matches may be returned if there is no same strand
        match. """

    min_start = min(start, end)
    max_end = max(start, end)
    query_interval = [min_start, max_end]

    query = """ SELECT gene_ID,
                       chromosome,
                       MIN(start) AS start,
                       MAX(end) AS end,
                       strand
                FROM temp_gene
                WHERE (chromosome = '%s') AND
                      ((start <= %d AND end >= %d) OR
                      (start >= %d AND end <= %d) OR
                      (start >= %d AND start <= %d) OR
                      (end >= %d AND end <= %d))
                 GROUP BY gene_ID;"""    


    cursor.execute(query % (chromosome, min_start, max_end,
                            min_start, max_end, min_start, max_end, min_start,
                            max_end))
    matches = cursor.fetchall()

    if len(matches) == 0:
        return None, None
    
    # Among multiple matches, preferentially return the same-strand gene with 
    # the greatest amount of overlap
    same_strand_matches = len([ x for x in matches if x["strand"] == strand ])

    if strand == "+" and same_strand_matches > 0 or \
        strand == "-" and same_strand_matches == 0:

        matches = [ x for x in matches if x["strand"] == "+" ]
        best_match = get_best_match(matches, query_interval)

    else:   
        matches = [ x for x in matches if x["strand"] == "-" ]
        best_match = get_best_match(matches, query_interval)

    return best_match['gene_ID'], best_match['strand']

def get_best_match(matches, query_interval):
    """ Given a set of gene matches and a query interval, return the match
        that has the greatest amount of overlap with the query."""

    max_overlap = 0
    best_match = None

    for match in matches:
        match_interval = [ match['start'], match['end'] ]
        overlap = get_overlap(query_interval, match_interval)
        if overlap >= max_overlap:
            max_overlap = overlap
            best_match = match

    return best_match            

def get_overlap(a, b):
    """ Computes the amount of overlap between two intervals.
        Returns 0 if there is no overlap. The function treats the start and
        ends of each interval as inclusive, meaning that if a = b = [10, 20],
        the overlap reported would be 11, not 10.
        Args:
            a: First interval, formattted as a list
            b: Second interval, formatted as a list
    """
    overlap = max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
    return overlap


def search_for_transcript(edge_IDs, transcript_dict):
    """ Given the edge IDs that make up a query transcript (5' to 3' order), 
        look for a match in the transcript dict. 
        Return gene ID and transcript ID if found, and None if not. """

    try:
        transcript = transcript_dict[edge_IDs]
        gene_ID = transcript["gene_ID"]
        return gene_ID, transcript

    except:
        return None, None

def process_FSM(edge_IDs, vertex_IDs, transcript_dict, run_info ):
    """ Given a transcript, try to find an FSM gene and transcript match for it. 
        Also delineate its type(s) of 5' and 3' novelty if applicable, which 
        will later be added to the database. In the case of an FSM with end 
        novelty, a novel transcript will be created. Returns None if no FSM
        matches are found."""

    novelty = []

    # Look for exact FSM first
    gene_ID, transcript_match = search_for_transcript(edge_IDs, transcript_dict)
    if transcript_match != None:
        transcript_ID = transcript_match['transcript_ID']
        return gene_ID, transcript_ID, []

    # At this point, return None for a monoexonic transcript, because different
    # ends are not allowed
    if len(edge_IDs) == 1:
        return None, None, novelty
    
    # Next, look for FSM with different ends. This triggers a novel transcript
    # that has novel 5' or 3' ends
    gene_ID, transcript_matches = search_without_transcript_ends(edge_IDs,
                                                               transcript_dict)
    if transcript_matches != None:
        # Create a new transcript
        novel_transcript = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                             transcript_dict, run_info)

        transcript_IDs = ",".join([str(match["transcript_ID"]) for match in transcript_matches])
        novelty.append((novel_transcript['transcript_ID'], run_info.idprefix, 
                        "TALON", "FSM_transcript", "TRUE"))
        novelty.append((novel_transcript['transcript_ID'], run_info.idprefix, 
                        "TALON", "related_transcript_IDs", transcript_IDs))
                
        return gene_ID, novel_transcript["transcript_ID"], novelty

    return None,None,novelty

def process_ISM(edge_IDs, vertex_IDs, transcript_dict, run_info ):
    """ Given a transcript, find the best gene and transcript match for it. 
        Also delineate its type(s) of novelty, which will later be added to 
        the database. If the transcript is identified as an ISM, a novel 
        transcript will be created."""

    suffix_ISM = False
    prefix_ISM = False

    gene_ID, ISM_matches = search_for_ISM(edge_IDs, transcript_dict)
    if ISM_matches == None:
        return None,None,[]

    # Look for suffix ISM novelty
    suffix_gene_ID, suffix_matches = search_for_transcript_suffix(edge_IDs, 
                                                               transcript_dict)
    if suffix_matches != None:
        suffix_ISM = True

    # Look for prefix ISM novelty
    prefix_gene_ID, prefix_matches = search_for_transcript_prefix(edge_IDs,
                                                               transcript_dict)
    if prefix_matches != None:
        prefix_ISM = True

    # Create a new transcript. Must wait until here to do it because otherwise 
    # the transcript can find suffix or prefix matches to itself!
    novel_transcript = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]
    ISM_IDs = [str(match[0]) for match in ISM_matches]
    novelty = [(novel_transcript, run_info.idprefix, "TALON", "ISM_transcript","TRUE")]  
 
    # Add novelty types
    if suffix_ISM == True:
        suffix_IDs = [str(match[0]) for match in suffix_matches]
        ISM_IDs += suffix_IDs
        novelty.append(((novel_transcript, run_info.idprefix, "TALON",
                         "ISM-suffix_transcript","TRUE"))) 
    if prefix_ISM == True:
        prefix_IDs = [str(match[0]) for match in prefix_matches]
        ISM_IDs += prefix_IDs
        novelty.append(((novel_transcript, run_info.idprefix, "TALON",
                         "ISM-prefix_transcript","TRUE")))
    
    ISM_IDs = ",".join(set(ISM_IDs))
    novelty.append(((novel_transcript, run_info.idprefix, "TALON",
                         "related_transcript_IDs", ISM_IDs)))

    return gene_ID, novel_transcript, novelty
    
def process_NIC(edge_IDs, vertex_IDs, strand, transcript_dict, vertex_2_gene, run_info):
    """ For a transcript that has been determined to be novel in catalog, find
        the proper gene match (documenting fusion event if applicable). To do 
        this, look up each vertex in the vertex_2_gene dict, and keep track of all
        same-strand genes. """

    gene_matches = []
    
    for vertex in vertex_IDs:
        if vertex in vertex_2_gene:
            curr_matches = vertex_2_gene[vertex]

            # Make sure the gene is on the correct strand
            gene_matches += [ x[0] for x in list(curr_matches) if x[1] == strand ]

    # Now count up how often we see each gene
    gene_tally = dict((x,gene_matches.count(x)) for x in set(gene_matches))

    # TODO: deal with fusions

    # For the main assignment, pick the gene that is observed the most
    gene_ID = max(gene_tally, key=gene_tally.get)

    # Create a new transcript of that gene
    novel_transcript = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)
    novelty = [(novel_transcript, run_info.idprefix, "TALON",
                         "NIC_transcript","TRUE")]

    return gene_ID, novel_transcript["transcript_ID"], novelty    


def find_gene_match_on_vertex_basis(vertex_IDs, strand, vertex_2_gene):
    """ Use vertices in a transcript to try to pinpoint the gene it belongs to.
    """

    gene_matches = []
    for vertex in vertex_IDs: 
        if vertex in vertex_2_gene:
            curr_matches = vertex_2_gene[vertex]

            # Make sure the gene is on the correct strand
            gene_matches += [ x[0] for x in curr_matches if x[1] == strand ]

    if len(gene_matches) == 0:
        return None

    # Now count up how often we see each gene
    gene_tally = dict((x,gene_matches.count(x)) for x in set(gene_matches))
    
    # TODO: deal with fusions

    # For the main assignment, pick the gene that is observed the most
    gene_ID = max(gene_tally, key=gene_tally.get)

    return gene_ID

def update_vertex_2_gene(gene_ID, vertex_IDs, strand, vertex_2_gene):
    """ Add all vertices with gene pairings to vertex_2_gene dict """

    for vertex in vertex_IDs:
        if vertex in vertex_2_gene:
            vertex_2_gene[vertex].add((gene_ID, strand))
        else:
            vertex_2_gene[vertex] = set()
            vertex_2_gene[vertex].add((gene_ID, strand))

    return
            

def identify_transcript(chrom, positions, strand, cursor, location_dict, edge_dict,
                        transcript_dict, vertex_2_gene, run_info):
    """ Inputs:
        - Information about the query transcript
          - chromosome
          - list of positions
          - strand
        - Data structures
          - location_dict (position --> vertex)
          - edge_dict (v1_v2_edgetype --> edge)
          - transcript_dict
          - vertex_2_gene (maps vertices to the gene(s) they are part of)
          - run_info

       Outputs:
          - Assigned gene ID
          - Assigned transcript ID
          - gene and transcript novelty entries (to be added to database)
          - IDs of start and end vertices
          - 5' and 3' deltas from assigned start/end vertices
    """
    gene_novelty = []
    transcript_novelty = []
    n_exons = len(positions)/2.0
 
    # Get vertex matches for the transcript positions
    vertex_IDs, v_novelty, diff_5p, diff_3p = match_all_transcript_vertices(
                                                                 chrom, 
                                                                 positions,
                                                                 strand,
                                                                 location_dict,
                                                                 run_info)

    # Get edge matches for transcript exons and introns based on the vertices
    edge_IDs, e_novelty = match_all_transcript_edges(vertex_IDs, strand,
                                                     edge_dict, run_info)

    # Check novelty of exons and splice jns. This will help us categorize 
    # what type of novelty the transcript has
    all_SJs_known = check_all_SJs_known(e_novelty)
    all_exons_known = check_all_exons_known(e_novelty)
    splice_vertices_known = (sum(v_novelty[1:-1]) == 0)
    all_exons_novel = (reduce(operator.mul, e_novelty, 1) == 1)

    # Look for FSM or ISM. This includes monoexonic cases where the exon is 
    # known
    if all_SJs_known or (n_exons == 1 and all_exons_known):
        #print("Transcript is either an FSM or an ISM")
        # Look for FSM first
        gene_ID, transcript_ID, transcript_novelty = process_FSM(edge_IDs, vertex_IDs, 
                                                             transcript_dict,
                                                             run_info) 
        if gene_ID == None:
            # Look for ISM
            #print("No FSM found, so looking for ISM...")
            gene_ID, transcript_ID, transcript_novelty = process_ISM(edge_IDs, vertex_IDs,
                                                             transcript_dict,
                                                             run_info)   
        # There are rare NIC cases where all of the edges are known, but the
        # transcript arrangement is still novel (i.e. Map2k4 case)   
        if gene_ID == None:
            #print("No ISM found, so looking for NIC...")
            gene_ID, transcript_ID, transcript_novelty = process_NIC(edge_IDs,
                                                                 vertex_IDs,
                                                                 strand,
                                                                 transcript_dict,
                                                                 vertex_2_gene,
                                                                 run_info)
           
    # Novel in catalog transcripts have known splice donors and acceptors,
    # but new connections between them. 
    elif splice_vertices_known and n_exons > 1 and not(all_exons_novel):
        #print("Transcript is definitely Novel in Catalog (NIC)")
        gene_ID, transcript_ID, transcript_novelty = process_NIC(edge_IDs, 
                                                                 vertex_IDs, 
                                                                 strand, 
                                                                 transcript_dict, 
                                                                 vertex_2_gene, 
                                                                 run_info)
    
    # Antisense transcript with splice junctions matching known gene
    elif splice_vertices_known and n_exons > 1:
        if strand == "+":
            anti_strand = "-"
        else:
            anti_strand = "+"
        anti_gene_ID = find_gene_match_on_vertex_basis(vertex_IDs, anti_strand, 
                                                       vertex_2_gene)
        gene_ID = create_gene(chrom, positions[0], positions[-1], 
                              strand, cursor, run_info) 
        transcript_ID = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]

        gene_novelty = [(gene_ID, run_info.idprefix, "TALON",
                         "antisense_gene","TRUE")]
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "related_gene_IDs", anti_gene_ID))
        transcript_novelty = [(transcript_ID, run_info.idprefix, "TALON", 
                               "antisense_transcript", "TRUE")]

    # Novel not in catalog transcripts contain new splice donors/acceptors
    # and contain at least one splice junction. They may belong to an existing
    # gene, but not necessarily.
    elif not(splice_vertices_known) and n_exons > 1 and not(all_exons_novel): 
        #print("Transcript is definitely Novel Not in Catalog (NNC)")
        gene_ID = find_gene_match_on_vertex_basis(vertex_IDs, strand, vertex_2_gene)
        transcript_ID = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]
        transcript_novelty = [(transcript_ID, run_info.idprefix, "TALON",
                               "NNC_transcript", "TRUE")]

    # Transcripts that don't match the previous categories end up here
    else:
        #print("Transcript is genomic and/or antisense")
        gene_ID, match_strand = search_for_overlap_with_gene(chrom, positions[0],
                                                             positions[1], strand, 
                                                             cursor, run_info)
        if gene_ID == None:
            gene_ID = create_gene(chrom, positions[0], positions[-1],
                              strand, cursor, run_info)
            gene_novelty = [(gene_ID, run_info.idprefix, "TALON",
                         "intergenic_novel","TRUE")]

        elif match_strand != strand:
            anti_gene_ID = gene_ID
            gene_ID = create_gene(chrom, positions[0], positions[-1],
                              strand, cursor, run_info)
            gene_novelty = [(gene_ID, run_info.idprefix, "TALON",
                         "antisense_gene","TRUE")]
            gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "related_gene_IDs",anti_gene_ID))

        transcript_ID = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]
        if match_strand != strand:
            transcript_novelty = [(transcript_ID, run_info.idprefix, "TALON",
                                  "antisense_transcript", "TRUE")]
        else:
            transcript_novelty = [(transcript_ID, run_info.idprefix, "TALON",
                                  "genomic_transcript", "TRUE")]
       
    # Add all novel vertices to vertex_2_gene now that we have the gene ID
    update_vertex_2_gene(gene_ID, vertex_IDs, strand, vertex_2_gene)
 
    # Process 5' and 3' end novelty on relevant transcripts
    if v_novelty[0] == 1:
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                  "5p_novel", "TRUE"))
        transcript_novelty.append((transcript_ID,"5p_novelty", None))
    if v_novelty[-1] == 1:
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                  "3p_novel", "TRUE"))

    # For novel genes and transcripts, add names to novelty entries
    if len(gene_novelty) > 0:
        gene_name = run_info.idprefix + "-gene_%d" % gene_ID
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_name", gene_name, None))
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_id", gene_name, None))
    if len(transcript_novelty) > 0:
        transcript_name = run_info.idprefix + "-transcript_%d" % transcript_ID
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_name", transcript_name, None))
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_id", transcript_name, None))

    # Add annotation entries for any novel exons
    exon_novelty = []
    if not all_exons_known:
        for exon,is_novel in zip(edge_IDs, e_novelty):
            if is_novel:
                exon_novelty.append((exon, run_info.idprefix, "TALON", 
                                     "exon_status", "NOVEL", None))

    # Package up information for output
    annotations = {'gene_ID': gene_ID,
                   'transcript_ID': transcript_ID,
                   'gene_novelty': gene_novelty,
                   'transcript_novelty': transcript_novelty,
                   'exon_novelty': exon_novelty,
                   'start_vertex': vertex_IDs[0],
                   'end_vertex': vertex_IDs[-1],
                   'start_delta': diff_5p,
                   'end_delta': diff_3p}
    
    return annotations

def check_inputs(options):
    """ Checks the input options provided by the user and makes sure that
        they are valid. Throw an error with descriptive help message if not."""
    # TODO: add tests to suite

    # Make sure that the genome build exists in the provided TALON database.
    conn = sqlite3.connect(options.database)
    cursor = conn.cursor()
    cursor.execute(""" SELECT DISTINCT name FROM genome_build """)
    builds = [ str(x[0]) for x in cursor.fetchall() ]
    if options.build not in builds:
        build_names = ", ".join(list(annot_builds))
        raise ValueError("Please specify a genome build that exists in the" +
                          " database. The choices are: " + build_names)

    # Make sure that each input dataset is not already in the database, and
    # also make sure that each dataset name is unique
    sam_files = []
    dataset_metadata = []
    curr_datasets = []

    cursor.execute(""" SELECT dataset_name FROM dataset """)
    existing_datasets = [ str(x[0]) for x in cursor.fetchall() ]

    with open(options.config_file, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            curr_sam = line[3]
            if len(line) != 4:
                raise ValueError('Incorrect number of comma-separated fields'+ \
                                 ' in config file. There should be four: ' + \
                                 '(dataset name, sample description, ' + \
                                 'platform, associated sam file).')

            metadata = (line[0], line[1], line[2])
            dataname = metadata[0]
            if dataname in existing_datasets:
                warnings.warn("Ignoring dataset with name '" + dataname + \
                              "' because it is already in the database.")
            elif dataname in curr_datasets:
                warnings.warn("Skipping duplicated instance of dataset '" + \
                               dataname + "'.")
            elif curr_sam in sam_files:
                warnings.warn("Skipping duplicated instance of sam file '" + \
                               curr_sam  + "'.")
            else:
                dataset_metadata.append(metadata)
                curr_datasets.append(dataname)
                if not curr_sam.endswith(".sam"):
                    raise ValueError('Last field in config file must be a .sam file')
                sam_files.append(curr_sam)      

    conn.close()
    return sam_files, dataset_metadata


def init_run_info(cursor, genome_build, min_coverage = 0.9, min_identity = 0):
    """ Initializes a dictionary that keeps track of important run information
        such as the desired genome build, the prefix for novel identifiers,
        and the novel counters for the run. """

    run_info = dstruct.Struct()
    run_info.build = genome_build
    run_info.min_coverage = min_coverage
    run_info.min_identity = min_identity

    # Fetch information from run_info table
    cursor.execute("""SELECT * FROM run_info""")
    for info in cursor.fetchall():
        info_name = info['item']
        value = info['value']
        if info_name != "idprefix":
            value = int(value)
        run_info[info_name] = value

    # Fetch counters
    query = "SELECT * FROM counters WHERE category != 'genome_build'"
    cursor.execute(query)

    for counter in cursor.fetchall():
        counter_name = counter['category']
        run_info[counter_name] = counter['count']

    return run_info

def prepare_data_structures(cursor, build, min_coverage, min_identity):
    """ Initializes data structures needed for the run and organizes them
        in a dictionary for more ease of use when passing them between functions
    """

    make_temp_novel_gene_table(cursor, build)
    run_info = init_run_info(cursor, build, min_coverage, min_identity)
    location_dict = make_location_dict(build, cursor)
    edge_dict = make_edge_dict(cursor)
    transcript_dict = make_transcript_dict(cursor)
    vertex_2_gene = make_vertex_2_gene_dict(cursor)
   
    struct_collection = dstruct.Struct() 
    struct_collection['run_info'] = run_info
    struct_collection['location_dict'] = location_dict
    struct_collection['edge_dict'] = edge_dict
    struct_collection['transcript_dict'] = transcript_dict
    struct_collection['vertex_2_gene'] = vertex_2_gene 

    return struct_collection

def process_all_sam_files(sam_files, dataset_list, cursor, struct_collection, 

                          outprefix):
    """ Iterates over the provided sam files. """

    novel_datasets = []
    all_observed_transcripts = []
    all_gene_annotations = []
    all_transcript_annotations = []
    all_exon_annotations = []
    all_abundance = []

    # Initialize QC output file
    qc_file = outprefix + "_talon_QC.log"
    o = open(qc_file, 'w')
    o.write("# TALON run filtering settings:\n")
    o.write("# Fraction read aligned: " + \
            str(struct_collection.run_info.min_coverage) + "\n")
    o.write("# Min read identity to reference: " + \
            str(struct_collection.run_info.min_identity) + "\n")
    o.write("# Min transcript length: " + \
            str(struct_collection.run_info.min_length) + "\n")
    o.write("-------------------------------------------\n")
    o.write("\t".join(["dataset", "read_ID", "passed_QC", "primary_mapped", 
                       "read_length", "fraction_aligned", "identity"]) + "\n")

    for sam, d_metadata in zip(sam_files, dataset_list):

        # Create annotation entry for this dataset
        struct_collection.run_info['dataset'] += 1
        d_id = struct_collection.run_info['dataset']     
        novel_datasets += [(d_id, d_metadata[0], d_metadata[1], d_metadata[2])]

        # Now process the current sam file
        observed_transcripts, gene_annotations, transcript_annotations, \
        exon_annotations, abundance = annotate_sam_transcripts(sam, d_id, cursor, struct_collection, o)
 
        # Consolidate the outputs
        all_observed_transcripts += observed_transcripts
        all_gene_annotations += gene_annotations
        all_transcript_annotations += transcript_annotations
        all_exon_annotations += exon_annotations
        all_abundance += abundance

    o.close()
    #print(all_transcript_annotations)
    #print(all_gene_annotations)
    
    return novel_datasets, all_observed_transcripts, all_gene_annotations, \
           all_transcript_annotations, all_abundance

def annotate_sam_transcripts(sam_file, dataset, cursor, struct_collection, QC_file):
    """ Process SAM transcripts and annotate the ones that pass QC """

    observed_transcripts = []
    gene_annotations = []
    transcript_annotations = []
    exon_annotations = []
    abundance = {}

    with open(sam_file) as sam:
        for line in sam:
            line = line.strip()

            # Ignore header
            if line.startswith("@"):
                continue

            # Check whether we should try annotating this read or not
            qc_metrics = check_read_quality(line, struct_collection)
            passed_qc = qc_metrics[1]
            QC_file.write("\t".join([str(x) for x in [dataset] + qc_metrics]) \
                          + "\n")

            if not passed_qc:
                continue

            # For transcripts that pass QC, parse the attributes to 
            # determine the chromosome, positions, and strand of the transcript
            try:
                read_ID, chrom, positions, strand, read_length = parse_transcript(line) 
            except:
                message = "Problem parsing transcript '%s'. Skipping.." \
                           % line.split("\t")[0]
                warnings.warn(message)
                continue

            # Now identify the transcript
            location_dict = struct_collection.location_dict
            edge_dict = struct_collection.edge_dict
            transcript_dict = struct_collection.transcript_dict
            vertex_2_gene = struct_collection.vertex_2_gene
            run_info = struct_collection.run_info
            try:
                annotation_info = identify_transcript(chrom, positions, strand, 
                                                      cursor, location_dict, 
                                                      edge_dict, transcript_dict, 
                                                      vertex_2_gene, run_info)
            except:
                warnings.warn("Problem identifying transcript '%s'. Skipping.."\
                               % read_ID)
                exit()
                continue
                            
            # Now that transcript has been annotated, unpack values and 
            # create an observed entry and abundance record
            gene_ID = annotation_info['gene_ID']
            transcript_ID = annotation_info['transcript_ID']
            gene_novelty = annotation_info['gene_novelty']
            transcript_novelty = annotation_info['transcript_novelty']
            exon_novelty = annotation_info['exon_novelty']
            start_vertex = annotation_info['start_vertex']
            end_vertex = annotation_info['end_vertex']
            start_delta = annotation_info['start_delta']
            end_delta = annotation_info['end_delta']

            struct_collection.run_info['observed'] += 1
            obs_ID = struct_collection.run_info['observed'] 
            observed = (obs_ID, gene_ID, transcript_ID, read_ID, dataset, 
                        start_vertex, end_vertex, start_delta, 
                        end_delta, read_length)
            observed_transcripts.append(observed)

            # Also add transcript to abundance dict
            try:
                abundance[transcript_ID] += 1
            except:
                abundance[transcript_ID] = 1

            # Update annotation records
            gene_annotations += gene_novelty
            transcript_annotations += transcript_novelty
            exon_annotations += exon_novelty

    # Before returning abundance, reformat it as database rows
    abundance_rows = []
    for transcript, count in abundance.items():
        curr_row = (transcript, dataset, count)
        abundance_rows.append(curr_row)
    
    return observed_transcripts, gene_annotations, transcript_annotations, \
           exon_annotations, abundance_rows
            

def parse_transcript(sam_read):
    """ Given a SAM transcript, parse the entry to obtain:
           - read ID
           - chromosome
           - positions (start, splice sites, end) ordered from 5' to 3'
           - strand
           - read length
    """
    sam = sam_read.split("\t")
    read_ID = sam[0]
    flag = int(sam[1])
    chromosome = sam[2]
    sam_start = int(sam[3]) # Start is earliest alignment position
    cigar = sam[5]
    seq = sam[9]
    read_length = len(seq) 
    other_fields = sam[11:]

    # Compute attributes
    sam_end = tutils.compute_transcript_end(sam_start, cigar) 
    intron_list = tutils.get_introns(other_fields, sam_start, cigar)
    # Adjust intron positions by 1 to get splice sites
    splice_sites = [x+1 if i%2==1 else x-1 for i,x in enumerate(intron_list)]
    positions = [sam_start] + splice_sites + [sam_end]

    # Flip the positions order if the read is on the minus strand
    if flag in [16, 272]:
        strand = "-"
        positions = positions[::-1]
    else:
        strand = "+"

    return read_ID, chromosome, positions, strand, read_length
    

def check_read_quality(sam_read, struct_collection):
    """ Process an individual sam read and return quality attributes. """

    sam = sam_read.split("\t")
    read_ID = sam[0]
    flag = sam[1]
    cigar = sam[5]
    seq = sam[9]
    read_length = len(seq)

    # Only use uniquely mapped transcripts
    if flag not in ["0", "16"]:
        return [read_ID, 0, 0, read_length, "NA", "NA"]

    # Only use reads that are greater than or equal to length threshold
    if read_length < struct_collection.run_info.min_length:
        return [read_ID, 0, 1, read_length, "NA", "NA"]

    # Locate the MD field of the sam transcript
    try:
        md_index = [i for i, s in enumerate(sam) if s.startswith('MD:Z:')][0]
    except:
        raise ValueError("SAM transcript %s lacks an MD tag" % read_ID)

    # Only use reads where alignment coverage and identity exceed
    # cutoffs
    coverage = tutils.compute_alignment_coverage(cigar)
    identity = tutils.compute_alignment_identity(sam[md_index], seq)   
 
    if coverage < struct_collection.run_info.min_coverage or \
       identity < struct_collection.run_info.min_identity:
        return [read_ID, 0, 1, read_length, coverage, identity]

    # At this point, the read has passed the quality control
    return [read_ID, 1, 1, read_length, coverage, identity]
    
def update_database(cursor, batch_size, novel_datasets, observed_transcripts, 
                    gene_annotations, transcript_annotations, abundance, 
                    struct_collection):
    """ Adds new entries to the database. """

   
    print("Updating abundance table....")
    batch_add_abundance(cursor, abundance, batch_size) 

    print("Updating observed transcript table...")
    batch_add_observed(cursor, observed_transcripts, batch_size)

    print("Updating gene, transcript, and exon annotations...")
    batch_add_annotations(cursor, gene_annotations, "gene", batch_size)
    batch_add_annotations(cursor, transcript_annotations, "transcript", batch_size)
    batch_add_annotations(cursor, exon_annotations, "exon", batch_size)

def batch_add_annotations(cursor, annotations, annot_type, batch_size):
    """ Add gene/transcript/exon annotations to the appropriate annotation table
    """

    if annot_type not in ["gene", "transcript", "exon"]:
        raise ValueError("When running batch annot update, must specify " + \
                         "annot_type as 'gene', 'exon', or 'transcript'.")

    index = 0
    while index < len(annotations):
        try:
            batch = annotations[index:index + batch_size]
        except:
            batch = annotations[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["ID", "annot_name", "source", "attribute", "value"]]) + ") "
            command = 'INSERT INTO "' + annot_type + '_annotations" ' + cols + \
                      "VALUES " + '(?,?,?,?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)
    return

def batch_add_observed(cursor, observed, batch_size):
    """ Adds observed tuples (obs_ID, gene_ID, transcript_ID, read_name,
        dataset, start_vertex_ID, end_vertex_ID, start_delta, end_delta, 
        read_length) to observed table of database. """

    index = 0
    while index < len(observed):
        try:
            batch = observed[index:index + batch_size]
        except:
            batch = observed[index:]
        index += batch_size

        # Add to database
        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["obs_ID", "gene_ID", "transcript_ID", "read_name",
                    "dataset", "start_vertex_ID", "end_vertex_ID",
                    "start_delta", "end_delta", "read_length"]]) + ") "
            command = 'INSERT INTO "observed"' + cols + \
                      "VALUES " + '(?,?,?,?,?,?,?,?,?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)
    return

def batch_add_abundance(cursor, abundances, batch_size):
    """ Adds abundance tuples (transcript_ID, dataset, count) to the abundance
        table of the database """

    index = 0
    while index < len(abundances):
        try:
            batch = abundances[index:index + batch_size]
        except:
            batch = abundances[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["transcript_id", "dataset", "count"]]) + ") "
            command = 'INSERT INTO "abundance"' + cols + "VALUES " + '(?,?,?)'
            cursor.executemany(command, batch)
        except Exception as e:
            print(e)
    return 


def main():
    """ Runs program """

    options = get_args()
    sam_files, dataset_list = check_inputs(options)

    # Fire up the database connection
    database = options.database
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Input parameters
    build = options.build
    min_coverage = float(options.min_coverage)
    min_identity = float(options.min_identity)
    outprefix = options.outprefix

    # Prepare data structures from the database content
    print("Processing annotation...")
    struct_collection = prepare_data_structures(cursor, build, min_coverage, 
                                                min_identity)

    # TODO: Read and annotate input sam files. Also, write output files.
    print("Processing SAM files...")
    novel_datasets, observed_transcripts, gene_annotations, \
    transcript_annotations, abundance = process_all_sam_files(sam_files, 
                                                                   dataset_list, 
                                                                         cursor, 
                                                              struct_collection, 
                                                                      outprefix)
    # TODO: Update database
    batch_size = 10000
    update_database(cursor, batch_size, novel_datasets, observed_transcripts, 
                    gene_annotations, transcript_annotations, abundance, 
                    struct_collection)

    # Validate database

    # conn.commit()

    # Generate output files if desired

    conn.close()

if __name__ == '__main__':
    main()
