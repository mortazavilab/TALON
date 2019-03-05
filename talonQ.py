# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

import cProfile
import argparse
from functools import reduce
import sqlite3
import sys
import dstruct
import operator
#import os
from pathlib import Path
import warnings
import transcript_utils as tutils
import query_utils as qutils
#talon_path = (os.path.abspath(__file__))
#print(sys.path)
#main_path = "/".join(talon_path.split("/")[0:-2])
#sys.path.append(main_path + "/post-TALON_tools")
#print(sys.path)
#import summarize_datasets

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

def make_gene_start_and_end_dict(cursor):
    """ Format of dicts:
            Key: gene ID from database
            Value: set object containing start vertices (or end vertices) of 
                   KNOWN transcripts from that gene
    """
    gene_starts = {}
    gene_ends = {}
    query = """SELECT DISTINCT gene_ID, 
                               start_vertex, 
                               end_vertex 
               FROM transcripts
               LEFT JOIN transcript_annotations as ta 
                   ON ta.ID = transcripts.transcript_ID
	       WHERE ta.attribute = 'transcript_status' 
                     AND ta.value = 'KNOWN'"""

    cursor.execute(query)
    for entry in cursor.fetchall():
        gene_ID = entry['gene_ID']
        start = entry['start_vertex']
        end = entry['end_vertex']
        try:
            gene_starts[gene_ID].add(start)
        except:
            gene_starts[gene_ID] = set()
            gene_starts[gene_ID].add(start)

        try:
            gene_ends[gene_ID].add(end)
        except:
            gene_ends[gene_ID] = set()
            gene_ends[gene_ID].add(end)

    return gene_starts, gene_ends
           

def make_transcript_dict(cursor, build):
    """ Format of dict:
            Key: tuple consisting of edges in transcript path
            Value: SQLite3 row from transcript table
    """
    transcript_dict = {}
    #query = """SELECT * FROM transcripts """
    query = """SELECT t.*,
	         loc1.chromosome as chromosome,
                 loc1.position as start_pos,
                 loc2.position as end_pos	
	       FROM transcripts AS t
	       LEFT JOIN location as loc1 ON t.start_vertex = loc1.location_ID 
	       LEFT JOIN location as loc2 ON t.end_vertex = loc2.location_ID
	       WHERE loc1.genome_build = ? AND loc2.genome_build = ?;
           """

    cursor.execute(query, [build, build])
    for transcript in cursor.fetchall():
        transcript_path = transcript["path"].split(",")
        transcript_path = frozenset([ int(x) for x in transcript_path ])
        transcript_dict[transcript_path] = transcript

    return transcript_dict

def make_location_dict(genome_build, cursor):
    """ Format of dict:
        chromosome -> dict(position -> SQLite3 row from location table)

        old:
            Key: chromosome, pos
            Value: SQLite3 row from location table
    """
    location_dict = {}
    query = """SELECT * FROM location WHERE genome_build = ? """
    cursor.execute(query, [genome_build])
    for location in cursor.fetchall():
        chromosome = location["chromosome"]
        position = location["position"]
        try:
            location_dict[chromosome][position] = location
        except:
            location_dict[chromosome] = {position: location} 

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
    try:
        return location_dict[chromosome][position]
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
    if chromosome in locations and position in locations[chromosome]:
        match = locations[chromosome][position]
        dist = 0
        return match, dist

    if pos_type != "start" and pos_type != "end":
        raise ValueError("Please set pos_type to either 'start' or 'end'.") 
    if strand != "+" and strand != "-":
        raise ValueError("Invalid strand specified: %s" % strand)

    # If there is no strict match, look for vertices that are
    #     (1) On the correct chromosome
    #     (2) Are not located at or past the adjoining splice junction
    #     (3) Are within the permitted cutoff distance

    if pos_type == "start": 
        max_dist = run_info.cutoff_5p
    else:
        max_dist = run_info.cutoff_3p

    # For + strand start and - strand end, look for a vertex with a smaller
    # position first (since degradtion is more biologically likely).
    # For the + strand, this would be a negative delta, and for the - strand,
    # it would be a positive delta
    if (strand == "+" and pos_type == "start") or \
       (strand == "-" and pos_type == "end"): 
        direction_priority = -1
        search_window_start = position - max_dist
        search_window_end = sj_pos
    else:
        direction_priority = 1
        search_window_start = sj_pos
        search_window_end = position + max_dist

    for dist in range(1,max_dist):
        curr_pos = position + dist*direction_priority
        if curr_pos > search_window_start and curr_pos < search_window_end: 
            match = search_for_vertex_at_pos(chromosome, curr_pos, locations)
            if match != None:
                if strand == "+":
                    if curr_pos < position:
                        return match, dist*(-1)  
                    else:
                        return match, dist  
                else:
                    if curr_pos < position:
                        return match, dist
                    else:
                        return match, dist*(-1)

        curr_pos = position - dist*direction_priority
        if curr_pos > search_window_start and curr_pos < search_window_end:
            match = search_for_vertex_at_pos(chromosome, curr_pos, locations)
            if match != None:
                if strand == "+":
                    if curr_pos < position:
                        return match, dist*(-1)
                    else:
                        return match, dist
                else:
                    if curr_pos < position:
                        return match, dist
                    else:
                        return match, dist*(-1)
    return None, None       
            
def create_vertex(chromosome, position, run_info, location_dict):
    """ Creates a novel vertex and adds it to the location data structure. """
    run_info.vertex += 1
    new_ID = run_info.vertex
    new_vertex = {'location_ID': new_ID,
                  'genome_build': run_info.build,
                  'chromosome': chromosome,
                  'position': position}

    try:
        location_dict[chromosome][position] = new_vertex
    except:
        location_dict[chromosome] = { position: new_vertex }

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
    cols = '("gene_ID", "chromosome", "start", "end", "strand")' 
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
                      'end_vertex': vertex_IDs[-1],
                      'n_exons': (len(edge_IDs) + 1)/2 }

    transcript_dict[frozenset(edge_IDs)] = new_transcript

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

#def search_for_transcript_suffix(edge_IDs, transcript_dict):
#    """ Given a list of edges in a query transcript, determine whether it is
#        a suffix for any transcript in the dict. It is OK for the final exon ID
#        to be different (and the first one in case there is a 5' end difference), 
#        but the splice junctions must match.
#    """  

#    if len(edge_IDs) > 1:
#        #suffix_matches = 
#        suffix_matches = list(filter(lambda t: edge_IDs[1:-1] == t[-len(edge_IDs) + 1:-1],
#                                     list(transcript_dict.keys())))
#    else:
#        suffix_matches = list(filter(lambda t: edge_IDs[0] == t[-1],
#                                     list(transcript_dict.keys())))

#    if len(suffix_matches) == 0:
#        return None, None

#    gene_ID = transcript_dict[suffix_matches[0]]["gene_ID"]
#    return gene_ID, suffix_matches

#def search_for_transcript_prefix(edge_IDs, transcript_dict):
#    """ Given a list of edges in a query transcript, determine whether it is
#        a prefix for any transcript in the dict. It is OK for the first and 
#        last exon IDs to be different, but the splice junctions must match.
#    """

#    if len(edge_IDs) > 1:
#        prefix_matches = list(filter(lambda t: edge_IDs[1:-1] == t[1:len(edge_IDs)-1],
#                                     list(transcript_dict.keys())))
#    else: 
#        prefix_matches = list(filter(lambda t: edge_IDs[0] == t[0],
#                                     list(transcript_dict.keys())))

#    if len(prefix_matches) == 0:
#        return None, None

#    gene_ID = transcript_dict[prefix_matches[0]]["gene_ID"]
#    return gene_ID, prefix_matches


def search_without_transcript_ends(edge_IDs, transcript_dict):
    """ Search for the body of the query transcript (i.e. leave out the 3' and 
        5' exons). Number of edges in the query and match must be the same. 
    """
    edge_IDs = frozenset(edge_IDs[1:-1])  

    try:
        matches = [ transcript_dict[x] for x in transcript_dict if edge_IDs.issubset(x)]

        gene_ID = transcripts[0]["gene_ID"]
        return gene_ID, transcripts

    except:
        return None,None
    

def search_for_ISM(edge_IDs, transcript_dict):
    """ Given a list of edges in a query transcript, determine whether it is an
        incomplete splice match (ISM) of any transcript in the dict."""               


    if len(edge_IDs) > 1:
        edges = frozenset(edge_IDs[1:-1])
    else:
        edges = frozenset(edge_IDs)

    ISM_matches = [ transcript_dict[x] for x in transcript_dict if edges.issubset(x)]

    if len(ISM_matches) > 0:
        return ISM_matches
    else:
        return None

 
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
    """ Given the edge IDs (in set data structure) that make up a query 
        transcript, look for a match in the transcript dict. 
        Return gene ID and transcript ID if found, and None if not. """

    try:
        transcript = transcript_dict[edge_IDs]
        gene_ID = transcript["gene_ID"]
        return gene_ID, transcript

    except:
        return None, None

def process_FSM_or_ISM(edge_IDs, vertex_IDs, transcript_dict, gene_starts,
                       gene_ends, run_info):
    """ Given a transcript, try to find an FSM gene and transcript match for it.
        In the case of an FSM with end novelty, a novel transcript will be
        created. Returns None if no FSM matches are found."""

    gene_ID = None
    transcript_ID = None
    novelty = []
    n_exons = (len(edge_IDs) + 1)/2.0

    # Look for exact FSM first because it's a fast search
    full_edge_set = frozenset(edge_IDs)
    gene_ID, transcript_match = search_for_transcript(full_edge_set, transcript_dict)
    if transcript_match != None:
        transcript_ID = transcript_match['transcript_ID']
        return gene_ID, transcript_ID, []

    # Now extract ISM matches of all kinds and then characterize them
    all_matches = search_for_ISM(edge_IDs, transcript_dict)     

    if all_matches == None:
        return None, None, []


    #FSM = []
    ISM = []
    suffix = []
    prefix = []
    gene_ID = all_matches[0]['gene_ID']

    # Check if the start and end vertices match known vertices for this gene.
    # If so, the transcript is eligible to be an NIC instead of ISM
    start_known = vertex_IDs[0] in gene_starts[gene_ID]
    end_known = vertex_IDs[-1] in gene_ends[gene_ID]

    for match in all_matches:
        transcript_ID = run_info.transcripts + 1
        
        # Transcript is FSM
        if match['n_exons'] == n_exons:
            gene_ID = match['gene_ID']
            #FSM.append(str(match['transcript_ID']))

            transcript_ID = match['transcript_ID']
            novelty = []
            return gene_ID, transcript_ID, novelty

        # Transcript is NIC
        elif start_known and end_known:
            novel_transcript = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)
            transcript_ID = novel_transcript['transcript_ID']
            novelty = [(transcript_ID, run_info.idprefix, "TALON",
                        "NIC_transcript", "TRUE")]  
            return gene_ID, transcript_ID, novelty   

        else:
            # Add ISM
            ISM.append(str(match['transcript_ID']))

            # Single-exon case
            if n_exons == 1:
                match_path = match['path']
                exon = str(edge_IDs[0])
                # Look for prefix
                if match_path.startswith(exon):
                    prefix.append(str(match['transcript_ID']))
                # Look for suffix
                if match_path.endswith(exon):
                    suffix.append(str(match['transcript_ID']))
                    #if len(FSM) == 0:
                    gene_ID = match['gene_ID']   
                continue     
 
            # Multi-exon case
            edge_str = ",".join([str(x) for x in edge_IDs[1:-1]])
            match_without_start = ",".join((match['path']).split(",")[1:])
            match_without_end = ",".join((match['path']).split(",")[:-1])

            # Look for prefix
            if match_without_start.startswith(edge_str):
                prefix.append(str(match['transcript_ID']))

            # Look for suffix    
            if match_without_end.endswith(edge_str):
                #if len(FSM) == 0:
                gene_ID = match['gene_ID']
                suffix.append(str(match['transcript_ID']))
           

    novel_transcript = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)
    transcript_ID = novel_transcript['transcript_ID']

    # 
    #if FSM != []:
    #    FSM_str = ",".join(FSM)
    #    novelty.append((transcript_ID, run_info.idprefix, "TALON", 
    #                    "FSM_transcript", "TRUE"))
    #    novelty.append((transcript_ID, run_info.idprefix, "TALON",
    #                    "FSM_to_IDs", FSM_str))
    #    return gene_ID, transcript_ID, novelty
    #else:
    ISM_str = ",".join(ISM)
    novelty.append((transcript_ID, run_info.idprefix, "TALON", 
                    "ISM_transcript", "TRUE"))
    novelty.append((transcript_ID, run_info.idprefix, "TALON",
                    "ISM_to_IDs", ISM_str))
    if prefix != []:
        prefix_str = ",".join(prefix)
        novelty.append((transcript_ID, run_info.idprefix, "TALON", 
                        "ISM-prefix_transcript", "TRUE"))
        novelty.append((transcript_ID, run_info.idprefix, "TALON",
                        "ISM-prefix_to_IDs", prefix_str))
    if suffix != []:
        suffix_str = ",".join(suffix)
        novelty.append((transcript_ID, run_info.idprefix, "TALON",
                        "ISM-suffix_transcript", "TRUE"))
        novelty.append((transcript_ID, run_info.idprefix, "TALON",
                        "ISM-suffix_to_IDs", suffix_str))   
 
    return gene_ID, transcript_ID, novelty
    
    
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
    transcript_ID = novel_transcript["transcript_ID"]
    novelty = [(transcript_ID, run_info.idprefix, "TALON",
                         "NIC_transcript","TRUE")]

    return gene_ID, transcript_ID, novelty    


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
                        transcript_dict, vertex_2_gene, gene_starts, gene_ends,
                        run_info):
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
          - gene_starts (maps gene IDs to known start vertices)
          - gene_ends (maps gene IDs to known end vertices) 
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

        # Look for FSM first
        gene_ID, transcript_ID, transcript_novelty = process_FSM_or_ISM(edge_IDs, 
                                                                vertex_IDs, 
                                                                transcript_dict,
                                                                gene_starts,
                                                                gene_ends, 
                                                                run_info)

        if gene_ID == None:
            gene_ID, transcript_ID, transcript_novelty = process_NIC(edge_IDs,
                                                                 vertex_IDs,
                                                                 strand,
                                                                 transcript_dict,
                                                                 vertex_2_gene,
                                                                 run_info)
           
    # Novel in catalog transcripts have known splice donors and acceptors,
    # but new connections between them. 
    elif splice_vertices_known and n_exons > 1 and not(all_exons_novel):
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

        # Handle gene annotations
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                             "antisense_gene","TRUE"))
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_antisense_to_IDs",anti_gene_ID))

        # Handle transcript annotations
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON", 
                                   "antisense_transcript", "TRUE"))

    # Novel not in catalog transcripts contain new splice donors/acceptors
    # and contain at least one splice junction.
    elif not(splice_vertices_known) and n_exons > 1 and not(all_exons_novel): 
        gene_ID = find_gene_match_on_vertex_basis(vertex_IDs, strand, vertex_2_gene)
        transcript_ID = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]

        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                               "NNC_transcript", "TRUE"))

    # Transcripts that don't match the previous categories end up here
    else:
        gene_ID, match_strand = search_for_overlap_with_gene(chrom, positions[0],
                                                             positions[1], strand, 
                                                             cursor, run_info)
        if gene_ID == None:
            gene_ID = create_gene(chrom, positions[0], positions[-1],
                              strand, cursor, run_info)

            gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "intergenic_novel","TRUE"))
            transcript_ID = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]
            transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                  "intergenic_transcript", "TRUE"))

        elif match_strand != strand:
            anti_gene_ID = gene_ID
            gene_ID = create_gene(chrom, positions[0], positions[-1],
                              strand, cursor, run_info)
            transcript_ID = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]

            gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "antisense_gene","TRUE"))
            gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_antisense_to_IDs",anti_gene_ID))
            transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                  "antisense_transcript", "TRUE"))
        else:
            transcript_ID = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                         transcript_dict, run_info)["transcript_ID"]
            transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                  "genomic_transcript", "TRUE"))

    # Add all novel vertices to vertex_2_gene now that we have the gene ID
    # TODO: we might be able to run this operation fewer times by screening novelty
    update_vertex_2_gene(gene_ID, vertex_IDs, strand, vertex_2_gene)
 
    # Process 5' and 3' end novelty on relevant transcripts
    #if v_novelty[0] == 1:
    #    transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
    #                              "5p_novel", "TRUE"))
    #if v_novelty[-1] == 1:
    #    transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
    #                              "3p_novel", "TRUE"))

    # For novel genes and transcripts, add names to novelty entries
    if len(gene_novelty) > 0:
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                                   "gene_status", "NOVEL"))
        gene_name = run_info.idprefix + "-gene_%d" % gene_ID
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_name", gene_name))
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_id", gene_name))
    if len(transcript_novelty) > 0:
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                   "transcript_status", "NOVEL"))
        transcript_name = run_info.idprefix + "-transcript_%d" % transcript_ID
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_name", transcript_name))
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_id", transcript_name))

    # Add annotation entries for any novel exons
    exon_novelty = []
    if not all_exons_known:
        for exon,is_novel in zip(edge_IDs, e_novelty):
            if is_novel:
                exon_novelty.append((exon, run_info.idprefix, "TALON", 
                                     "exon_status", "NOVEL"))

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

    # Make sure that the input database exists!
    database = options.database
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)

    # Make sure that the genome build exists in the provided TALON database.
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute(""" SELECT DISTINCT name FROM genome_build """)
    builds = [ str(x[0]) for x in cursor.fetchall() ]
    if options.build not in builds:
        build_names = ", ".join(list(builds))
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

            # Make sure that the sam file exists
            if not Path(curr_sam).exists():
                raise ValueError("SAM file '%s' does not exist!" % curr_sam)

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
    transcript_dict = make_transcript_dict(cursor, build)
    vertex_2_gene = make_vertex_2_gene_dict(cursor)
    gene_starts, gene_ends = make_gene_start_and_end_dict(cursor)   

    struct_collection = dstruct.Struct() 
    struct_collection['run_info'] = run_info
    struct_collection['location_dict'] = location_dict
    struct_collection['edge_dict'] = edge_dict
    struct_collection['transcript_dict'] = transcript_dict
    struct_collection['vertex_2_gene'] = vertex_2_gene 
    struct_collection['gene_starts'] = gene_starts
    struct_collection['gene_ends'] = gene_ends

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
        print("Working with dataset %s..." % d_metadata[0])

        # Create annotation entry for this dataset
        struct_collection.run_info['dataset'] += 1
        d_id = struct_collection.run_info['dataset']     
        d_name = d_metadata[0]
        novel_datasets += [(d_id, d_name, d_metadata[1], d_metadata[2])]

        # Now process the current sam file
        observed_transcripts, gene_annotations, transcript_annotations, \
        exon_annotations, abundance = annotate_sam_transcripts(sam, d_name, cursor, struct_collection, o)
 
        # Consolidate the outputs
        all_observed_transcripts.extend(observed_transcripts)
        all_gene_annotations.extend(gene_annotations)
        all_transcript_annotations.extend(transcript_annotations)
        all_exon_annotations.extend(exon_annotations)
        all_abundance.extend(abundance)

    o.close()
    
    return novel_datasets, all_observed_transcripts, all_gene_annotations, \
           all_transcript_annotations, all_exon_annotations, all_abundance

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
            gene_starts = struct_collection.gene_starts
            gene_ends = struct_collection.gene_ends
            run_info = struct_collection.run_info
            #try:
            annotation_info = identify_transcript(chrom, positions, strand, 
                                                      cursor, location_dict, 
                                                      edge_dict, transcript_dict, 
                                                      vertex_2_gene, 
                                                      gene_starts, gene_ends,
                                                      run_info)
            #except:
            #    warnings.warn("Problem identifying transcript '%s'. Skipping.."\
            #                   % read_ID)    
            #    continue
                            
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
    
def update_database(cursor, batch_size, datasets, observed_transcripts, 
                    gene_annotations, transcript_annotations, exon_annotations,
                    abundance, struct_collection):
    """ Adds new entries to the database. """

    print("Adding novel genes to database...")
    add_genes(cursor)

    print("Adding novel transcripts to database...")
    batch_add_transcripts(cursor, struct_collection.transcript_dict, batch_size) 

    print("Adding novel exons/introns to database...")
    batch_add_edges(cursor, struct_collection.edge_dict, batch_size)

    print("Adding novel vertices/locations to database...")
    batch_add_locations(cursor, struct_collection.location_dict, batch_size)

    print("Updating gene-vertex assignments...")
    batch_add_vertex2gene(cursor, struct_collection.vertex_2_gene, batch_size)

    print("Adding %d dataset record(s) to database..." % len(datasets))
    add_datasets(cursor, datasets)  
 
    print("Updating abundance table....")
    batch_add_abundance(cursor, abundance, batch_size) 

    print("Adding %d transcript observation(s) to database..." % len(observed_transcripts))
    batch_add_observed(cursor, observed_transcripts, batch_size)

    print("Updating counters...")
    update_counter(cursor, struct_collection.run_info)

    print("Updating gene, transcript, and exon annotations...")
    batch_add_annotations(cursor, gene_annotations, "gene", batch_size)
    batch_add_annotations(cursor, transcript_annotations, "transcript", batch_size)
    batch_add_annotations(cursor, exon_annotations, "exon", batch_size)

    return
 
def update_counter(cursor, run_info):
    # Update the database counter
    
    update_g = 'UPDATE "counters" SET "count" = ? WHERE "category" = "genes"'
    cursor.execute(update_g,[run_info.genes])
    
    update_t = 'UPDATE "counters" SET "count" = ? WHERE "category" = "transcripts"'
    cursor.execute(update_t,[run_info.transcripts])

    update_e = 'UPDATE "counters" SET "count" = ? WHERE "category" = "edge"'
    cursor.execute(update_e,[run_info.edge])

    update_v = 'UPDATE "counters" SET "count" = ? WHERE "category" = "vertex"'
    cursor.execute(update_v,[run_info.vertex])

    update_d = 'UPDATE "counters" SET "count" = ? WHERE "category" = "dataset"'
    cursor.execute(update_d,[run_info.dataset])

    update_o = 'UPDATE "counters" SET "count" = ? WHERE "category" = "observed"'
    cursor.execute(update_o,[run_info.observed])

    return

def batch_add_vertex2gene(cursor, vertex_2_gene, batch_size):
    """ Add new vertex-gene relationships to the vertex table """
    entries = []
    for vertex_ID, gene_set in vertex_2_gene.items():
        for gene in gene_set:
            entries.append((vertex_ID, gene[0]))

    index = 0
    while index < len(entries):
        try:
            batch = entries[index:index + batch_size]
        except:
            batch = entries[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["vertex_ID", "gene_ID"]]) + ") "
            command = 'INSERT OR IGNORE INTO "vertex"' + cols + "VALUES " + \
                      '(?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)
            sys.exit(1)
    return

def batch_add_locations(cursor, location_dict, batch_size):
    """ Add new locations to database """
    location_entries = []
    for chrom_dict in location_dict.values():
        for loc in list(chrom_dict.values()):
            if type(loc) is dict:
                location_entries.append((loc['location_ID'],
                                         loc['genome_build'],
                                         loc['chromosome'],
                                         loc['position']))

    index = 0
    while index < len(location_entries):
        try:
            location_batch = location_entries[index:index + batch_size]
        except:
            location_batch = location_entries[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["location_ID", "genome_build", "chromosome", "position"]]) + ") "
            command = 'INSERT INTO "location"' + cols + "VALUES " + \
                      '(?,?,?,?)'
            cursor.executemany(command, location_batch)

        except Exception as e:
            print(e)
            sys.exit(1)
    return
    

def batch_add_edges(cursor, edge_dict, batch_size):
    """ Add new edges to database """
    edge_entries = []
    for edge in list(edge_dict.values()):
        if type(edge) is dict:
            edge_entries.append((edge['edge_ID'],
                                 edge['v1'],
                                 edge['v2'],
                                 edge['edge_type']))

    index = 0
    while index < len(edge_entries):
        try:
            batch = edge_entries[index:index + batch_size]
        except:
            batch = edge_entries[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["edge_ID", "v1", "v2", "edge_type"]]) + ") "
            command = 'INSERT INTO "edge"' + cols + "VALUES " + '(?,?,?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)
            sys.exit(1)

    return

def batch_add_transcripts(cursor, transcript_dict, batch_size):
    """ Add new transcripts to database """

    transcript_entries = []
    for transcript in list(transcript_dict.values()):
        if type(transcript) is dict:
            transcript_entries.append((transcript['transcript_ID'],
                                       transcript['gene_ID'],
                                       transcript['path'],
                                       transcript['start_vertex'],
                                       transcript['end_vertex'],
                                       transcript['n_exons']))

    index = 0
    while index < len(transcript_entries):
        try:
            transcript_batch = transcript_entries[index:index + batch_size]
        except:
            transcript_batch = transcript_entries[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["transcript_id", "gene_id", "path", "start_vertex",
                     "end_vertex", "n_exons"]]) + ") "
            command = 'INSERT INTO "transcripts"' + cols + "VALUES " + '(?,?,?,?,?,?)'
            cursor.executemany(command, transcript_batch)
        except Exception as e:
            print(e)
            sys.exit(1)

    return

def add_genes(cursor):
    """ Extract gene entries from the temporary table """

    query = "INSERT or IGNORE INTO genes SELECT gene_ID, strand FROM temp_gene;"
    cursor.execute(query)
    return

def add_datasets(cursor, datasets):
    """ Add dataset records to database """

    try:
        cols = " (" + ", ".join([str_wrap_double(x) for x in
               ["dataset_ID", "dataset_name", "sample", "platform"]]) + ") "
        command = 'INSERT INTO "dataset"' + cols + \
                  "VALUES " + '(?,?,?,?)'
        cursor.executemany(command, datasets)

    except Exception as e:
        print(e)
        sys.exit(1)
    return

def batch_add_annotations(cursor, annotations, annot_type, batch_size):
    """ Add gene/transcript/exon annotations to the appropriate annotation table
    """
    batch_size = 1
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
            command = 'INSERT OR IGNORE INTO "' + annot_type + '_annotations" ' + cols + \
                      "VALUES " + '(?,?,?,?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)
            sys.exit(1)
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
            sys.exit(1)
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
            sys.exit(1)
    return 

def check_database_integrity(cursor):
    """ Perform some checks on the database. Run before committing changes"""

    print("Validating database........")
    # For each category, check that the number of table entries matches the counter
    counter_query = "SELECT * FROM counters"
    cursor.execute(counter_query)
    counters = cursor.fetchall()

    for table_name, curr_counter in counters:
        curr_counter = int(curr_counter)

        # Vertex case needs to be handled differently
        if table_name == "vertex":
            table_name = "location"
             
        query = "select COUNT(*) from " + table_name
        cursor.execute(query)
        actual_count = int(cursor.fetchone()[0])

        if actual_count != curr_counter:
            print("table_count: "  + str(actual_count))
            print("counter_value: " + str(curr_counter))
            raise ValueError("Database counter for '" + table_name + \
                  "' does not match the number of entries in the table." + \
                  " Discarding changes to database and exiting...")

    return

def write_counts_log_file(cursor, outprefix):
    """ Create a log file with the following columns:
            - dataset name
            - Number of reads annotated
            - Number of known genes detected (total)
            - Number of novel genes detected (total)
            - Number of known transcripts detected (total)
            - Number of novel transcripts detected (total)
            Breakdowns by category
            - Number of antisense genes detected
            - Number of intergenic genes detected
            - Number of known transcripts
            - Number of FSM transcripts detected (perfect + with novelty)
            - Number of total ISM transcripts detected
            - Number of suffix ISMs detected
            - Number of antisense transcripts detected
            - Number of genomic transcripts detected
    """

    # Run utility
    summarize_datasets.write_counts_file(cursor, outprefix)

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

    # Read and annotate input sam files. Also, write output QC log file.
    print("Processing SAM files...")
    print("-------------------------------------")
    datasets, observed_transcripts, gene_annotations, transcript_annotations, \
    exon_annotations, abundance = process_all_sam_files(sam_files, dataset_list, cursor, 
                                      struct_collection, outprefix)
    print("-------------------------------------")

    # Update database
    batch_size = 10000
    update_database(cursor, batch_size, datasets, observed_transcripts,
                    gene_annotations, transcript_annotations, exon_annotations,
                    abundance, struct_collection)

    # Validate database
    check_database_integrity(cursor) 
    conn.commit()

    # TODO: output files
    # Write a file enumerating how many known/novel genes and transcripts
    # were detected in each dataset
    #write_counts_log_file(cursor, outprefix)
    conn.close()


if __name__ == '__main__':
    #pr = cProfile.Profile()
    #pr.enable()
    main()
    #pr.disable()
    #pr.print_stats()
