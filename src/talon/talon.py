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
from . import dstruct
import operator
#import os
from pathlib import Path
import warnings
from . import transcript_utils as tutils
from . import query_utils as qutils
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

def make_gene_start_and_end_dict(cursor, build):
    """ Format of dicts:
            Key: gene ID from database
            Value: dict mapping positions to start vertices (or end vertices) of 
                   KNOWN transcripts from that gene
    """
    gene_starts = {}
    gene_ends = {}
    query = """SELECT  gene_ID,
                               start_vertex,
                               end_vertex,
                               loc1.position as start,
                               loc2.position as end
               FROM transcripts
               LEFT JOIN transcript_annotations as ta
                   ON ta.ID = transcripts.transcript_ID
               LEFT JOIN location as loc1
                   ON transcripts.start_vertex = loc1.location_ID
               LEFT JOIN location as loc2
                   ON transcripts.end_vertex = loc2.location_ID
               WHERE ta.attribute = 'transcript_status'
                     AND ta.value = 'KNOWN'
                     AND loc1.genome_build = '%s'
                     AND loc2.genome_build = '%s'
              """
    cursor.execute(query % (build, build))
    for entry in cursor.fetchall():
        gene_ID = entry['gene_ID']
        start_vertex = entry['start_vertex']
        end_vertex = entry['end_vertex']
        start_pos = entry['start']
        end_pos = entry['end']
        try:
            gene_starts[gene_ID][start_pos] = start_vertex
        except:
            gene_starts[gene_ID] = {}
            gene_starts[gene_ID][start_pos] = start_vertex

        try:
            gene_ends[gene_ID][end_pos] = end_vertex
        except:
            gene_ends[gene_ID] = {}
            gene_ends[gene_ID][end_pos] = end_vertex

    return gene_starts, gene_ends
           

def make_transcript_dict(cursor, build):
    """ Format of dict:
            Key: tuple consisting of edges in transcript path
            Value: SQLite3 row from transcript table
    """
    transcript_dict = {}
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
        transcript_path = transcript["jn_path"]
        if transcript_path != None:
            transcript_path = transcript_path.split(",") + \
                              [transcript["start_exon"], transcript["end_exon"]]
            transcript_path = frozenset([ int(x) for x in transcript_path])
        else:
            transcript_path = frozenset([transcript["start_exon"]])
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

def make_temp_monoexonic_transcript_table(cursor, build):
    """ Attaches a temporary database with a table that has the following fields:
            - gene_ID
            - transcript_ID
            - chromosome
            - start (min position)
            - end (max position)
            - strand
        The purpose is to allow location-based matching for monoexonic query
        transcripts. """

    command = """ CREATE TEMPORARY TABLE IF NOT EXISTS temp_monoexon AS
                  SELECT t.gene_ID, 
                         t.transcript_ID, 
			 loc1.chromosome, 
			 loc1.position as start,
                         loc2.position as end,
			 genes.strand,
                         t.start_vertex,
                         t.end_vertex,
                         t.start_exon as exon_ID
                  FROM transcripts as t
	          LEFT JOIN location as loc1 ON loc1.location_ID = t.start_vertex
	          LEFT JOIN location as loc2 ON loc2.location_ID = t.end_vertex
	          LEFT JOIN genes ON genes.gene_ID = t.gene_ID
                  WHERE n_exons = 1 AND loc1.genome_build = '%s' 
	          AND loc2.genome_build = '%s' """
    cursor.execute(command % (build, build))
    
    return

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

def match_monoexon_vertices(chromosome, positions, strand, location_dict,
                                  run_info):
    """ Given the start and end of a single-exon transcript, this function looks 
        for a matching vertex for each position. Also returns a list where each 
        index indicates whether that vertex is novel to the data structure 
        (0 for known, 1 for novel) """

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

        if vertex_match == None:
            # If no vertex matches the position, one is created.
            vertex_match = create_vertex(chromosome, position, run_info,
                                         location_dict)["location_ID"]
            novelty.append(1)
        else:
            novelty.append(0)

        # Add to running list of matches
        vertex_matches.append(vertex_match)

    return tuple(vertex_matches), tuple(novelty), diff_5p, diff_3p

def match_splice_vertices(chromosome, positions, strand, location_dict,
                                  run_info):
    """ Given a chromosome and a list of positions from the transcript in 5' to
        3' end order, this function looks for a matching vertex for each splice
        junction position (so it ignores the ends). Also returns a list where 
        each index indicates whether that vertex is novel to the data structure 
        (0 for known, 1 for novel) """

    # Returned by function
    vertex_matches = []
    novelty = []
  
    # Iterate over positions
    for curr_index in range(1,len(positions)-1):
        position = positions[curr_index]

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

    return vertex_matches, novelty
   

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

def permissive_match_with_gene_priority(chromosome, position, strand, sj_pos, 
                                        pos_type, gene_ID, gene_locs, locations, run_info):
    """ Tries to match a position to a known start/end vertex from the same 
        gene. If none is found, the normal permissive match procedure is
        invoked.
    """
    # Check inputs
    if pos_type != "start" and pos_type != "end":
        raise ValueError("Please set pos_type to either 'start' or 'end'.")
    if strand != "+" and strand != "-":
        raise ValueError("Invalid strand specified: %s" % strand)

    # Try exact match first
    if chromosome in locations and position in locations[chromosome]:
        match = locations[chromosome][position]
        dist = 0
        if gene_ID in gene_locs:
            if position in gene_locs[gene_ID]:
                return match['location_ID'], dist, 1
            else:
                return match['location_ID'], dist, 0
        else:
            return match['location_ID'], dist, 0   
 
    # This approach only works when there are known starts/ends for this gene
    if gene_ID in gene_locs:

        # Get cutoff distance
        if pos_type == "start":
            max_dist = run_info.cutoff_5p
        else:
            max_dist = run_info.cutoff_3p

        if (strand == "+" and pos_type == "start") or \
           (strand == "-" and pos_type == "end"):
            search_window_start = position - max_dist
            search_window_end = sj_pos
        else:
            search_window_start = sj_pos
            search_window_end = position + max_dist

        # Compute distance to all of the gene positions on file
        min_abs_dist = max_dist + 1
        best_dist = None
        closest_vertex = None
        for known_location in gene_locs[gene_ID]:
            if known_location < search_window_start or known_location > search_window_end:
                continue
        
            curr_dist = compute_delta(known_location, position, strand)
            if abs(curr_dist) < min_abs_dist:
                best_dist = curr_dist
                min_abs_dist = abs(curr_dist)
                closest_vertex = gene_locs[gene_ID][known_location]

        # If a valid match is found, return it
        if min_abs_dist <= max_dist:
            return closest_vertex, best_dist, 1

    # Otherwise, revert to permissive match approach.
    match, dist = permissive_vertex_search(chromosome, position, strand, 
                                           sj_pos, pos_type,
                                           locations, run_info)
    return match, dist, 0
    

def permissive_vertex_search(chromosome, position, strand, sj_pos, pos_type,
                             locations, run_info):
    """ Given a position, this function tries to find a vertex match within the
        cutoff distance that also comes before the splice junction begins.
        If no vertex is found, the function returns None. """

    # Try a strict match first
    if chromosome in locations and position in locations[chromosome]:
        match = locations[chromosome][position]
        dist = 0
        return match['location_ID'], dist

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
                dist = compute_delta(curr_pos, position, strand)
                return match['location_ID'], dist

        curr_pos = position - dist*direction_priority
        if curr_pos > search_window_start and curr_pos < search_window_end:
            match = search_for_vertex_at_pos(chromosome, curr_pos, locations)
            if match != None:
                dist = compute_delta(curr_pos, position, strand)
                return match['location_ID'], dist

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

def create_transcript(chromosome, start_pos, end_pos, gene_ID, edge_IDs, vertex_IDs, 
                      transcript_dict, run_info):
    """Creates a novel transcript and adds it to the transcript data structure.
    """
    run_info.transcripts += 1
    new_ID = run_info.transcripts
    if len(edge_IDs) > 1:
        jn_path = ",".join(map(str, edge_IDs[1:-1]))
    else:
        jn_path = None
        
    new_transcript = {'transcript_ID': new_ID,
                      'gene_ID': gene_ID,
                      'jn_path': jn_path,
                      'start_exon': edge_IDs[0],
                      'end_exon': edge_IDs[-1],
                      'start_vertex': vertex_IDs[0],
                      'end_vertex': vertex_IDs[-1],
                      'n_exons': (len(edge_IDs) + 1)/2,
                      'chromosome': chromosome,
                      'start_pos': start_pos,
                      'end_pos': end_pos }

    path_key = frozenset(edge_IDs)
    transcript_dict[path_key] = new_transcript

    return new_transcript

def check_all_exons_known(novelty):
    """ Given a list in which each element represents the novelty (1) or
        known-ness of a transcript edge (0), determine whether all of the
        exons are known or not. Return True if all are known, and False
        otherwise. Input should not include first or last exon. """

    if len(novelty) == 1:
        return novelty[0] == 0

    exons = novelty[1::2]

    if sum(exons) > 0:
        return False
    else:
        return True

def check_all_SJs_known(novelty):
    """ Given a list in which each element represents the novelty (1) or 
        known-ness of a transcript edge (0), determine whether all of the
        introns are known or not. Return True if all are known, and False 
        otherwise. Input should not include first or last exon. If there is
        only one entry, then that means there is one splice junction (two exons)"""

    if len(novelty) == 1:
        return novelty[0] == 0

    introns = novelty[::2]
        
    if sum(introns) > 0:
        return False
    else:
        return True

def match_all_splice_edges(vertices, strand, edge_dict, run_info):
    """ Given a list of splice junction-only vertex IDs from the transcript in 5' to
        3' end order, this function looks for a matching edge ID for each
        position. If none exists, it creates one. """

    edge_matches = []
    novelty = []
    edge_type = "intron"

    for index_1 in range(0, len(vertices) - 1):
        index_2 = index_1 + 1

        if index_1 % 2 != 0:
            edge_type = "exon"
        else:
            edge_type = "intron"

        vertex_1 = vertices[index_1]
        vertex_2 = vertices[index_2]

        edge_match, curr_novelty = match_or_create_edge(vertex_1, vertex_2, 
                                                        edge_type, strand, 
                                                        edge_dict, run_info)
        edge_matches.append(edge_match)
        novelty.append(curr_novelty)

    return edge_matches, novelty

def match_or_create_edge(vertex_1, vertex_2, edge_type, strand, edge_dict, run_info):
    """ Searches for edge match to provided set of vertices. If none found, 
        creates a new edge. """
    novelty = 0
    edge_match = search_for_edge(vertex_1, vertex_2, edge_type, edge_dict)

    if edge_match == None:
        # If no edge matches the position, one is created.
        edge_match = create_edge(vertex_1, vertex_2, edge_type, strand,
                                     edge_dict, run_info)
        novelty = 1
    return edge_match["edge_ID"], novelty

def match_all_transcript_edges(vertices, strand, edge_dict, run_info):
    """ Given a list of vertex IDs from the transcript in 5' to
        3' end order, this function looks for a matching edge ID for each
        position. If none exists, it creates one. Only used for monoexon case"""

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
        
        edge_match, curr_novelty = match_or_create_edge(vertex_1, vertex_2, 
                                                        edge_type, strand,
                                                        edge_dict, run_info) 
        edge_matches.append(edge_match)
        novelty.append(curr_novelty)
                                                
    return tuple(edge_matches), tuple(novelty)


def search_for_ISM(edge_IDs, transcript_dict):
    """ Given a list of edges in a query transcript, determine whether it is an
        incomplete splice match (ISM) of any transcript in the dict. Will also
        return FSM matches if they're there"""               

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

def process_FSM(chrom, positions, strand, edge_IDs, vertex_IDs, all_matches, gene_starts, gene_ends,
                edge_dict, locations, run_info):
    """ Given a transcript, try to find an FSM match for it """
    gene_ID = None
    transcript_ID = None
    novelty = []
    start_vertex = None
    end_vertex = None
    start_exon = None
    end_exon = None
    diff_5p = None 
    diff_3p = None   
    jn_path = None

    # Check if any of the matches have the same number of exons as the query.
    # Such a match should be prioritized because it's an FSM
    n_exons = len(positions)/2
    FSM_matches = [ x for x in all_matches if x['n_exons'] == n_exons ]

    if len(FSM_matches) == 0:
        return None, None, [], None

    else:
        transcript_match = FSM_matches[0]
        gene_ID = transcript_match['gene_ID']
        transcript_ID = transcript_match['transcript_ID']
 
        # Check whether the query's 5' and 3' ends are within range of those of
        # the match. If not, perform a permissive match assignment
        curr_5p_diff = compute_delta(transcript_match['start_pos'], positions[0],
                                     strand)
        curr_3p_diff = compute_delta(transcript_match['end_pos'], positions[-1],
                                     strand) 
        # ---------------- 5' end ---------------------------------
        if abs(curr_5p_diff) <= run_info.cutoff_5p:
            start_vertex = transcript_match['start_vertex']
            start_exon = transcript_match['start_exon']
            diff_5p = curr_5p_diff
            start_novelty = 0
        else:
            # First get a permissively matched start vertex
            start_vertex, start_exon, start_novelty, known_start, diff_5p = process_5p(chrom, 
                                                                   positions, strand,
                                                                   vertex_IDs, 
                                                                   gene_ID, gene_starts, 
                                                                   edge_dict,
                                                                   locations, run_info)
        # ---------------- 3' end ---------------------------------
        if abs(curr_3p_diff) <= run_info.cutoff_3p:
            end_vertex = transcript_match['end_vertex']
            end_exon = transcript_match['end_exon']
            diff_3p = curr_3p_diff
            end_novelty = 0
        else:
            # First get a permissively matched end vertex
            end_vertex, end_exon, end_novelty, known_end, diff_3p = process_3p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_ends,
                                                                   edge_dict,
                                                                   locations, run_info)

        edge_IDs = [start_exon] + edge_IDs + [end_exon]
        vertex_IDs = [start_vertex] + vertex_IDs + ["end_vertex"]

    # Package information for output
    start_end_info = {"start_vertex": start_vertex,
                      "end_vertex": end_vertex,
                      "start_exon": start_exon,
                      "end_exon": end_exon,
                      "diff_5p": diff_5p,
                      "diff_3p": diff_3p,
                      "start_novelty": start_novelty,
                      "end_novelty": end_novelty,
                      "vertex_IDs": vertex_IDs,
                      "edge_IDs": edge_IDs}

    return gene_ID, transcript_ID, novelty, start_end_info

def process_5p(chrom, positions, strand, vertex_IDs, gene_ID, gene_starts, edge_dict,
               locations, run_info):
    """ Conduct permissive match for 5' end and return assigned vertex,
        edge, and distance """

    # First get a permissively matched start vertex
    start_vertex, diff_5p, known_start = permissive_match_with_gene_priority(chrom,
                                         positions[0], strand, positions[1],
                                         "start", gene_ID, gene_starts,
                                         locations, run_info)
    if start_vertex == None:
        start_vertex = create_vertex(chrom, positions[0], run_info,
                                             locations)['location_ID']

    # Then get the start exon
    start_exon, start_novelty = match_or_create_edge(start_vertex,
                                                     vertex_IDs[0],
                                                     "exon", strand,
                                                     edge_dict, run_info)

    # If known_start == 1, the start vertex is a known startpoint of this gene.
    #  start novelty refers to the novelty of the first exon (1 if yes, 0 if no)
    return start_vertex, start_exon, start_novelty, known_start, diff_5p


def process_3p(chrom, positions, strand, vertex_IDs, gene_ID, gene_ends, edge_dict,
               locations, run_info):
    """ Conduct permissive match for 3' end and return assigned vertex,
        edge, and distance """

    # First get a permissively matched end vertex
    end_vertex, diff_3p, known_end = permissive_match_with_gene_priority(chrom,
                                          positions[-1], strand, positions[-2],
                                          "end", gene_ID, gene_ends,
                                          locations, run_info)
    if end_vertex == None:
        end_vertex = create_vertex(chrom, positions[-1], run_info,
                                   locations)['location_ID']
    # Then get the end exon
    end_exon, end_novelty = match_or_create_edge(vertex_IDs[-1],
                                                 end_vertex,
                                                 "exon", strand,
                                                  edge_dict, run_info)
    # If known_end == 1, the end vertex is a known endpoint of this gene.
    # end novelty refers to the novelty of the final exon (1 if yes, 0 if no)
    return end_vertex, end_exon, end_novelty, known_end, diff_3p


def process_ISM(chrom, positions, strand, edge_IDs, vertex_IDs, all_matches, transcript_dict,
                gene_starts, gene_ends, edge_dict, locations, run_info):
    """ Given a transcript, try to find an ISM match for it. If the
        best match is an ISM with known ends, that will be promoted to NIC. """

    gene_ID = None
    transcript_ID = None
    novelty = []
    start_end_info = {}
    n_exons = len(positions)/2

    ISM = []
    suffix = []
    prefix = []
    gene_ID = all_matches[0]['gene_ID']
    # Get matches for the ends
    if n_exons > 1:
        start_vertex, start_exon, start_novelty, known_start, diff_5p = process_5p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_starts,
                                                                   edge_dict,
                                                                   locations, run_info) 
        end_vertex, end_exon, end_novelty, known_end, diff_3p = process_3p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_ends,
                                                                   edge_dict,
                                                                   locations, run_info)
        # Update info
        edge_IDs = [start_exon] + edge_IDs + [end_exon]
        vertex_IDs = [start_vertex] + vertex_IDs + [end_vertex]
        
        start_end_info["start_vertex"] = start_vertex
        start_end_info["end_vertex"] = end_vertex
        start_end_info["start_exon"] = start_exon
        start_end_info["end_exon"] = end_exon
        start_end_info["start_novelty"] = start_novelty
        start_end_info["end_novelty"] = end_novelty
        start_end_info["diff_5p"] = diff_5p
        start_end_info["diff_3p"] = diff_3p
        start_end_info["edge_IDs"] = edge_IDs
        start_end_info["vertex_IDs"] = vertex_IDs
    else:
        known_start = 0
        known_end = 0   

    # If the 5' and 3' vertex sites are known for this gene, return NIC
    if known_start and known_end:
        novel_transcript = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)
        transcript_ID = novel_transcript['transcript_ID']
        novelty = [(transcript_ID, run_info.idprefix, "TALON",
                        "NIC_transcript", "TRUE")]
        return gene_ID, transcript_ID, novelty, start_end_info    

    # Iterate over matches to characterize ISMs
    for match in all_matches:
        transcript_ID = run_info.transcripts + 1

        # Add ISM
        ISM.append(str(match['transcript_ID']))

        # Single-exon case
        if n_exons == 1:
            if match["n_exons"] == 1:
                # Return FSM
                gene_ID = match["gene_ID"]
                transcript_ID = match["transcript_ID"]
                novelty = []
                return gene_ID, transcript_ID, novelty, start_end_info

            match_path = match['jn_path']
            exon = str(edge_IDs[0])
            # Look for prefix
            if match_path.startswith(exon):
                prefix.append(str(match['transcript_ID']))
            # Look for suffix
            if match_path.endswith(exon):
                suffix.append(str(match['transcript_ID']))
                gene_ID = match['gene_ID']
            continue

        # Multi-exon case
        edge_str = ",".join([str(x) for x in edge_IDs[1:-1]])

        # Look for prefix
        if match['jn_path'].startswith(edge_str):
            prefix.append(str(match['transcript_ID']))

        # Look for suffix
        if match['jn_path'].endswith(edge_str):
            gene_ID = match['gene_ID']
            suffix.append(str(match['transcript_ID'])) 

    novel_transcript = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)

    transcript_ID = novel_transcript['transcript_ID']

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

    return gene_ID, transcript_ID, novelty, start_end_info

    
def process_NIC(chrom, positions, strand, edge_IDs, vertex_IDs, transcript_dict,
                gene_starts, gene_ends, edge_dict, locations, vertex_2_gene, run_info):
    """ For a transcript that has been determined to be novel in catalog, find
        the proper gene match (documenting fusion event if applicable). To do 
        this, look up each vertex in the vertex_2_gene dict, and keep track of all
        same-strand genes. """

    gene_matches = []
    start_end_info = {}
    
    for vertex in vertex_IDs:
        if vertex in vertex_2_gene:
            curr_matches = vertex_2_gene[vertex]

            # Make sure the gene is on the correct strand
            gene_matches += [ x[0] for x in list(curr_matches) if x[1] == strand ]

  
    # Now count up how often we see each gene
    gene_tally = dict((x,gene_matches.count(x)) for x in set(gene_matches))

    # TODO: deal with fusions

    # For the main assignment, pick the gene that is observed the most
    if len(gene_tally) == 0:
        return None, None, [], None

    gene_ID = max(gene_tally, key=gene_tally.get)

    # Get matches for the ends
    start_vertex, start_exon, start_novelty, known_start, diff_5p = process_5p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_starts,
                                                                   edge_dict,
                                                                   locations, run_info)
    end_vertex, end_exon, end_novelty, known_end, diff_3p = process_3p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_ends,
                                                                   edge_dict,
                                                                   locations, run_info)
    # Update info
    edge_IDs = [start_exon] + edge_IDs + [end_exon]
    vertex_IDs = [start_vertex] + vertex_IDs + [end_vertex]
    start_end_info["start_vertex"] = start_vertex
    start_end_info["end_vertex"] = end_vertex
    start_end_info["start_exon"] = start_exon
    start_end_info["end_exon"] = end_exon
    start_end_info["start_novelty"] = start_novelty
    start_end_info["end_novelty"] = end_novelty
    start_end_info["diff_5p"] = diff_5p
    start_end_info["diff_3p"] = diff_3p
    start_end_info["edge_IDs"] = edge_IDs
    start_end_info["vertex_IDs"] = vertex_IDs

    # Create a new transcript of that gene
    novel_transcript = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)
    transcript_ID = novel_transcript["transcript_ID"]
    novelty = [(transcript_ID, run_info.idprefix, "TALON",
                         "NIC_transcript","TRUE")]

    return gene_ID, transcript_ID, novelty, start_end_info    


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

def process_NNC(chrom, positions, strand, edge_IDs, vertex_IDs, transcript_dict,
                gene_starts, gene_ends, edge_dict, locations, vertex_2_gene, run_info):
    """ Novel not in catalog case """

    novelty = []
    start_end_info = {}

    gene_ID = find_gene_match_on_vertex_basis(vertex_IDs, strand, vertex_2_gene)
    if gene_ID == None:
        return None, None, [], None

    # Get matches for the ends
    start_vertex, start_exon, start_novelty, known_start, diff_5p = process_5p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_starts,
                                                                   edge_dict,
                                                                   locations, run_info)
    end_vertex, end_exon, end_novelty, known_end, diff_3p = process_3p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_ends,
                                                                   edge_dict,
                                                                   locations, run_info)
    # Update info
    edge_IDs = [start_exon] + edge_IDs + [end_exon]
    vertex_IDs = [start_vertex] + vertex_IDs + [end_vertex]
    start_end_info["start_vertex"] = start_vertex
    start_end_info["end_vertex"] = end_vertex
    start_end_info["start_exon"] = start_exon
    start_end_info["end_exon"] = end_exon
    start_end_info["start_novelty"] = start_novelty
    start_end_info["end_novelty"] = end_novelty
    start_end_info["diff_5p"] = diff_5p
    start_end_info["diff_3p"] = diff_3p
    start_end_info["edge_IDs"] = edge_IDs
    start_end_info["vertex_IDs"] = vertex_IDs
    
    transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]

    novelty.append((transcript_ID, run_info.idprefix, "TALON",
                               "NNC_transcript", "TRUE"))

    return gene_ID, transcript_ID, novelty, start_end_info

def process_spliced_antisense(chrom, positions, strand, edge_IDs, vertex_IDs, transcript_dict,
                gene_starts, gene_ends, edge_dict, locations, vertex_2_gene, run_info, cursor):
    """ Annotate a transcript as antisense with splice junctions """

    gene_novelty = []
    transcript_novelty = []
    start_end_info = {}

    if strand == "+":
        anti_strand = "-"
    else:
        anti_strand = "+"
    anti_gene_ID = find_gene_match_on_vertex_basis(vertex_IDs, anti_strand,
                                                       vertex_2_gene)
    if anti_gene_ID == None:
        return None, None, gene_novelty, transcript_novelty, start_end_info

    # Take care of ends
    start_vertex, start_exon, start_novelty, known_start, diff_5p = process_5p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   anti_gene_ID, gene_ends,
                                                                   edge_dict,
                                                                   locations, run_info)
    end_vertex, end_exon, end_novelty, known_end, diff_3p = process_3p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   anti_gene_ID, gene_starts,
                                                                   edge_dict,
                                                                   locations, run_info)
    # Update info
    edge_IDs = [start_exon] + edge_IDs + [end_exon]
    vertex_IDs = [start_vertex] + vertex_IDs + [end_vertex]
    start_end_info["start_vertex"] = start_vertex
    start_end_info["end_vertex"] = end_vertex
    start_end_info["start_exon"] = start_exon
    start_end_info["end_exon"] = end_exon
    start_end_info["start_novelty"] = start_novelty
    start_end_info["end_novelty"] = end_novelty
    start_end_info["diff_5p"] = diff_5p
    start_end_info["diff_3p"] = diff_3p
    start_end_info["edge_IDs"] = edge_IDs
    start_end_info["vertex_IDs"] = vertex_IDs    

    gene_ID = create_gene(chrom, positions[0], positions[-1],
                              strand, cursor, run_info)
    transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]

    # Handle gene annotations
    gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                             "antisense_gene","TRUE"))
    gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_antisense_to_IDs",anti_gene_ID))

    # Handle transcript annotations
    transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                   "antisense_transcript", "TRUE"))

    return gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info

def process_remaining_mult_cases(chrom, positions, strand, edge_IDs, vertex_IDs, transcript_dict,
                gene_starts, gene_ends, edge_dict, locations, vertex_2_gene, run_info, cursor):
    """ This function is a catch-all for multiexonic transcripts that were not
        FSM, ISM, NIC, NNC, or spliced antisense.
    """
    gene_novelty = []
    transcript_novelty = []
    start_end_info = {}

    gene_ID, match_strand = search_for_overlap_with_gene(chrom, positions[0],
                                                         positions[1], strand,
                                                         cursor, run_info)

    # We don't care about the gene when making these assignments
    start_vertex, start_exon, start_novelty, known_start, diff_5p = process_5p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_starts,
                                                                   edge_dict,
                                                                   locations, run_info)
    end_vertex, end_exon, end_novelty, known_end, diff_3p = process_3p(chrom,
                                                                   positions, strand,
                                                                   vertex_IDs,
                                                                   gene_ID, gene_ends,
                                                                   edge_dict,
                                                                   locations, run_info)
    # Update info
    edge_IDs = [start_exon] + edge_IDs + [end_exon]
    vertex_IDs = [start_vertex] + vertex_IDs + [end_vertex]
    start_end_info["start_vertex"] = start_vertex
    start_end_info["end_vertex"] = end_vertex
    start_end_info["start_exon"] = start_exon
    start_end_info["end_exon"] = end_exon
    start_end_info["start_novelty"] = start_novelty
    start_end_info["end_novelty"] = end_novelty
    start_end_info["diff_5p"] = diff_5p
    start_end_info["diff_3p"] = diff_3p
    start_end_info["edge_IDs"] = edge_IDs
    start_end_info["vertex_IDs"] = vertex_IDs

    if gene_ID == None:
        gene_ID = create_gene(chrom, positions[0], positions[-1],
                          strand, cursor, run_info)

        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                     "intergenic_novel","TRUE"))

        transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                              "intergenic_transcript", "TRUE"))

    elif match_strand != strand:
        anti_gene_ID = gene_ID
        gene_ID = create_gene(chrom, positions[0], positions[-1],
                          strand, cursor, run_info)
        transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]

        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                     "antisense_gene","TRUE"))
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                     "gene_antisense_to_IDs",anti_gene_ID))
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                              "antisense_transcript", "TRUE"))
    else:
        transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                              "genomic_transcript", "TRUE"))
    
    return gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info

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
    gene_ID = None
 
    # Get vertex matches for the transcript positions
    vertex_IDs, v_novelty = match_splice_vertices(chrom, positions, strand,
                                                   location_dict, run_info)

    # Get edge matches for transcript exons and introns based on the vertices
    edge_IDs, e_novelty = match_all_splice_edges(vertex_IDs, strand,
                                                     edge_dict, run_info)

    # Check novelty of exons and splice jns. This will help us categorize 
    # what type of novelty the transcript has
    all_SJs_known = check_all_SJs_known(e_novelty)
    all_exons_known = check_all_exons_known(e_novelty)
    splice_vertices_known = (sum(v_novelty) == 0)
    all_exons_novel = (reduce(operator.mul, e_novelty, 1) == 1)

    # Look for FSM or ISM. 
    if all_SJs_known:
        # Get all FSM/ISM matches
        all_matches = search_for_ISM(edge_IDs, transcript_dict)
        if all_matches != None:
            # Look for FSM first
            gene_ID, transcript_ID, transcript_novelty, start_end_info = process_FSM(chrom,
                                                            positions, strand, edge_IDs,
                                                            vertex_IDs, all_matches,
                                                            gene_starts, gene_ends,
                                                            edge_dict,
                                                            location_dict, run_info)
            if gene_ID == None:
                # Now look for ISM
                gene_ID, transcript_ID, transcript_novelty, start_end_info = process_ISM(chrom,
                                                            positions,
                                                            strand, edge_IDs,
                                                            vertex_IDs,
                                                            all_matches,
                                                            transcript_dict,
                                                            gene_starts, gene_ends,
                                                            edge_dict, location_dict,
                                                            run_info)
        # Look for NIC
        if gene_ID == None:
            gene_ID, transcript_ID, transcript_novelty, start_end_info = process_NIC(chrom,
                                                            positions,
                                                            strand, edge_IDs,
                                                            vertex_IDs, transcript_dict,
                                                            gene_starts, gene_ends,
                                                            edge_dict, location_dict,
                                                            vertex_2_gene, run_info)
            
    # Novel in catalog transcripts have known splice donors and acceptors,
    # but new connections between them. 
    elif splice_vertices_known and gene_ID == None:
        gene_ID, transcript_ID, transcript_novelty, start_end_info = process_NIC(chrom,
                                                            positions,
                                                            strand, edge_IDs,
                                                            vertex_IDs, transcript_dict,
                                                            gene_starts, gene_ends,
                                                            edge_dict, location_dict,
                                                            vertex_2_gene, run_info)
    
    # Antisense transcript with splice junctions matching known gene
    if splice_vertices_known and gene_ID == None:
        gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info = \
                                      process_spliced_antisense(chrom, positions,
                                                                  strand, edge_IDs,
                                                                  vertex_IDs,
                                                                  transcript_dict,
                                                                  gene_starts,
                                                                  gene_ends,
                                                                  edge_dict, location_dict,
                                                                  vertex_2_gene, run_info,
                                                                  cursor)

    # Novel not in catalog transcripts contain new splice donors/acceptors
    # and contain at least one splice junction.
    elif not(splice_vertices_known): 
        gene_ID, transcript_ID, transcript_novelty, start_end_info = process_NNC(chrom,
                                                            positions,
                                                            strand, edge_IDs,
                                                            vertex_IDs, transcript_dict,
                                                            gene_starts, gene_ends,
                                                            edge_dict, location_dict,
                                                            vertex_2_gene, run_info)

    # Transcripts that don't match the previous categories end up here
    if gene_ID == None:
        gene_ID, transcript_ID, gene_novelty, transcript_novelty, start_end_info = \
                             process_remaining_mult_cases(chrom, positions,
                                                                strand, edge_IDs,
                                                                vertex_IDs,
                                                                transcript_dict,
                                                                gene_starts, gene_ends,
                                                                edge_dict, location_dict,
                                                                vertex_2_gene, run_info,
                                                                cursor)

    # Add all novel vertices to vertex_2_gene now that we have the gene ID
    vertex_IDs = start_end_info["vertex_IDs"]
    edge_IDs = start_end_info["edge_IDs"]
    e_novelty = [start_end_info["start_novelty"]] + e_novelty + \
                [start_end_info["end_novelty"]]

    update_vertex_2_gene(gene_ID, vertex_IDs, strand, vertex_2_gene)
 
    # For novel genes and transcripts, add names to novelty entries
    talon_gene_name, talon_transcript_name = construct_names(gene_ID, 
                                                             transcript_ID,
                                                             run_info.idprefix,
                                                             run_info.n_places)
    if len(gene_novelty) > 0:
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                                   "gene_status", "NOVEL"))
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_name", talon_gene_name))
        gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_id", talon_gene_name))
    if len(transcript_novelty) > 0:
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                   "transcript_status", "NOVEL"))
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_name", talon_transcript_name))
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_id", talon_transcript_name))

    # Add annotation entries for any novel exons
    exon_novelty = []
    exons = edge_IDs[::2]
    e_novelty = e_novelty[::2]

    if sum(e_novelty) > 0:
        for exon,is_novel in zip(exons, e_novelty):
            if is_novel:
                exon_novelty.append((exon, run_info.idprefix, "TALON", 
                                     "exon_status", "NOVEL"))

    # Package up information for output
    annotations = {'gene_ID': gene_ID,
                   'transcript_ID': transcript_ID,
                   'gene_novelty': gene_novelty,
                   'transcript_novelty': transcript_novelty,
                   'exon_novelty': exon_novelty,
                   'start_vertex': start_end_info["start_vertex"],
                   'end_vertex': start_end_info["end_vertex"],
                   'start_exon': start_end_info["start_exon"],
                   'end_exon': start_end_info["end_exon"],
                   'start_delta': start_end_info["diff_5p"],
                   'end_delta': start_end_info["diff_3p"]}
    
    return annotations

def construct_names(gene_ID, transcript_ID, prefix, n_places):
    """ Create a gene and transcript name using the TALON IDs.
        The n_places variable indicates how many characters long the numeric
        part of the name should be. """

    gene_ID_str = str(gene_ID).zfill(n_places)
    gene_name = prefix + "G" + gene_ID_str

    transcript_ID_str = str(transcript_ID).zfill(n_places)
    transcript_name = prefix + "T" + transcript_ID_str

    return gene_name, transcript_name

def check_inputs(options):
    """ Checks the input options provided by the user and makes sure that
        they are valid. Throw an error with descriptive help message if not."""

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
    make_temp_monoexonic_transcript_table(cursor, build)
    run_info = init_run_info(cursor, build, min_coverage, min_identity)
    location_dict = make_location_dict(build, cursor)
    edge_dict = make_edge_dict(cursor)
    transcript_dict = make_transcript_dict(cursor, build)
    vertex_2_gene = make_vertex_2_gene_dict(cursor)
    gene_starts, gene_ends = make_gene_start_and_end_dict(cursor, build)   

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


def compute_delta(orig_pos, new_pos, strand):
    """ Given a starting position and a new position, compute the distance
        between them. The sign indicates whether the second point is 
        upstream or downstream of the original with respect to strand. """

    abs_dist = abs(orig_pos - new_pos)
    if strand == "+":
        if new_pos < orig_pos:
            return -1*abs_dist
        else:
            return abs_dist
    elif strand == "-":
        if new_pos < orig_pos:
            return abs_dist
        else:
            return -1*abs_dist
    else:
        raise ValueError("Strand must be either + or -")
    

def identify_monoexon_transcript(chrom, positions, strand, cursor, location_dict,
                                 edge_dict, transcript_dict, vertex_2_gene,
                                 gene_starts, gene_ends, run_info):

    gene_novelty = []
    transcript_novelty = []
    exon_novelty = []

    cutoff_5p = run_info.cutoff_5p
    cutoff_3p = run_info.cutoff_3p

    start = positions[0]
    end = positions[-1]
    # First, look for a monoexonic FSM transcript match that meets the cutoff
    # distance criteria
    query = """ SELECT * 
                    FROM temp_monoexon AS tm
                    WHERE tm.chromosome = '%s'
                    AND tm.strand = '%s'
                    AND tm.start >= %d
                    AND tm.start <= %d
                    AND tm.end >= %d
                    AND tm.end <= %d  
               """
    lower_start_bound = start - cutoff_5p
    upper_start_bound = start + cutoff_5p
    lower_end_bound = end - cutoff_3p
    upper_end_bound = end + cutoff_3p
    cursor.execute(query % (chrom, strand, lower_start_bound, upper_start_bound,
                            lower_end_bound, upper_end_bound))
    matches = cursor.fetchall()

    # If there is more than one match, apply a tiebreaker (pick the one with 
    # the most overlap
    if len(matches) > 0:
        best_overlap = 0
        best_match = None
        for match in matches:
            # get overlap and compare
            match_interval = [match['start'], match['end']]
            overlap = get_overlap([start, end], match_interval)
            if overlap >= best_overlap:
                best_overlap = overlap
                best_match = match

        gene_ID = best_match['gene_ID']
        transcript_ID = best_match['transcript_ID']
        vertex_IDs = (best_match['start_vertex'], best_match['end_vertex'])
        edge_IDs = [best_match['exon_ID']]
        diff_5p = compute_delta(best_match['start'], start, strand)
        diff_3p = compute_delta(best_match['end'], end, strand)

    # If there is no match, proceed to genomic/antisense style matching.
    else:
        # Start by performing vertex match               
        vertex_IDs, v_novelty, diff_5p, diff_3p = match_monoexon_vertices(
                                                                 chrom,
                                                                 positions,
                                                                 strand,
                                                                 location_dict,
                                                                 run_info)

        # Get edge match (or create new edge)
        edge_IDs, e_novelty = match_all_transcript_edges(vertex_IDs, strand,
                                                     edge_dict, run_info)

        # If the exon is known, then this transcript must be ISM or NIC
        gene_ID = None
        if e_novelty[0] == 0:
            all_matches = search_for_ISM(edge_IDs, transcript_dict)

            if all_matches != None:
                gene_ID, transcript_ID, transcript_novelty, info = process_ISM(chrom, positions, 
                                                                     strand, edge_IDs, 
                                                                     vertex_IDs, all_matches, 
                                                                     transcript_dict,
                                                                     gene_starts, gene_ends, 
                                                                     edge_dict, location_dict, 
                                                                     run_info)
        if gene_ID == None:
            # Find best gene match using overlap search if the ISM/NIC check didn't work
            gene_ID, match_strand = search_for_overlap_with_gene(chrom, positions[0],
                                                             positions[1], strand,
                                                             cursor, run_info) 
            # Intergenic case
            if gene_ID == None:
                gene_ID = create_gene(chrom, positions[0], positions[-1],
                                      strand, cursor, run_info)

                gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                                     "intergenic_novel","TRUE"))
                transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]
                transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                      "intergenic_transcript", "TRUE"))     
            # Antisense case
            elif match_strand != strand:
                anti_gene_ID = gene_ID
                gene_ID = create_gene(chrom, positions[0], positions[-1],
                                      strand, cursor, run_info)
                transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]

                gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                                     "antisense_gene","TRUE"))
                gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                                        "gene_antisense_to_IDs",anti_gene_ID))
                transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                          "antisense_transcript", "TRUE"))

            # Same strand
            else:
                transcript_ID = create_transcript(chrom, positions[0], positions[-1],
                                              gene_ID, edge_IDs, vertex_IDs,
                                              transcript_dict, run_info)["transcript_ID"]
                transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                  "genomic_transcript", "TRUE"))

        # Add all novel vertices to vertex_2_gene now that we have the gene ID
        update_vertex_2_gene(gene_ID, vertex_IDs, strand, vertex_2_gene)

        talon_gene_name, talon_transcript_name = construct_names(gene_ID,
                                                             transcript_ID,
                                                             run_info.idprefix,
                                                             run_info.n_places)

        # Add novel gene annotation attributes
        if len(gene_novelty) > 0:
            gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_status", "NOVEL"))
            gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_name", talon_gene_name))
            gene_novelty.append((gene_ID, run_info.idprefix, "TALON",
                         "gene_id", talon_gene_name))

        # Add novel transcript annotation attributes
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                                   "transcript_status", "NOVEL"))
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_name", talon_transcript_name))
        transcript_novelty.append((transcript_ID, run_info.idprefix, "TALON",
                         "transcript_id", talon_transcript_name))

        # Add annotation entries for any novel exons
        if e_novelty[0] == 1:
            exon_novelty.append((edge_IDs[0], run_info.idprefix, "TALON",
                                     "exon_status", "NOVEL"))
        
        # Add the novel transcript to the temporary monoexon table
        new_mono = ( gene_ID, transcript_ID, chrom, start, end, strand, 
                     vertex_IDs[0], vertex_IDs[-1], edge_IDs[0] )
        cols = '("gene_ID", "transcript_ID", "chromosome", "start", "end",' + \
                 '"strand", "start_vertex", "end_vertex", "exon_ID")'
        command = 'INSERT INTO temp_monoexon ' + cols + ' VALUES ' + '(?,?,?,?,?,?,?,?,?)'
        cursor.execute(command, new_mono)

    # Package annotation information
    annotations = {'gene_ID': gene_ID,
                   'transcript_ID': transcript_ID,
                   'gene_novelty': gene_novelty,
                   'transcript_novelty': transcript_novelty,
                   'exon_novelty': exon_novelty,
                   'start_vertex': vertex_IDs[0],
                   'end_vertex': vertex_IDs[-1],
                   'start_exon': edge_IDs[0],
                   'end_exon': edge_IDs[-1],
                   'start_delta': diff_5p,
                   'end_delta': diff_3p}
    return annotations

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
            n_exons = len(positions)/2
            if n_exons > 1:
                annotation_info = identify_transcript(chrom, positions, strand, 
                                     cursor, location_dict, 
                                     edge_dict, transcript_dict, 
                                     vertex_2_gene, 
                                     gene_starts, gene_ends,
                                     run_info)
            else:
                annotation_info = identify_monoexon_transcript(chrom, positions, strand,
                                                  cursor, location_dict,
                                                  edge_dict, transcript_dict,
                                                  vertex_2_gene,
                                                  gene_starts, gene_ends,
                                                  run_info)
            #except Exception as e:
            #    print(e)
            #    warnings.warn("Problem identifying transcript '%s'. Skipping.."\
            #                   % read_ID) 
            #    sys.exit(1)
                #continue
                            
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
            start_exon = annotation_info['start_exon']
            end_exon = annotation_info['end_exon']

            struct_collection.run_info['observed'] += 1
            obs_ID = struct_collection.run_info['observed'] 
            observed = (obs_ID, gene_ID, transcript_ID, read_ID, dataset, 
                        start_vertex, end_vertex, start_exon, end_exon,
                        start_delta, end_delta, read_length)
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
                                 edge['edge_type'],
                                 edge['strand']))

    index = 0
    while index < len(edge_entries):
        try:
            batch = edge_entries[index:index + batch_size]
        except:
            batch = edge_entries[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["edge_ID", "v1", "v2", "edge_type", "strand"]]) + ") "
            command = 'INSERT INTO "edge"' + cols + "VALUES " + '(?,?,?,?,?)'
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
                                       transcript['start_exon'],
                                       transcript['jn_path'],
                                       transcript['end_exon'],
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
                   ["transcript_id", "gene_id", "start_exon", "jn_path", 
                     "end_exon", "start_vertex", "end_vertex", "n_exons"]]) + ") "
            command = 'INSERT INTO "transcripts"' + cols + "VALUES " + '(?,?,?,?,?,?,?,?)'
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
        dataset, start_vertex_ID, end_vertex_ID, start_exon, end_exon, 
        start_delta, end_delta, read_length) to observed table of database. """

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
                    "dataset", "start_vertex", "end_vertex",
                    "start_exon", "end_exon",
                    "start_delta", "end_delta", "read_length"]]) + ") "
            command = 'INSERT INTO "observed"' + cols + \
                      "VALUES " + '(?,?,?,?,?,?,?,?,?,?,?,?)'
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
