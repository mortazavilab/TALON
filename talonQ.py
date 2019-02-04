# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

import argparse
import sqlite3
import dstruct

def get_args():
    """ Fetches the arguments for the program """

    program_desc = """TALON takes transcripts from one or more long read
                      datasets (SAM format) and assigns them transcript and gene 
                      identifiers based on a database-bound annotation. 
                      Novel events are assigned new identifiers."""
    parser = argparse.ArgumentParser(description=program_desc)
    parser.add_argument('--db', dest = 'database', metavar='FILE,', type=str,
               help='TALON database. Created using build_talon_annotation.py')
    parser.add_argument('--build', dest = 'build', metavar='STRING,', type=str,
               help='Genome build (i.e. hg38) to use. Must be in the database.')
    parser.add_argument('--idprefix', dest = 'idprefix', metavar='STRING,', 
               type=str, default = "TALON",
               help='Optional: a prefix to use when creating novel IDs')
    args = parser.parse_args()
    return args

def make_transcript_dict(cursor):
    """ Format of dict:
            Key: tuple consisting of edges in transcript path
            Value: SQLite3 row from transcript table
         Note: if we later want to add novelty designation or other annotation
         info to the transcript table, this is a very good place to do that.
    """
    transcript_dict = {}
    query = """SELECT * FROM transcripts """
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

def match_all_transcript_vertices(chromosome, positions, location_dict, 
                                  run_info):
    """ Given a chromosome and a list of positions from the transcript in 5' to
        3' end order, this function looks for a matching vertex for each 
        position. """ 
 
    vertex_matches = []
 
    # Start and ends require permissive matching approach 
    start = positions.pop(0)
    end = positions.pop(-1)

    # Remaining mid-transcript positions go through strict matching process
    for position in positions:
        vertex_match = search_for_vertex_at_pos(chromosome, position, 
                                                location_dict)
        if vertex_match == None:
            # If no vertex matches the position, one is created.
            vertex_match = create_vertex(chromosome, position, run_info, 
                                         location_dict)

        # Add to running list of matches
        vertex_matches.append(vertex_match['location_ID'])

    return vertex_matches

def permissive_vertex_search(chromosome, position, strand, sj_pos, pos_type,
                             locations, cutoff, run_info):
    """ Given a position, this function tries to find a vertex match within the
        cutoff distance that also comes before the splice junction begins.
        If no vertex is found, the function returns None. """

    # Try a strict match first
    if (chromosome, position) in locations:
        match = locations[(chromosome, position)]['location_ID']
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
            closest = locations[candidate_match]['location_ID']
        # TODO: consider tiebreaker case.
        #elif abs(dist) == abs(min_dist)
        #    # Tiebreaker: 

    if closest == None:
        closest = create_vertex(chromosome, position, 
                                run_info, locations)['location_ID']
        min_dist = "NA"

    return closest, min_dist

            
def create_vertex(chromosome, position, run_info, location_dict):
    """ Creates a novel vertex and adds it to the location data structure. """
    run_info.vertex += 1
    new_ID = "%s-%d" % (run_info.prefix, run_info.vertex)
    new_vertex = {'location_ID': new_ID,
                  'genome_build': run_info.build,
                  'chromosome': chromosome,
                  'position': position}

    location_dict[(chromosome, position)] = new_vertex

    return new_vertex

def create_edge(vertex_1, vertex_2, edge_type, strand, edge_dict, run_info):
    """ Creates a novel edge and adds it to the edge data structure. """
    run_info.edge += 1
    new_ID = "%s-%d" % (run_info.prefix, run_info.edge)
    new_edge = {'edge_ID': new_ID,
                'v1': vertex_1,
                'v2': vertex_2,
                'edge_type': edge_type,
                'strand': strand }
    
    edge_dict[(vertex_1, vertex_2, edge_type)] = new_edge

    return new_edge              

def match_all_transcript_edges(vertices, strand, edge_dict, run_info):
    """ Given a list of vertex IDs from the transcript in 5' to
        3' end order, this function looks for a matching edge ID for each
        position. If none exists, it creates one. """

    edge_matches = []
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

        # Add to running list of matches
        edge_matches.append(edge_match['edge_ID'])

    return edge_matches

def search_for_transcript_suffix(edge_IDs, transcript_dict):
    """ Given a list of edges in a query transcript, determine whether it is
        a suffix for any transcript in the dict. We're looking for the gene ID
        here rather than worrying about exactly which transcript it came from.
    """
  
    if type(edge_IDs) is list:
        edge_IDs = tuple(edge_IDs)
    try:
        suffix_match = list(filter(lambda t: edge_IDs == t[-len(edge_IDs):],
                                     list(transcript_dict.keys())))[0]
        transcript = transcript_dict[suffix_match]
        gene_ID = transcript["gene_ID"]
        return gene_ID

    except:
        return None                 
    
def search_for_overlap_with_gene(chromosome, start, end, strand, 
                                 cursor, run_info):
    """ Given a start and an end value for an interval, query the database to
        determine whether the interval overlaps with a gene. If it overlaps
        more than one, prioritize same-strand first and foremost. """
    
    min_start = min(start, end)
    max_end = max(start, end)

    query = """ SELECT g.gene_ID,
                      v.vertex_ID,
                      loc.chromosome,
                      loc.position,
                      g.strand 
                   FROM genes as g
                   LEFT JOIN vertex as v ON g.gene_ID = v.gene_ID 
                   LEFT JOIN location as loc ON loc.location_ID = v.vertex_ID
                   WHERE loc.genome_build = '%s' 
                   AND loc.chromosome = '%s' 
                   AND loc.position >= %d 
                   AND loc.position <= %d """
    cursor.execute(query % (run_info.build, chromosome, min_start, max_end))
    matches = cursor.fetchall()

    if len(matches) == 0:
        return None
    
    # Among multiple matches, preferentially return the same-strand gene with 
    # the closest match to the 3' end 
    same_strand_matches = [ x for x in matches if x["strand"] == strand ]

    if strand == "+" and len(same_strand_matches) > 0 or \
        strand == "-" and len(same_strand_matches) == 0:
        matches = sorted(same_strand_matches, key = lambda x: x["position"], 
                         reverse = True )
        return(matches[0]["gene_ID"])

    else:   
        matches = sorted(same_strand_matches, key = lambda x: x["position"],
                         reverse = False )
        return(matches[0].gene_ID)


def search_for_transcript(edge_IDs, transcript_dict):
    """ Given the edge IDs that make up a query transcript (5' to 3' order), 
        look for a match in the transcript dict. 
        Return gene ID and transcript ID if found, and None if not. """

    if type(edge_IDs) is list:
        edge_IDs = tuple(edge_IDs)

    try:
        transcript = transcript_dict[edge_IDs]
        gene_ID = transcript["gene_ID"]
        transcript_ID = transcript["transcript_ID"]
        return gene_ID, transcript_ID

    except:
        return None, None

def check_inputs(options):
    """ Checks the input options provided by the user and makes sure that
        they are valid. Throw an error with descriptive help message if not."""

    # Make sure that the genome build exists in the provided TALON database.
    conn = sqlite3.connect(options.database)
    cursor = conn.cursor()
    cursor.execute(""" SELECT DISTINCT name FROM genome_build """)
    builds = cursor.fetchone()
    if options.build not in builds:
        build_names = ", ".join(list(annot_builds))
        raise ValueError("Please specify a genome build that exists in the" +
                          " database. The choices are: " + build_names)
    annot_builds = cursor.fetchall()

def init_run_info(cursor, genome_build, idprefix):
    """ Initializes a dictionary that keeps track of important run information
        such as the desired genome build, the prefix for novel identifiers,
        and the novel counters for the run. """

    run_info = dstruct.Struct()
    run_info.build = genome_build
    run_info.prefix = idprefix
    run_info.cutoff_5p = 500
    run_info.cutoff_3p = 300

    query = "SELECT * FROM counters WHERE category != 'genome_build'"
    cursor.execute(query)

    for counter in cursor.fetchall():
        counter_name = counter['category']
        run_info[counter_name] = counter['count']

    return run_info

 
def main():
    """ Runs program """

    options = get_args()
    check_inputs(options)

    # Fire up the database connection
    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Prepare data structures from the information stored in the database
    run_info = init_run_info(cursor, options.build, options.idprefix)
    location_dict = make_location_dict(options.build, cursor)
    edge_dict = make_edge_dict(cursor)
    transcript_dict = make_transcript_dict(cursor)

    chrom = "chr1"
    #vertex_IDs = [ 11, 12, 13, 14, 15, 16]
    strand = "+"
    #edge_IDs = match_all_transcript_edges(vertex_IDs, strand,
    #                                                    edge_dict, run_info)

    #gene_ID, transcript_ID = search_for_transcript(edge_IDs, transcript_dict)
    #edge_IDs = [11,12,13]
    #search_for_transcript_suffix(edge_IDs, transcript_dict)
    gene_ID = search_for_overlap_with_gene(chrom, 910, 1010, strand,
                                 cursor, run_info)
    print(gene_ID)
    conn.close()

if __name__ == '__main__':
    main()
