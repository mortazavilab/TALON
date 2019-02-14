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
            closest = locations[candidate_match]#['location_ID']

    if closest == None:
        min_dist = None

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

def create_transcript(gene_ID, edge_IDs, vertex_IDs, transcript_dict, run_info):
    """Creates a novel transcript and adds it to the transcript data structure.
    """
    run_info.transcripts += 1
    new_ID = "%s-%d" % (run_info.prefix, run_info.transcripts)
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
        return novelty[0] == 1

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
        We're looking for the gene ID here rather than worrying about exactly 
        which transcript it came from.
    """  

    try:
        suffix_match = list(filter(lambda t: edge_IDs[1:-1] == t[-len(edge_IDs) + 1:-1],
                                     list(transcript_dict.keys())))[0]
        transcript = transcript_dict[suffix_match]
        gene_ID = transcript["gene_ID"]
        return gene_ID

    except:
        return None


def search_for_transcript_prefix(edge_IDs, transcript_dict):
    """ Given a list of edges in a query transcript, determine whether it is
        a prefix for any transcript in the dict. It is OK for the first and 
        last exon IDs to be different, but the splice junctions must match.
        We're looking for the gene ID here rather than worrying about exactly
        which transcript it came from (since in some cases it could be more 
        than one).
    """

    try:
        prefix_match = list(filter(lambda t: edge_IDs[1:-1] == t[1:len(edge_IDs)-1],
                                     list(transcript_dict.keys())))[0]
        transcript = transcript_dict[prefix_match]
        gene_ID = transcript["gene_ID"]
        return gene_ID

    except:
        return None

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

    edge_IDs = edge_IDs[1:-1]

    ISM_matches = list(filter(lambda t: bytes(edge_IDs) in bytes(t),
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
        overlap. Antisense matches may be returned if there is no sam strand
        option. """

    min_start = min(start, end)
    max_end = max(start, end)
    query_interval = [min_start, max_end]

    query = """ SELECT gene_ID, 
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
                       AND loc.chromosome = '%s'
                       GROUP BY g.gene_ID)
                 WHERE (start <= %d AND end >= %d) OR
                       (start >= %d AND end <= %d) OR
                       (start >= %d AND start <= %d) OR
                       (end >= %d AND end <= %d);"""

    cursor.execute(query % (run_info.build, chromosome, min_start, max_end,
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

def process_FSM(edge_IDs, vertex_IDs, vertex_novelties, transcript_dict,
                       run_info ):
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
        return gene_ID, transcript_ID, None
    
    # Next, look for FSM with different ends. This triggers a novel transcript
    # that has novel 5' or 3' ends
    gene_ID, transcript_matches = search_without_transcript_ends(edge_IDs,
                                                               transcript_dict)
    if transcript_matches != None:
        # Create a new transcript
        novel_transcript = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                             transcript_dict, run_info)

        transcript_IDs = ",".join([str(match[0]) for match in transcript_matches])
        novelty.append((novel_transcript['transcript_ID'],"FSM", transcript_IDs))
        if vertex_novelties[0] == 1:
            novelty.append(((novel_transcript["transcript_ID"],"5p", None)))
        if vertex_novelties[-1] == 1:
            novelty.append(((novel_transcript["transcript_ID"],"3p", None)))
                
        return gene_ID, novel_transcript["transcript_ID"], novelty

    return None,None,None

def process_ISM(edge_IDs, vertex_IDs, vertex_novelties, transcript_dict, 
                       run_info ):
    """ Given a transcript that is known to be a full splice match (FSM) or an
        incomplete splice match (ISM) on the basis of its junctions,
        find the best gene and transcript match for it. Also delineate its 
        type(s) of novelty, which will later be added to the database. 
        In the case of an ISM, a novel transcript will be created."""

    novelty = []

    
    #gene_ID, transcript_matches = search_for_ISM(edge_IDs, transcript_dict)
    #    if transcript_matches != None:

            # Create a new transcript
    #        novel_transcript = create_transcript(gene_ID, edge_IDs, vertex_IDs,
                                                 #transcript_dict, run_info)

    #    else:
            # If transcript is not an ISM, throw an error
    #        raise ValueError("Transcript is supposedly an FSM or ISM, " + \
                                # "but no match was found")
    #        novelty = [(novel_transcript["transcript_ID"], "ISM", 
    #                    transcript_IDs)]

        # Look for suffix ISM
    #    gene_ID, transcript_matches = search_for_transcript_suffix(edge_IDs, 
                                         #                      transcript_dict)
    #    if transcript_matches != None:
    #        transcript_IDs = ",".join([match[0] for match in transcript_matches])
    #        novelty.append(((novel_transcript["transcript_ID"], "ISM_suffix", 
    #                         transcript_IDs)))

        # Look for prefix ISM
    #    gene_ID, transcript_matches = search_for_transcript_prefix(edge_IDs,
                                          #                     transcript_dict)
    #    if transcript_matches != None:
    #        transcript_IDs = ",".join([match[0] for match in transcript_matches])
    #        novelty.append(((novel_transcript["transcript_ID"],
    #                         "ISM_prefix", transcript_IDs))
            

    # Check on ends:
    #if vertex_novelties[0] == 1:
    #    novelty.append(((novel_transcript["transcript_ID"],"5p", None))
    #if vertex_novelties[-1] == 1:
    #    novelty.append(((novel_transcript["transcript_ID"],"3p", None))

    #return gene_ID, novel_transcript["transcript_ID"], novelty
    

def identify_transcript(read_ID, chrom, positions, strand, location_dict, edge_dict,
                        transcript_dict, run_info):
    """ Inputs:
        - Information about the query transcript
          - chromosome
          - list of positions
          - strand
        - Data structures
          - location_dict (position --> vertex)
          - edge_dict (v1_v2_edgetype --> edge)
          - possibly a gene dict
          - transcript_dict
          - run_info

       Outputs:
          - Assigned gene ID
          - Assigned transcript ID
          - 'observed' entry (to be added to database) 
          - novelty entries (to be added to database)
    """
    transcript_novelty = set()
    n_exons = (len(positions) + 1)/2.0
 
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
    ends_novel = (e_novelty[0] + e_novelty[-1]) > 0

    if len(edge_IDs) == 1:
        # Handle monoexonic case. Check for full monoexonic match
        pass

    elif all_SJs_known:
        print("Transcript is either an FSM or an ISM")
        gene_ID, transcript_ID, novelty = process_FSM(edge_IDs, vertex_IDs, 
                                                             v_novelty, 
                                                             transcript_dict,
                                                             run_info) 
        
        return gene_ID, transcript_ID, novelty
    elif all_exons_known and not(all_SJs_known):
        print("Transcript is definitely Novel in Catalog (NIC)")
    elif not(all_exons_known) and not(all_SJs_known) and n_exons > 1:
        print("Transcript is definitely Novel Not in Catalog (NNC)")
    else:
        print("Transcript is genomic and/or antisense")

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
    strand = "+"
    read_ID = "toy_read"
    positions = ( 1, 100, 500, 600, 900, 1500 )
    #gene_ID, transcript_ID, novelty = identify_transcript(read_ID,
    #                                           chrom, positions, strand,
    #                                           location_dict, edge_dict,
    #                                           transcript_dict, run_info)

    edge_IDs = (1, 2, 3, 4, 500) # Last edge is novel
    vertex_IDs = (1, 2, 3, 4, 5, 500) # Last vertex is novel
    v_novelty = (0, 0, 0, 0, 0, 1)

    print(transcript_dict)
    gene_ID, transcript_ID, novelty = process_FSM(edge_IDs, vertex_IDs,
                                                      v_novelty, transcript_dict,
                                                      run_info)
    print(gene_ID)
    print(novelty)
    conn.close()

if __name__ == '__main__':
    main()
