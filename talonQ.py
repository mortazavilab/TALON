# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

import argparse
import sqlite3

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
        transcript_path = tuple(transcript["path"].split(","))
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
            run_info['vertex'] += 1
            new_ID = "%s-%d" % (run_info['prefix'], run_info['vertex'])
            new_vertex = {'location_ID': new_ID, 
                          'genome_build': run_info['build'],
                          'chromosome': chromosome,
                          'position': position}
            vertex_match = new_vertex
            location_dict[(chromosome, position)] = new_vertex

        # Add to running list of matches
        vertex_matches.append(vertex_match['location_ID'])

    return vertex_matches

# TODO: validation of input options
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

def init_run_info(cursor, build, prefix):
    """ Initializes a dictionary that keeps track of important run information
        such as the desired genome build, the prefix for novel identifiers,
        and the novel counters for the run. """

    run_info = {'build': build, 'prefix':prefix}
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
    conn.close()
    #chrom = "chr1"
    #pos1 = 500
    #pos2 = 600
    #v1 = search_for_vertex_at_pos(chrom, pos1, location_dict)["location_ID"]
    #v2 = search_for_vertex_at_pos(chrom, pos2, location_dict)["location_ID"]

    #edge_match = search_for_edge(v1, v2, "exon", edge_dict)

if __name__ == '__main__':
    main()
