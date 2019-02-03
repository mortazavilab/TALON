# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

import argparse
import sqlite3
import warnings

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

    args = parser.parse_args()
    return args

def make_location_dict(genome_build, cursor):
    """ Format of dict:
            Key: chromosome_pos
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

    return(location_dict)
    
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
        
    return vertex_matches

def search_for_edge(v1, v2, edge_type, edge_dict):
    """ """
    query_key = (v1, v2, edge_type)
    try:
        return edge_dict[query_key]
    except:
        return None

def search_for_all_transcript_vertices(position_pairs, cursor):
    """ Given a list of (chromosome, position) tuples, this function performs 
         a batch query of the location table in the provided database to 
         look for a matching vertex for each position. """ 
    
    query = """SELECT * FROM location
                   WHERE chromosome = ?
                   AND position = ? """

    for pair in position_pairs:
        cursor.execute(query, pair)
        vertex_matches = cursor.fetchall()


# TODO: validation of input options
def check_inputs(options):

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

def main():
    """ Runs program """

    options = get_args()
    check_inputs(options)
    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    location_dict = make_location_dict(options.build, cursor)
    edge_dict = make_edge_dict(cursor)
    conn.close()

    chrom = "chr1"
    pos1 = 500
    pos2 = 600
    v1 = search_for_vertex_at_pos(chrom, pos1, location_dict)["location_ID"]
    v2 = search_for_vertex_at_pos(chrom, pos2, location_dict)["location_ID"]

    edge_match = search_for_edge(v1, v2, "exon", edge_dict)

if __name__ == '__main__':
    main()
