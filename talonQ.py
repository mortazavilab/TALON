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
        key = "%s_%d" % (chromosome, position)
        location_dict[key] = location

    return(location_dict)
    

def search_for_vertex_at_pos(chromosome, position, location_dict):
    """ Given a chromosome and a position (1-based), this function queries the 
        location dict to determine whether a vertex 
        fitting those criteria exists. Returns the row if yes, and __ if no.
    """
    query_key = "%s_%d" % (chromosome, position)
    try:
        return location_dict[query_key]
    except:
        return None
        
    #query = """SELECT * FROM location 
    #               WHERE genome_build = ?
    #               AND chromosome = ? 
    #               AND position = ? """
    
    #cursor.execute(query, [genome_build, chromosome, pos])
    #vertex_matches = cursor.fetchall()
    return vertex_matches

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
        print(vertex_matches)


# TODO: validation of input options

def main():
    """ Runs program """

    options = get_args()
    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    make_location_dict("hg38", cursor)

    conn.close()

if __name__ == '__main__':
    main()
