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

def search_for_vertex_at_pos(genome_build, chromosome, pos, cursor):
    """ Given a chromosome and a position (1-based), this function queries the 
        location table in the provided database to determine whether a vertex 
        fitting those criteria exists. Returns the row(s) if yes, and __ if no.
    """
    query = """SELECT * FROM location 
                   WHERE genome_build = ?
                   AND chromosome = ? 
                   AND position = ? """
    
    cursor.execute(query, [genome_build, chromosome, pos])
    vertex_matches = cursor.fetchall()
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

    
    pairs = [("chr1",0), ("chr1", 1), ("chr1", 1000)]
    search_for_all_transcript_vertices(pairs, cursor)
    #search_for_vertex_at_pos("chr1", 1, cursor)
    #search_for_vertex_at_pos("chr1", 0, cursor)

    conn.close()

if __name__ == '__main__':
    main()
