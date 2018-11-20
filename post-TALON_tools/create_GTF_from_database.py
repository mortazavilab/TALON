# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# create_GTF_from_database.py is a utility that outputs the genes, transcripts,
# and exons stored a TALON database into a GTF annotation file. 

import itertools
import operator
from optparse import OptionParser
import sqlite3
import warnings
import filter_talon_transcripts as filt
import sys
sys.path.append("..")
import talon as TALON
import os

def getOptions():
    parser = OptionParser()

    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")

    parser.add_option("--build", "-b", dest = "build",
        help = "Genome build to use. Note: must be in the TALON database.",
        type = "string")

    parser.add_option("--annot", "-a", dest = "annot",
        help = """Which annotation version to use. Will determine which
                  annotation transcripts are considered known or novel
                  relative to. Note: must be in the TALON database.""",
        type = "string")

    parser.add_option("--filter", dest ="filtering", action='store_true',
                      help = "If this option is set, the transcripts in the  \
                      database will be filtered prior to GTF creation \
                      (for more information, see filter_talon_transcripts.py)")

    parser.add_option("--observed", dest ="observed", action='store_true',
                      help = "If this option is set, the GTF file will only  \
                      include transcripts that were observed in at least one \
                      dataset. Redundant if 'filter' option is set.")   
 
    parser.add_option("--pairings", "-p",  dest = "pairings_file",
        help = """Optional (only relevant if filter = true): A file indicating 
                  which datasets should be considered together when filtering 
                  novel transcripts (i.e. biological replicates).
                  Format: Each line of the file constitutes a group, with
                  member datasets separated by commas.
                  If no file is provided, then novel transcripts appearing in
                  any two datasets will be accepted.""",
        metavar = "FILE", type = "string", default = None)

    parser.add_option("--o", dest = "outprefix", help = "Prefix for output GTF",
        metavar = "FILE", type = "string")


    (options, args) = parser.parse_args()
    return options

def handle_filtering(options):
    """ Determines which transcripts to allow in the GTF file. This can be done
        in two different ways. If filtering is turned off, then all of the 
        transcripts in the database go on the whitelist (modified by 'observed'
        option). If filtering is turned on, known transcripts and novel 
        transcripts appearing in any two datasets will be accepted. This can be 
        tuned further by providing a dataset pairing file, but this is optional.
    """

    database = options.database
    annot = options.annot
    pairings = options.pairings_file

    # If filtering is turned on, run filtering step to get whitelisted transcripts
    if options.filtering == True:
        print "--------------------------------------"
        print "Running transcript filtering process:"
        if pairings != None:
            whitelist = filt.filter_talon_transcripts(database, annot, dataset_pairings = pairings)
        else:
            whitelist = filt.filter_talon_transcripts(database, annot)
        print "--------------------------------------"

    # If "observed" option is on, all observed transcripts are included
    elif options.observed == True:

        conn = sqlite3.connect(database)
        cursor = conn.cursor()
        query = """ 
                SELECT
                    t.gene_ID,
                    t.transcript_ID,
                    COUNT(*)
                FROM transcripts t
                INNER JOIN abundance ON t.transcript_ID = abundance.transcript_ID
                GROUP BY t.transcript_ID;
                """ 
        cursor.execute(query)
        whitelist = [(str(x[0]),str(x[1])) for x in cursor.fetchall()]

        conn.close()

    # Otherwise, the whitelist will be every transcript ID in the database
    else:

        conn = sqlite3.connect(database)
        cursor = conn.cursor()

        cursor.execute("SELECT gene_ID,transcript_ID FROM transcripts")
        whitelist = [(str(x[0]),str(x[1])) for x in cursor.fetchall()]

        conn.close()

    # Sort the whitelist tuples on gene ID
    sorted_whitelist = sorted(whitelist, key=lambda x: x[0])

    return sorted_whitelist

def create_outname(options):
    """ Creates filename for the output GTF that reflects the input options that
        were used. """
    
    outname = options.outprefix + "_talon"
    if options.filtering == True:
        outname = "_".join([ outname, "filtered" ])
    elif options.observed == True: 
        outname = "_".join([ outname, "observedOnly" ])
    
    outname += ".gtf"
    return outname

def create_gtf(database, annot, genome_build, whitelist):

    # Divide tuples into separate sublists based on gene
    gene_groups = []
    for key,group in itertools.groupby(whitelist,operator.itemgetter(0)):
        gene_groups.append(list(group))

    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    
    exon_tree = read_edges(cursor, genome_build, "exon")
    intron_tree = read_edges(cursor, genome_build, "intron")
    transcripts = read_transcripts(cursor, exon_tree, intron_tree)
 



def main():
    options = getOptions()
    database = options.database
    annot = options.annot
    build = options.build
    outfile = create_outname(options)

    # Determine which transcripts to include
    transcript_whitelist = handle_filtering(options)

    create_gtf(database, annot, build, transcript_whitelist)    
    

if __name__ == '__main__':
    main()
