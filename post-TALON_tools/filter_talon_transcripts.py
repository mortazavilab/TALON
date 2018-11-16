# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# filter_talon_transcripts.py is a utility that filters the transcripts inside
# a TALON database to produce a transcript whitelist. This list can then be 
# used by downstream analysis tools to determine which transcripts and other
# features should be reported (for example in a GTF file).

from optparse import OptionParser
import sqlite3

def filter_talon_transcripts(database, annot, dataset_pairings = None
                                              known_filtered = False, 
                                              novel_filtered = True,
                                              novel_multiexon_reqmt = True):

    # Create a set to keep track of whitelisted transcripts
    transcript_whitelist = set()

    # Connect to the database
    conn = sqlite3.connect(annot)
    cursor = conn.cursor()

    # Get known/novel status of transcripts

    # query that joins transcript table to transcript status
    #SELECT 
    #	transcripts.gene_ID,
    #	transcripts.transcript_ID,
    #	transcripts.path,
    #	transcript_annotations.ID,
    #	transcript_annotations.value
    #FROM transcripts
    #LEFT JOIN transcript_annotations ON transcripts.transcript_ID = transcript_annotations.ID 
    #WHERE transcript_annotations.annot_name = "KRT17_test" AND transcript_annotations.attribute = "transcript_status";


    # Disconnect from database
    conn.close()

    return

def getOptions():
    parser = OptionParser()
    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")
    parser.add_option("--annot", "-a", dest = "annot",
        help = "Which annotation version to use. Will determine which annotation transcripts are considered known or novel relative to. Note: must be in the TALON database.", type = "string")
    parser.add_option("--o", dest = "outfile", help = "Outfile name",
        metavar = "FILE", type = "string")


    (options, args) = parser.parse_args()
    return options

def main()
    options = getOptions()
    filter_talon_transcripts(options.database, options.annot)


if __name__ == '__main__':
    main()
