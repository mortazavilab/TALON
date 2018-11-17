# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# filter_talon_transcripts.py is a utility that filters the transcripts inside
# a TALON database to produce a transcript whitelist. This list can then be 
# used by downstream analysis tools to determine which transcripts and other
# features should be reported (for example in a GTF file).

from optparse import OptionParser
import sqlite3

def filter_talon_transcripts(database, annot, dataset_pairings = None,
                                              known_filtered = False, 
                                              novel_filtered = True,
                                              novel_multiexon_reqmt = True):

    # Create a set to keep track of whitelisted transcripts
    transcript_whitelist = set()

    # Connect to the database
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    # Query that joins transcript table to gene and transcript status
    query = """SELECT 
                   t.gene_ID,
                   t.transcript_ID,
                   t.path,
                   ga.value,
                   ta.value
               FROM transcripts t
               LEFT JOIN gene_annotations ga ON t.gene_ID = ga.ID
               LEFT JOIN transcript_annotations ta ON t.transcript_ID = ta.ID 
               WHERE (ga.annot_name = %s OR ga.annot_name = "talon_run")
                   AND ga.attribute = "gene_status" 
                   AND (ta.annot_name = %s  OR ta.annot_name = "talon_run") 
                   AND ta.attribute = "transcript_status";
            """
   
    annot_str = '"' + annot + '"'
    cursor.execute(query % (annot_str, annot_str))
    transcripts = cursor.fetchall()
    #print transcripts

    # Iterate over transcripts

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

def main():
    options = getOptions()
    filter_talon_transcripts(options.database, options.annot)


if __name__ == '__main__':
    main()
