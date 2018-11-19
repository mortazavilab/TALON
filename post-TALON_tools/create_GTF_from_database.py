# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# create_GTF_from_database.py is a utility that outputs the genes, transcripts,
# and exons stored a TALON database into a GTF annotation file. 

from optparse import OptionParser
import sqlite3
import warnings

def getOptions():
    parser = OptionParser()
    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")

    parser.add_option("--annot", "-a", dest = "annot",
        help = """Which annotation version to use. Will determine which
                  annotation transcripts are considered known or novel
                  relative to. Note: must be in the TALON database.""",
        type = "string")

    parser.add_option("--filter", dest ="filter", action='store_true',
                      help = "If this option is set, the transcripts in the  \
                      database will be filtered prior to GTF creation \
                      (for more information, see filter_talon_transcripts.py)")
    
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

def main():
    options = getOptions()

if __name__ == '__main__':
    main()
