# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

from optparse import OptionParser
import intervaltree

def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "infiles", 
                      help = "Comma-delimited list of input SAM files",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--o", dest = "outfile",
                      help = "output file", metavar = "FILE", type = "string", default = "out.txt")
    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    
    with open(options.infile, 'r') as f:
        pass

if __name__ == '__main__':
    main()
