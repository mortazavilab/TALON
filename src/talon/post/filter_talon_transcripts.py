# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# filter_talon_transcripts.py is a utility that filters the transcripts inside
# a TALON database to produce a transcript whitelist. This list can then be 
# used by downstream analysis tools to determine which transcripts and other
# features should be reported (for example in a GTF file).

from optparse import OptionParser
import sqlite3
import sys
import warnings
import os
from pathlib import Path
script_path = os.path.abspath(__file__)
main_path = "/".join(script_path.split("/")[0:-2])
sys.path.append(main_path)
import query_utils as qutils

def filter_talon_transcripts(database, annot, dataset_pairings = None,
                                              known_filtered = False,
                                              novel_filtered = True,
                                              novel_multiexon_reqmt = True):
    # Create a set to keep track of whitelisted transcripts
    # Each entry is a gene-transcript tuple
    transcript_whitelist = set()

    # Connect to the database
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    # If dataset pairings are not provided, simply make the pairing set
    # a list of every dataset in the database
    if dataset_pairings == None:
        cursor.execute("SELECT dataset_name FROM dataset")
        datasets = [str(x[0]) for x in cursor.fetchall()]
        pairing_list = [datasets]
    else:
        pairing_list = dataset_pairings

    # Filter transcripts separately for each dataset group
    for datasets in pairing_list:
        if len(datasets) <= 1 and novel_filtered == True:
            print("""Warning: Only one dataset in group. This means that
                   "only known transcripts will pass the filter 
                    for this group.""")
        else:
            print("Group: %s" % ", ".join([str(x) for x in datasets]))

        # First, accept all known transcripts and all NICs
        known_transcripts = qutils.fetch_known_transcripts_with_gene_label(cursor, datasets) 
        transcript_whitelist.update(known_transcripts)
        
        # Now, conditionally accept ISM, NNC, antisense, and intergenic transcripts 
        # (must be reproducible)
        NIC_transcripts = qutils.fetch_reproducible_NICs(cursor, datasets)
        reproducible_ISMs = qutils.fetch_reproducible_ISMs(cursor, datasets)
        reproducible_NNCs = qutils.fetch_reproducible_NNCs(cursor, datasets)
        reproducible_antisense = qutils.fetch_reproducible_antisense(cursor, datasets)
        reproducible_intergenic = qutils.fetch_reproducible_intergenic(cursor, datasets)
        transcript_whitelist.update(NIC_transcripts)
        transcript_whitelist.update(reproducible_ISMs)
        transcript_whitelist.update(reproducible_NNCs)
        transcript_whitelist.update(reproducible_antisense)
        transcript_whitelist.update(reproducible_intergenic)

    return transcript_whitelist


def process_pairings(pairings_file):
    """ Reads in pairings from the comma-delimited pairings file and creates 
        a list of lists """

    pairings = []
    with open(pairings_file, 'r') as f:
        for group in f:
            group = group.strip().split(',')
            pairings.append(group)
    return pairings

def getOptions():
    parser = OptionParser()
    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")
    parser.add_option("--annot", "-a", dest = "annot",
        help = """Which annotation version to use. Will determine which 
                  annotation transcripts are considered known or novel 
                  relative to. Note: must be in the TALON database.""", 
        type = "string")
    parser.add_option("--pairings", "-p",  dest = "pairings_file",
        help = """Optional: A file indicating which datasets should be 
                  considered together when filtering novel transcripts 
                  (i.e. biological replicates). 
                  Format: Each line of the file constitutes a group, with 
                  member datasets separated by commas. 
                  If no file is provided, then novel transcripts appearing in 
                  any two datasets will be accepted.""", 
        metavar = "FILE", type = "string", default = None)

    parser.add_option("--o", dest = "outfile", help = "Outfile name",
        metavar = "FILE", type = "string")


    (options, args) = parser.parse_args()
    return options

def check_annot_validity(annot, database):
    """ Make sure that the user has entered a correct annotation name """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    cursor.execute("SELECT DISTINCT annot_name FROM gene_annotations")
    annotations = [str(x[0]) for x in cursor.fetchall()]
    conn.close()

    if "TALON" in annotations:
        annotations.remove("TALON")

    if annot == None:
        message = "Please provide a valid annotation name. " + \
                  "In this database, your options are: " + \
                  ", ".join(annotations)
        raise ValueError(message)

    if annot not in annotations:
        message = "Annotation name '" + annot + \
                  "' not found in this database. Try one of the following: " + \
                  ", ".join(annotations)
        raise ValueError(message)

    return

def main():
    options = getOptions()
    database = options.database
    annot = options.annot

    # Make sure that the input database exists!
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)    

    # Make sure that the provided annotation name is valid
    check_annot_validity(annot, database)

    if options.pairings_file != None:
        pairings = process_pairings(options.pairings_file)
        whitelist = filter_talon_transcripts(database, annot,
                                             dataset_pairings = pairings)
    else:
        whitelist = filter_talon_transcripts(options.database, options.annot)

    # Write transcript IDs to file
    o = open(options.outfile, 'w')
    print("Writing whitelisted gene-transcript-category pairs to " + options.outfile + "...")
    for transcript in whitelist:
        o.write(",".join([str(x) for x in transcript]) + "\n")
    o.close()
    

if __name__ == '__main__':
    main()
