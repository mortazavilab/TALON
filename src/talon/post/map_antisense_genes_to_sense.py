# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# map_antisense_genes_to_sense.py is a utility that outputs the ID of the 
# corresponding sense gene for every antisense gene in the database

from optparse import OptionParser
import sqlite3
import warnings
import sys
import os
from pathlib import Path

def getOptions():
    parser = OptionParser()

    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")
    parser.add_option("--annot", "-a", dest = "annot",
        help = """Which annotation version to use. Will determine which
                  annotation is used to fetch gene names. Note: 
                  Must be in the TALON database.""", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output GTF",
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

def create_gene_name_dict(cursor, annot):
    """ Create a dictionary mapping TALON gene IDs to their names in the
        annot annotation"""

    cursor.execute("""SELECT ga.ID,
                             ga.value AS gene_name
                        FROM gene_annotations AS ga
                        WHERE ga.attribute = 'gene_name'
                        AND (ga.annot_name = '%s' OR ga.source = 'TALON')""" \
                        % (annot)) 

    gene_names = {}
    for entry in cursor.fetchall():
        ID = int(entry["ID"])
        name = entry["gene_name"]
        gene_names[ID] = name

    return gene_names

def main():
    options = getOptions()
    database = options.database
    annot = options.annot
    outfile = options.outprefix + "_antisense_mapping.csv"

    # Make sure that the input database exists!
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)

    if annot == None:
        raise ValueError("Please provide a valid annotation name")
    check_annot_validity(annot, database)
    
    # Connect to database and perform query
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    cursor.execute("""SELECT ga.ID As antisense_talon_ID,
                             ga.value AS sense_talon_ID
                        FROM gene_annotations AS ga
                        WHERE ga.attribute = 'gene_antisense_to_IDs'""")
    antisense_rows = cursor.fetchall()

    # Create a dict of gene names
    gene_name_dict = create_gene_name_dict(cursor, annot) 

    # Write antisense-sense pairs to file. When there is more than one sense match,
    # create separate lines
    o = open(outfile, 'w')
    o.write(",".join(["antisense_talon_ID", "sense_talon_ID", "sense_gene_name"]) + "\n")
    for entry in antisense_rows:
        sense_IDs = entry["sense_talon_ID"].split(",")
        for sense_ID in sense_IDs:
            o.write(",".join([str(entry["antisense_talon_ID"]),
                              str(sense_ID),
                              gene_name_dict[int(sense_ID)]]) + "\n")
          
    o.close() 
    
if __name__ == '__main__':
    main()
