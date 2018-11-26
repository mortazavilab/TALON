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
import pdb
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

def get_annotations(database, feat_type, annot, whitelist = None):
    """ Extracts annotations from the gene/transcript/exon annotation table of
        the database (depending on choice of feat_type). Limited to rows where
        the annot_name column matches the value of annot.

        Returns:
            annotation_dict: dictionary data structure in which the keys are
                             gene/transcript/exon TALON IDs (depending on 
                             choice of feat_type) and the value is a list of 
                             annotation tuples.
    """
    # Fetch the annotations
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    table_name = feat_type + "_annotations"

    if whitelist == None:
        query = "SELECT * FROM " + table_name + " WHERE annot_name = '" + annot + "'"
    else:
        whitelist_string = "(" + ','.join(whitelist) + ")"
        query = "SELECT * FROM " + table_name + " WHERE annot_name = '" + annot + \
                "' AND ID IN " + whitelist_string

    cursor.execute(query)
    annotation_tuples = cursor.fetchall()

    # Sort based on ID
    sorted_annotations = sorted(annotation_tuples, key=lambda x: x[0]) 

    # Group by ID and store in a dictionary
    ID_groups = {}
    for key,group in itertools.groupby(sorted_annotations,operator.itemgetter(0)):
        ID_groups[str(key)] = list(group)

    return ID_groups

def get_gene_2_transcripts(database, genome_build, whitelist):
    """ Creates a dictionary mapping gene IDs to the transcripts that belong to
        them. The columns in each tuple are:
            0: gene ID
            1: transcript ID
            2: chromosome
            3: start position (min of 5' and 3')
            4: end position (max of 5' and 3')
            5: strand
            6: edge path
 """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    whitelist_string = "(" + ','.join(whitelist) + ")"
    query = """
            SELECT 
               t.gene_ID,
               t.transcript_ID,
               loc1.chromosome,
               MIN(loc1.position,loc2.position),
               MAX(loc1.position,loc2.position),
               loc1.strand,
               t.path
           FROM transcripts t
           LEFT JOIN location loc1 ON t.start_vertex = loc1.location_ID
           LEFT JOIN location loc2 ON t.end_vertex = loc2.location_ID
           WHERE loc1.genome_build = '""" + genome_build + """' AND 
           loc2.genome_build = '""" + genome_build + \
           """' AND t.transcript_ID IN """ + whitelist_string 
           
    cursor.execute(query)
    transcript_tuples = cursor.fetchall()

    # Sort based on gene ID
    sorted_transcript_tuples = sorted(transcript_tuples, key=lambda x: x[0])

    # Group by gene ID and store in a dictionary
    gene_groups = {}
    for key,group in itertools.groupby(transcript_tuples,operator.itemgetter(0)):
        # Sort by transcript start position
        gene_groups[str(key)] = sorted(list(group), key=lambda x: x[3])

    conn.close()
    return gene_groups 

def fetch_exon_locations(database, genome_build):
    """ Queries the database to create a dictionary mapping exon IDs to 
        the chromosome, start, end, and strand of the exon """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    query = """
            SELECT 
                e.edge_ID,
                loc1.chromosome,
                MIN(loc1.position,loc2.position),
                MAX(loc1.position,loc2.position),
                loc1.strand
             FROM edge e
             LEFT JOIN location loc1 ON e.v1 = loc1.location_ID
             LEFT JOIN location loc2 ON e.v2 = loc2.location_ID
             WHERE loc1.genome_build = '""" + genome_build + """' AND
             loc2.genome_build = '""" + genome_build + \
             """' AND e.edge_type = 'exon';"""
    
    cursor.execute(query)
    exon_location_tuples = cursor.fetchall()

    # Create dictionary
    exon_locations = {}
    for loc_tuple in exon_location_tuples:
        exon_ID = str(loc_tuple[0])
        exon_locations[exon_ID] = loc_tuple[1:]

    conn.close() 
    return exon_locations

def create_gtf(database, annot, genome_build, whitelist):

    # Divide transcript tuples into separate sublists based on gene
    gene_groups = []
    for key,group in itertools.groupby(whitelist,operator.itemgetter(0)):
        gene_groups.append(list(group))

    # Get gene, transcript, and exon annotations
    gene_whitelist = [ x[0] for x in whitelist ]
    transcript_whitelist = [ x[1] for x in whitelist ]

    gene_annotations = get_annotations(database, "gene", annot, 
                                       whitelist = gene_whitelist)  
    transcript_annotations = get_annotations(database, "transcript", annot,
                                             whitelist = transcript_whitelist) 
    exon_annotations = get_annotations(database, "exon", annot)

    # Get transcript data from the database
    gene_2_transcripts = get_gene_2_transcripts(database, genome_build, 
                         transcript_whitelist)

    # Get exon location info from database
    exon_ID_2_location = fetch_exon_locations(database, genome_build)
 
    # -------------------------------------------------------------

    # Create a GTF entry for every gene
    for gene_ID, transcript_tuples in gene_2_transcripts.iteritems():

        gene_GTF = get_gene_GTF_entry(gene_ID, transcript_tuples, gene_annotations)
        print gene_GTF
    
        # Create a GTF entry for every transcript
        for transcript_entry in transcript_tuples:
            transcript_GTF_line = get_transcript_GTF_entry(transcript_entry, 
                                                           transcript_annotations)
            print transcript_GTF_line
            transcript_edges = str(transcript_entry[6]).split(",")
            

            # Create a GTF entry for every exon (skip introns)
            for exon_ID in transcript_edges[::2]:
                exon_GTF_line = get_exon_GTF_entry(exon_ID, exon_ID_2_location, 
                                                   exon_annotations)
                print exon_GTF_line   
                exit()
        exit()

def get_gene_GTF_entry(gene_ID, associated_transcript_tuples, gene_annotations):
    """ Creates a GTF annotation entry for the given gene """

    # Get attribute info for gene
    curr_annot = gene_annotations[gene_ID]

    # GTF fields
    chromosome = associated_transcript_tuples[0][2]
    source = [ str(x[-1]) for x in curr_annot if str(x[-2]) == "source"][0]
    feature = "gene"
    start = str(associated_transcript_tuples[0][3])
    end = str(associated_transcript_tuples[-1][4])
    score = "."
    strand = associated_transcript_tuples[0][5]
    frame = "."

    GTF = "\t".join([chromosome, source, feature, start, end, score, strand, frame])
    return GTF


def get_transcript_GTF_entry(transcript_entry, transcript_annotations):
    """ Creates a GTF annotation entry for the given transcript """

    transcript_ID = str(transcript_entry[1])
    curr_transcript_annot = transcript_annotations[transcript_ID]

    # GTF fields for transcript
    chromosome = str(transcript_entry[2])
    source = [ str(x[-1]) for x in curr_transcript_annot if str(x[-2]) == "source"][0]
    feature = "transcript"
    start = str(transcript_entry[3])
    end = str(transcript_entry[4])
    score = "."
    strand = transcript_entry[5]
    frame = "."

    GTF = "\t".join([chromosome, source, feature, start, end, score, strand, frame])
    return GTF

def get_exon_GTF_entry(exon_ID, exon_ID_2_location, exon_annotations):
    """ Creates a GTF annotation entry for the given exon """

    curr_exon_location = exon_ID_2_location[exon_ID]
    curr_exon_annot = exon_annotations[exon_ID]
    chromosome = str(curr_exon_location[0])
    source = [ str(x[-1]) for x in curr_exon_annot if str(x[-2]) == "source"][0]
    feature = "exon"
    start = str(curr_exon_location[1])
    end = str(curr_exon_location[2])
    score = "."
    strand = curr_exon_location[3]
    frame = "."

    GTF = "\t".join([chromosome, source, feature, start, end, score, strand, frame])
    return GTF

def main():
    options = getOptions()
    database = options.database
    annot = options.annot
    build = options.build
    outfile = create_outname(options)

    if build == None:
        raise ValueError("Please provide a valid genome build name")
    if annot == None:
        raise ValueError("Please provide a valid annotation name")


    # Determine which transcripts to include
    transcript_whitelist = handle_filtering(options)

    create_gtf(database, annot, build, transcript_whitelist)    
    

if __name__ == '__main__':
    main()
