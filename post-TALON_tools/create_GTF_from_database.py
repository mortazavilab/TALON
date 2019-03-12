# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# create_GTF_from_database.py is a utility that outputs the genes, transcripts,
# and exons stored a TALON database into a GTF annotation file. 

import copy
import itertools
import operator
from optparse import OptionParser
import sqlite3
import warnings
import filter_talon_transcripts as filt
import pdb
import sys
import os
from pathlib import Path
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, os.pardir)))

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
        print("--------------------------------------")
        print("Running transcript filtering process:")
        if pairings != None:
            pairings_list = []
            with open(pairings, 'r') as f:
                for pairing in f:
                    pairings_list.append(pairing.strip().split(","))

            whitelist = filt.filter_talon_transcripts(database, annot, 
                                  dataset_pairings = pairings_list)
        else:
            whitelist = filt.filter_talon_transcripts(database, annot)
        print("--------------------------------------")

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
        whitelist = [(x[0],x[1]) for x in cursor.fetchall()]

        conn.close()

    # Otherwise, the whitelist will be every transcript ID in the database
    else:

        conn = sqlite3.connect(database)
        cursor = conn.cursor()

        cursor.execute("SELECT gene_ID,transcript_ID FROM transcripts")
        whitelist = [(x[0],x[1]) for x in cursor.fetchall()]

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
        query = "SELECT * FROM " + table_name + " WHERE annot_name = '" + annot + \
         "' OR source = 'TALON'"
    else:
        whitelist_string = "(" + ','.join([str(x) for x in whitelist]) + ")"
        query = "SELECT * FROM " + table_name + " WHERE (annot_name = '" + annot + \
                "' OR source = 'TALON') AND ID IN " + whitelist_string

    cursor.execute(query)
    annotation_tuples = cursor.fetchall()

    # Sort based on ID
    sorted_annotations = sorted(annotation_tuples, key=lambda x: x[0]) 

    # Group by ID and store in a dictionary
    ID_groups = {}
    for key,group in itertools.groupby(sorted_annotations,operator.itemgetter(0)):
        ID_groups[key] = list(group)

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
            7. n_exons
 """

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    whitelist_string = "(" + ','.join([str(x) for x in whitelist]) + ")"
    query = """
            SELECT 
               t.gene_ID,
               t.transcript_ID,
               loc1.chromosome,
               MIN(loc1.position,loc2.position) AS min_pos,
               MAX(loc1.position,loc2.position) AS max_pos,
               genes.strand,
               t.jn_path,
               t.start_exon,
               t.end_exon,
               t.n_exons
           FROM transcripts t
           LEFT JOIN location loc1 ON t.start_vertex = loc1.location_ID
           LEFT JOIN location loc2 ON t.end_vertex = loc2.location_ID
           LEFT JOIN genes ON t.gene_ID = genes.gene_ID
           WHERE loc1.genome_build = '""" + genome_build + """' AND 
           loc2.genome_build = '""" + genome_build + \
           """' AND t.transcript_ID IN """ + whitelist_string 
           
    cursor.execute(query)
    transcript_tuples = cursor.fetchall()

    # Sort based on gene ID
    sorted_transcript_tuples = sorted(transcript_tuples, key=lambda x: x["gene_ID"])    

    gene_groups = {}
    for key,group in itertools.groupby(sorted_transcript_tuples,operator.itemgetter(0)):
        # Sort by transcript start position
        gene_groups[key] = sorted(list(group), key=lambda x: x["min_pos"])
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
                e.strand
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
        exon_ID = loc_tuple[0]
        exon_locations[exon_ID] = loc_tuple[1:]

    conn.close() 
    return exon_locations

def create_gtf(database, annot, genome_build, whitelist, outfile):

    # Create separate gene and transcript whitelists
    gene_whitelist = []
    transcript_whitelist = []
    for key,group in itertools.groupby(whitelist,operator.itemgetter(0)):
        gene_whitelist.append(key)
        for id_tuple in list(group):
            transcript_whitelist.append(id_tuple[1])

    # Get gene, transcript, and exon annotations
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

    o = open(outfile, 'w')

    # Create a GTF entry for every gene
    for gene_ID, transcript_tuples in gene_2_transcripts.items():
        curr_annot = gene_annotations[gene_ID]

        gene_annotation_dict = {}
        for annot in curr_annot:
            attribute = annot[3]
            value = annot[4]
            gene_annotation_dict[attribute] = value
        gene_GTF_line = get_gene_GTF_entry(gene_ID, transcript_tuples, 
                                      copy.copy(gene_annotation_dict))
        o.write(gene_GTF_line + "\n")
    
        # Create a GTF entry for every transcript of this gene
        for transcript_entry in transcript_tuples:
            transcript_ID = transcript_entry["transcript_ID"]
            curr_transcript_annot = transcript_annotations[transcript_ID]
            
            transcript_annotation_dict = {}
            for annot in curr_transcript_annot:
                attribute = annot[3]
                value = annot[4]
                transcript_annotation_dict[attribute] = value
            transcript_GTF_line = get_transcript_GTF_entry(transcript_entry, 
                                            copy.copy(gene_annotation_dict),
                                       copy.copy(transcript_annotation_dict))
            o.write(transcript_GTF_line + "\n")
            if transcript_entry["n_exons"] != None:
                transcript_edges = [str(transcript_entry["start_exon"])] + \
                                   str(transcript_entry["jn_path"]).split(",")+ \
                                   [str(transcript_entry["end_exon"])]
            else:
                transcript_edges = [transcript_entry["start_exon"]]

            # Create a GTF entry for every exon of this transcript (skip introns)
            exon_num = 1
            for exon_ID in transcript_edges[::2]:
                exon_ID = int(exon_ID)
                curr_exon_annot = exon_annotations[exon_ID]

                exon_annotation_dict = {}
                for annot in curr_exon_annot:
                    attribute = annot[3]
                    value = annot[4]
                    exon_annotation_dict[attribute] = value

                
                exon_GTF_line = get_exon_GTF_entry(gene_ID, transcript_ID, 
                                                   exon_ID, exon_num,
                                                   exon_ID_2_location,
                                                   copy.copy(gene_annotation_dict), 
                                                   copy.copy(transcript_annotation_dict),
                                                   exon_annotation_dict)
                o.write(exon_GTF_line + "\n")
                exon_num += 1
    o.close()
    return

def make_descriptor_string(attribute, value):
    """ Create a key-value string to form part of a GTF entry.
        Example:    gene_id and ENSG00000117676.13
                          becomes
                    gene_id "ENSG00000117676.13";
    """

    return str(attribute) + ' "' + str(value) + '";'

def format_GTF_tag_values_for_gene(gene_ID, annotation_dict):
    """ Parses the annotations for this gene, and supplements them where 
        necessary for novel transcripts """
 
    attributes = []

    # Mandatory: Gene ID
    if "gene_id" in annotation_dict:
        gene_ID_val = annotation_dict.pop("gene_id")
    else:
        gene_ID_val = gene_ID
    attributes.append(make_descriptor_string("gene_id", gene_ID_val))

    # Mandatory: Gene Name
    if "gene_name" in annotation_dict:
        gene_name = annotation_dict.pop("gene_name")
    else:
        gene_name = "TALON-" + str(gene_ID)
    attributes.append(make_descriptor_string("gene_name", gene_name))

    # Mandatory: Gene Status
    gene_status = annotation_dict.pop("gene_status")
    attributes.append(make_descriptor_string("gene_status", gene_status))
   
    # Gene type
    if "gene_type" in annotation_dict:
        gene_type = annotation_dict.pop("gene_type")
        attributes.append(make_descriptor_string("gene_type", gene_type))

    # Source
    if "source" in annotation_dict:
        source = annotation_dict.pop("source")
    else:
        source = "TALON"
        attributes.append(make_descriptor_string("source", source))
    
    # TALON Gene ID
    attributes.append(make_descriptor_string("talon_gene", gene_ID)) 

    # Add any remaining annotations
    for attribute,value in sorted(annotation_dict.items()):
        attributes.append(make_descriptor_string(attribute, value))

    return attributes

def format_GTF_tag_values_for_transcript(gene_ID, transcript_ID, gene_annot_dict,
                                         transcript_annot_dict):
    """ Parses the annotations for this transcript, and supplements them where
        necessary for novel transcripts """

    attributes = []

    # Mandatory: Gene ID
    if "gene_id" in gene_annot_dict:
        gene_ID_val = gene_annot_dict.pop("gene_id")
    else:
        gene_ID_val = gene_ID
    attributes.append(make_descriptor_string("gene_id", gene_ID_val))

    # Mandatory: Transcript ID
    if "transcript_id" in transcript_annot_dict:
        transcript_ID_val = transcript_annot_dict.pop("transcript_id")
    else:
        transcript_ID_val = transcript_ID
    attributes.append(make_descriptor_string("transcript_id", transcript_ID_val))

    # Mandatory: Gene Name
    if "gene_name" in gene_annot_dict:
        gene_name = gene_annot_dict.pop("gene_name")
    else:
        gene_name = "TALON_gene-" + str(gene_ID)
    attributes.append(make_descriptor_string("gene_name", gene_name))

    # Mandatory: Gene Status
    gene_status = gene_annot_dict.pop("gene_status")
    attributes.append(make_descriptor_string("gene_status", gene_status))

    # Gene Type
    if "gene_type" in gene_annot_dict:
        gene_type = gene_annot_dict.pop("gene_type")
        attributes.append(make_descriptor_string("gene_type", gene_type))
 
    # Transcript Type
    if "transcript_type" in transcript_annot_dict:
        transcript_type = transcript_annot_dict.pop("transcript_type")
        attributes.append(make_descriptor_string("transcript_type", transcript_type))

    # Mandatory: Transcript Status
    transcript_status = transcript_annot_dict.pop("transcript_status")
    attributes.append(make_descriptor_string("transcript_status", transcript_status))
    
    # Mandatory: Transcript Name
    if "transcript_name" in transcript_annot_dict:
        transcript_name = transcript_annot_dict.pop("transcript_name")
    else:
        transcript_name = "TALON_transcript-" + str(transcript_ID)
    attributes.append(make_descriptor_string("transcript_name", transcript_name))
    
    # TALON Gene ID
    attributes.append(make_descriptor_string("talon_gene", gene_ID))

    # TALON Transcript ID
    attributes.append(make_descriptor_string("talon_transcript", transcript_ID))

    # Add any remaining annotations
    for attribute,value in sorted(transcript_annot_dict.items()):
        attributes.append(make_descriptor_string(attribute, value)) 

    return attributes

def format_GTF_tag_values_for_exon(gene_ID, transcript_ID, exon_ID, exon_number,
                                   gene_annot_dict, transcript_annot_dict, 
                                   exon_annot_dict):
    """ Parses the annotations for this exon, and supplements them where
        necessary for novel exons """

    attributes = []

    # Mandatory: Gene ID
    if "gene_id" in gene_annot_dict:
        gene_ID_val = gene_annot_dict.pop("gene_id")
    else:
        gene_ID_val = gene_ID
    attributes.append(make_descriptor_string("gene_id", gene_ID_val))

    # Mandatory: Transcript ID
    if "transcript_id" in transcript_annot_dict:
        transcript_ID_val = transcript_annot_dict.pop("transcript_id")
    else:
        transcript_ID_val = transcript_ID
    attributes.append(make_descriptor_string("transcript_id", transcript_ID_val))

    # Gene Type
    if "gene_type" in gene_annot_dict:
        gene_type = gene_annot_dict.pop("gene_type")
        attributes.append(make_descriptor_string("gene_type", gene_type))

    # Mandatory: Gene Status
    gene_status = gene_annot_dict.pop("gene_status")
    attributes.append(make_descriptor_string("gene_status", gene_status))

    # Mandatory: Gene Name
    if "gene_name" in gene_annot_dict:
        gene_name = gene_annot_dict.pop("gene_name")
    else:
        gene_name = "TALON_gene-" + str(gene_ID)
    attributes.append(make_descriptor_string("gene_name", gene_name))

    # Transcript Type
    if "transcript_type" in transcript_annot_dict:
        transcript_type = transcript_annot_dict.pop("transcript_type")
        attributes.append(make_descriptor_string("transcript_type", transcript_type))

    # Mandatory: Transcript Status
    transcript_status = transcript_annot_dict.pop("transcript_status")
    attributes.append(make_descriptor_string("transcript_status", transcript_status))

    # Mandatory: Transcript Name
    if "transcript_name" in transcript_annot_dict:
        transcript_name = transcript_annot_dict.pop("transcript_name")
    else:
        transcript_name = "TALON_transcript-" + str(transcript_ID)
    attributes.append(make_descriptor_string("transcript_name", transcript_name))

    # Exon number
    attributes.append(make_descriptor_string("exon_number", exon_number))

    # Exon ID
    if "exon_id" in exon_annot_dict:
        exon_ID_val = exon_annot_dict.pop("exon_id")
    else:
        exon_ID_val = exon_ID
    attributes.append(make_descriptor_string("exon_id", exon_ID_val))

    # TALON Gene ID
    attributes.append(make_descriptor_string("talon_gene", gene_ID))

    # TALON Transcript ID
    attributes.append(make_descriptor_string("talon_transcript", transcript_ID))

    # TALON Exon ID
    attributes.append(make_descriptor_string("talon_exon", exon_ID))
    if "exon_number" in exon_annot_dict:
        exon_annot_dict.pop("exon_number")

    # Add any remaining annotations
    for attribute,value in sorted(exon_annot_dict.items()):
        attributes.append(make_descriptor_string(attribute, value))

    return attributes

def get_gene_GTF_entry(gene_ID, associated_transcript_tuples, annotation_dict):
    """ Creates a GTF annotation entry for the given gene """

    if "source" in annotation_dict:
        source = annotation_dict["source"]
    else:
        source = "TALON"

    # GTF fields
    chromosome = associated_transcript_tuples[0]["chromosome"]
    feature = "gene"
    start = str(associated_transcript_tuples[0]["min_pos"])
    end = str(associated_transcript_tuples[-1]["max_pos"])
    score = "."
    strand = associated_transcript_tuples[0]["strand"]
    frame = "."
    attributes = " ".join(format_GTF_tag_values_for_gene(gene_ID, annotation_dict))

    GTF = '\t'.join([chromosome, source, feature, start, end, score, strand, 
                     frame, attributes])
    return GTF


def get_transcript_GTF_entry(transcript_entry, curr_gene_annot_dict, curr_transcript_annot_dict):
    """ Creates a GTF annotation entry for the given transcript """

    if "source" in curr_transcript_annot_dict:
        source = curr_transcript_annot_dict["source"]
    else:
        source = "TALON"

    gene_ID = transcript_entry["gene_ID"]
    transcript_ID = transcript_entry["transcript_ID"]

    # GTF fields for transcript
    chromosome = str(transcript_entry["chromosome"])
    feature = "transcript"
    start = str(transcript_entry["min_pos"])
    end = str(transcript_entry["max_pos"])
    score = "."
    strand = transcript_entry["strand"]
    frame = "."
    attributes = " ".join(format_GTF_tag_values_for_transcript(gene_ID, 
                                                               transcript_ID, 
                                                               curr_gene_annot_dict,
                                                               curr_transcript_annot_dict))

    GTF = '\t'.join([chromosome, source, feature, start, end, score, strand,
                     frame, attributes])
    return GTF

def get_exon_GTF_entry(gene_ID, transcript_ID, exon_ID, exon_num, exon_ID_2_location, 
                       curr_gene_annot_dict, curr_transcript_annot_dict,
                       curr_exon_annot_dict):
    """ Creates a GTF annotation entry for the given exon """

    if "source" in curr_exon_annot_dict:
        source = curr_exon_annot_dict["source"]
    else:
        source = "TALON"

    curr_exon_location = exon_ID_2_location[exon_ID]
    chromosome = str(curr_exon_location[0])
    feature = "exon"
    start = str(curr_exon_location[1])
    end = str(curr_exon_location[2])
    score = "."
    strand = curr_exon_location[3]
    frame = "."
    attributes = " ".join(format_GTF_tag_values_for_exon(gene_ID, 
                                                         transcript_ID, 
                                                         exon_ID, exon_num,
                                                         curr_gene_annot_dict,
                                                         curr_transcript_annot_dict,
                                                         curr_exon_annot_dict))

    GTF = '\t'.join([chromosome, source, feature, start, end, score, strand, 
                     frame, attributes])
    return GTF

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
    build = options.build
    outfile = create_outname(options)


    if build == None:
        raise ValueError("Please provide a valid genome build name")
    if annot == None:
        raise ValueError("Please provide a valid annotation name")
    check_annot_validity(annot, database)

    # Make sure that the input database exists!
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)


    # Determine which transcripts to include
    transcript_whitelist = handle_filtering(options)

    create_gtf(database, annot, build, transcript_whitelist, outfile)    
    

if __name__ == '__main__':
    main()
