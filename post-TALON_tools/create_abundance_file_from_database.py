# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# create_abundance_file_from_database.py is a utility that outputs the abundance 
# for each transcript in the TALON database across datasets. Modified by 
# filtering option.

from optparse import OptionParser
import sqlite3
import sys
import os
import filter_talon_transcripts as filt
from pathlib import Path
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, os.pardir)))
main_path = "/".join(script_dir.split("/")[0:-1])
sys.path.append(main_path)
import dstruct as dstruct
import query_utils as qutils

def getOptions():
    parser = OptionParser()

    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")

    parser.add_option("--annot", "-a", dest = "annot",
        help = """Which annotation version to use. Will determine which
                  annotation transcripts are considered known or novel
                  relative to. Note: must be in the TALON database.""",
        type = "string")

    parser.add_option("--filter", dest ="filtering", action='store_true',
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

    parser.add_option("--o", dest = "outprefix", help = "Prefix for output file",
        metavar = "FILE", type = "string")


    (options, args) = parser.parse_args()
    return options

def handle_filtering(options):
    """ Determines which transcripts to allow in the abundance file. This can be done
        in two different ways. If filtering is turned off, then all of the
        transcripts in the database go on the whitelist. 
        If filtering is turned on, known transcripts and novel
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

    # Otherwise, the whitelist will be every transcript ID in the database
    else:

        conn = sqlite3.connect(database)
        cursor = conn.cursor()

        cursor.execute("SELECT gene_ID,transcript_ID FROM transcripts")
        whitelist = sorted(cursor.fetchall())

        conn.close()

    # Sort the whitelist transcript IDs
    whitelist = [str(x[1]) for x in whitelist]
    sorted_whitelist = sorted(whitelist)

    return sorted_whitelist

def create_outname(options):
    """ Creates filename for the output GTF that reflects the input options that
        were used. """

    outname = options.outprefix + "_talon_abundance"
    if options.filtering == True:
        outname = "_".join([ outname, "filtered" ])

    outname += ".tsv"
    return outname

def fetch_dataset_list(database):
    """ Gets a list of all datasets in the database """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    cursor.execute("SELECT dataset_name FROM dataset")
    datasets = [str(x[0]) for x in cursor.fetchall()]
   
    conn.close()
    return datasets

def create_abundance_dict(database, datasets):
    """Process the abundance table by dataset in order to create a dictionary
       data structure organized like this:
           transcript_ID -> dataset -> abundance in that dataset
    """
    abundance = {}    

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    for dataset in datasets:
        query = """ SELECT transcript_ID, count FROM abundance 
                    WHERE dataset = '%s' """ % dataset
        cursor.execute(query)

        for transcript in cursor.fetchall():
            transcript_ID = transcript["transcript_ID"]
            count = transcript["count"]  

            if transcript_ID in abundance:
                abundance[transcript_ID][dataset] = count
            else:
                abundance[transcript_ID] = {}
                abundance[transcript_ID][dataset] = count 

    conn.close()
    return abundance

def fetch_abundances(database, datasets, annot, whitelist):
    """Constructs a query to get the following information for every 
       whitelisted transcript:
           1) TALON gene ID
           2) TALON transcript ID
           3) Gene ID (from annotation specified in 'annot', None otherwise)
           4) Transcript ID (from annotation specified in 'annot', None otherwise)
           5) Gene name (from annotation specified in 'annot', None otherwise)
           6) Transcript name (from annotation specified in 'annot', None otherwise)
           7) Gene annotation status (KNOWN/NOVEL)
           8) Transcript annotation status (KNOWN/NOVEL)
           9) number of exons in transcript
           10+) Count of this transcript in each dataset in the database

        Returns a list of tuples (one tuple per transcript)
    """

    datasets = fetch_dataset_list(database)
    abundance = create_abundance_dict(database, datasets)

    col_query = """SELECT
	               t.gene_ID,
	               t.transcript_ID,
                       ga_id.value AS annot_gene_id,
                       ta_id.value AS annot_transcript_id,
	               ga_name.value AS annot_gene_name,
	               ta_name.value AS annot_transcript_name,
                       t.n_exons,
	               ga_status.value AS gene_status,
	               ta_status.value AS transcript_status"""

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    whitelist_string = "WHERE t.transcript_ID IN (" + ','.join(whitelist) + ");"

    name_status_query = """        
                FROM transcripts t
                LEFT JOIN gene_annotations ga_id ON t.gene_ID = ga_id.ID
                    AND ga_id.annot_name = '%s'
                    AND ga_id.attribute = 'gene_id'
                LEFT JOIN transcript_annotations ta_id ON t.transcript_ID = ta_id.ID
                    AND ta_id.annot_name = '%s'
                    AND ta_id.attribute = 'transcript_id'
                LEFT JOIN gene_annotations ga_name ON t.gene_ID = ga_name.ID 
	            AND ga_name.annot_name = '%s' 
                    AND ga_name.attribute = 'gene_name'
                LEFT JOIN transcript_annotations ta_name ON t.transcript_ID = ta_name.ID 
	            AND ta_name.annot_name = '%s' 
                    AND ta_name.attribute = 'transcript_name' 
                LEFT JOIN gene_annotations ga_status ON t.gene_ID = ga_status.ID 
	            AND ga_status.annot_name = '%s' 
                    AND ga_status.attribute = 'gene_status' 
                LEFT JOIN transcript_annotations ta_status ON t.transcript_ID = ta_status.ID 
	            AND ta_status.annot_name = '%s' 
                    AND ta_status.attribute = 'transcript_status'
                """ % (annot, annot, annot, annot, annot, annot)

    full_query = "\n".join([col_query, name_status_query, whitelist_string])

    try:
        abundance_tuples = (cursor.execute(full_query)).fetchall()
        colnames = list(abundance_tuples[0].keys()) + datasets
    except Exception as e:
        print(e)
        raise RuntimeError("Something went wrong with the database query")

    conn.close()

    # Now iterate over the query results to incorporate the abundance information
    final_abundance = []
    for entry in abundance_tuples:
        transcript_ID = entry["transcript_ID"]
        if transcript_ID not in abundance:
            continue

        # Get abundance of this transcript in each dataset
        dataset_counts = []
        for dataset in datasets:
            if dataset in abundance[transcript_ID]:
                dataset_counts.append(str(abundance[transcript_ID][dataset]))
            else: 
                dataset_counts.append("0")

        # Combine abundance info with rest of transcript information        
        combined_entry = list(entry) + dataset_counts
        final_abundance.append(combined_entry)

    return final_abundance, colnames

def write_abundance_file(abundances, col_names, prefix, datasets, novelty_types, outfile):
    """ Writes abundances to an output file """

    o = open(outfile, 'w')
    
    novelty_type_cols = ["antisense_gene", "intergenic_gene", "ISM_transcript",
                         "ISM_prefix_transcript", "ISM_suffix_transcript",
                         "NIC_transcript", "NNC_transcript", "antisense_transcript",
                         "intergenic_transcript", "genomic_transcript"]

    first_dataset_index = len(col_names) - len(datasets)
    first_colnames = col_names[0:first_dataset_index]
    dataset_colnames = col_names[first_dataset_index:]
    all_colnames = first_colnames + novelty_type_cols + dataset_colnames
    o.write("\t".join(all_colnames) + "\n")

    abundance_list = [list(elem) for elem in abundances]
    
    # Find indices of columns that may need 'None' replaced
    gene_ID_index = all_colnames.index("gene_ID")
    transcript_ID_index = all_colnames.index("transcript_ID")
    gene_name_index = all_colnames.index("annot_gene_name")
    transcript_name_index = all_colnames.index("annot_transcript_name")
    status_indices = [i for i, s in enumerate(all_colnames) if 'status' in s]
    dataset_indices = [i for i,s in enumerate(all_colnames) if s in set(datasets)]

    # Iterate over abundances, fixing Nones, and write to file
    for transcript in abundances:
        curr_novelty = get_gene_and_transcript_novelty_types(transcript[gene_ID_index], 
                                                             transcript[transcript_ID_index], 
                                                             novelty_types)
        transcript = list(transcript)
        transcript = transcript[0:first_dataset_index] + \
                     [ curr_novelty[x] for x in novelty_type_cols] + \
                     transcript[first_dataset_index:]

        if transcript[gene_name_index] == None:
            transcript[gene_name_index] = prefix + "-gene_" + \
                                          str(transcript[gene_ID_index])
        if transcript[transcript_name_index] == None:
            transcript[transcript_name_index] = prefix + "-transcript_" + \
                                            str(transcript[transcript_ID_index])
        for index in status_indices:
            if transcript[index] == None:
                transcript[index] = "NOVEL"
        for index in dataset_indices:
            if transcript[index] == None:
                transcript[index] = 0
        o.write("\t".join([str(x) for x in transcript]) + "\n")
        
    o.close()
    return

def get_gene_and_transcript_novelty_types(gene_ID, transcript_ID, novelty_type):
    """ Look up gene and transcript IDs in data structure to determine which types
        of novelty are present """

    curr_novel = {}
    curr_novel["antisense_gene"] = "antisense_gene" \
               if gene_ID in novelty_type.antisense_genes else "No"
    curr_novel["intergenic_gene"] = "intergenic_gene" \
               if gene_ID in novelty_type.intergenic_genes else "No"
    curr_novel["ISM_transcript"] = "ISM_transcript" \
               if transcript_ID in novelty_type.ISM_transcripts else "No"
    curr_novel["ISM_prefix_transcript"] = "ISM_prefix_transcript" \
               if transcript_ID in novelty_type.ISM_prefix else "No"
    curr_novel["ISM_suffix_transcript"] = "ISM_suffix_transcript" \
               if transcript_ID in novelty_type.ISM_suffix else "No"
    curr_novel["NIC_transcript"] = "NIC_transcript" \
               if transcript_ID in novelty_type.NIC_transcripts else "No"
    curr_novel["NNC_transcript"] = "NNC_transcript" \
               if transcript_ID in novelty_type.NNC_transcripts else "No"
    curr_novel["antisense_transcript"] = "antisense_transcript" \
               if transcript_ID in novelty_type.antisense_transcripts else "No"
    curr_novel["intergenic_transcript"] = "intergenic_transcript" \
               if transcript_ID in novelty_type.intergenic_transcripts else "No"
    curr_novel["genomic_transcript"] = "genomic_transcript" \
               if transcript_ID in novelty_type.genomic_transcripts else "No"

    return curr_novel

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

def make_novelty_type_struct(database, datasets):
    """ Create a data structure where it is possible to look up whether a gene
        or transcript belongs to a particular category of novelty"""

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row 
    cursor = conn.cursor() 

    novelty_type = dstruct.Struct()
    novelty_type.known_genes = set(qutils.fetch_all_known_genes_detected(cursor, datasets))
    novelty_type.antisense_genes = set(qutils.fetch_antisense_genes(cursor, datasets))
    novelty_type.intergenic_genes = set(qutils.fetch_intergenic_novel_genes(cursor, datasets))
    novelty_type.known_transcripts = set(qutils.fetch_all_known_transcripts_detected(cursor, datasets))
    novelty_type.ISM_transcripts = set(qutils.fetch_all_ISM_transcripts(cursor, datasets))
    novelty_type.ISM_prefix = set(qutils.fetch_prefix_ISM_transcripts(cursor, datasets))
    novelty_type.ISM_suffix = set(qutils.fetch_suffix_ISM_transcripts(cursor, datasets))
    novelty_type.NIC_transcripts = set(qutils.fetch_NIC_transcripts(cursor, datasets))
    novelty_type.NNC_transcripts = set(qutils.fetch_NNC_transcripts(cursor, datasets))
    novelty_type.antisense_transcripts = set(qutils.fetch_antisense_transcripts(cursor, datasets))
    novelty_type.intergenic_transcripts = set(qutils.fetch_intergenic_transcripts(cursor, datasets))
    novelty_type.genomic_transcripts = set(qutils.fetch_genomic_transcripts(cursor, datasets))

    conn.close()
    return novelty_type

def fetch_naming_prefix(database):
    """ Get naming prefix from the database run_info table """
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    cursor.execute("SELECT value FROM run_info WHERE item = 'idprefix'")
    prefix = cursor.fetchone()[0]

    conn.close()
    return prefix

def main():
    options = getOptions()
    database = options.database
    annot = options.annot
    outfile = create_outname(options)

    check_annot_validity(annot, database)

    # Make sure that the input database exists!
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)

    # Determine which transcripts to include
    transcript_whitelist = handle_filtering(options)

    # Create the abundance file
    datasets = datasets = fetch_dataset_list(database)
    novelty_type = make_novelty_type_struct(database, datasets)
    abundances, colnames = fetch_abundances(database, datasets, annot, transcript_whitelist)
    prefix = fetch_naming_prefix(database)
    write_abundance_file(abundances, colnames, prefix, datasets, novelty_type, outfile)

if __name__ == '__main__':
    main()
