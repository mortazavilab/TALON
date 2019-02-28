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
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, os.pardir)))

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

    if len(datasets) > 0:
        col_query = col_query + ","

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Create a mapping to allow us to escape possible special characters in
    # dataset names
    dataset_id_mapper = {}
    new_dataset_names = []
    count = 1
    for dataset in datasets:
        name = "abd_" + str(count)
        dataset_id_mapper[name] = dataset
        new_dataset_names.append(name)
        count += 1

    dataset_cols = ",\n".join([ name + ".count AS '" + dataset_id_mapper[name] + "'" for name in new_dataset_names ])
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

    # Create an abundance query for every dataset  
    abundance_queries = "\n".join(["LEFT JOIN abundance AS " + x + \
                                   " ON t.transcript_ID = " + x + \
                                   ".transcript_ID AND " + x + \
                                   ".dataset = '" + dataset_id_mapper[x] \
                                   + "'" for x in new_dataset_names])

    # Combine the subparts of the query into one and run it
    full_query = "\n".join([col_query, dataset_cols, name_status_query, 
                            abundance_queries, whitelist_string])

    try:
        abundance_tuples = (cursor.execute(full_query)).fetchall()
    except Exception as e:
        print(e)
        raise RuntimeError("Something went wrong with the database query")

    conn.close() 
    return abundance_tuples

def write_abundance_file(abundances, datasets, outfile):
    """ Writes abundances to an output file """

    o = open(outfile, 'w')
    col_names = abundances[0].keys()
    o.write("\t".join(col_names) + "\n")

    abundance_list = [list(elem) for elem in abundances]

    # Find indices of columns that may need 'None' replaced
    annot_indices = [i for i, s in enumerate(col_names) if 'annot' in s]
    status_indices = [i for i, s in enumerate(col_names) if 'status' in s]
    dataset_indices = [i for i,s in enumerate(col_names) if s in set(datasets)]

    # Iterate over abundances, fixing Nones, and write to file
    for transcript in abundances:
        transcript = list(transcript)
        for index in annot_indices:
            if transcript[index] == None:
                transcript[index] = "NA"
        for index in status_indices:
            if transcript[index] == None:
                transcript[index] = "NOVEL"
        for index in dataset_indices:
            if transcript[index] == None:
                transcript[index] = 0
        o.write("\t".join([str(x) for x in transcript]) + "\n")
        
    o.close()
    return

def check_annot_validity(annot, database):
    """ Make sure that the user has entered a correct annotation name """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    cursor.execute("SELECT DISTINCT annot_name FROM gene_annotations")
    annotations = [str(x[0]) for x in cursor.fetchall()]
    conn.close()

    if "talon_run" in annotations:
        annotations.remove("talon_run") 

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
    outfile = create_outname(options)

    check_annot_validity(annot, database)

    # Determine which transcripts to include
    transcript_whitelist = handle_filtering(options)

    # Create the abundance file
    datasets = datasets = fetch_dataset_list(database)
    abundances = fetch_abundances(database, datasets, annot, transcript_whitelist)
    write_abundance_file(abundances, datasets, outfile)

if __name__ == '__main__':
    main()
