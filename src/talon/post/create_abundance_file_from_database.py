# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# create_abundance_file_from_database.py is a utility that outputs the abundance
# for each transcript in the TALON database across datasets. Modified by
# filtering option.

import sqlite3
import itertools
import operator
from optparse import OptionParser
from pathlib import Path

from . import filter_talon_transcripts as filt
from .. import dstruct as dstruct
from .. import length_utils as lu
from . import post_utils as putils
from . import ab_utils as autils
from .. import query_utils as qutils
from .. import talon as talon


def getOptions():
    parser = OptionParser()

    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")

    parser.add_option("--annot", "-a", dest = "annot",
        help = """Which annotation version to use. Will determine which
                  annotation transcripts are considered known or novel
                  relative to. Note: must be in the TALON database.""",
        type = "string")

    parser.add_option("--whitelist", dest = "whitelist",
                      help = "Whitelist file of transcripts to include in the \
                              output. First column should be TALON gene ID, \
                              second column should be TALON transcript ID",
                      metavar = "FILE", type = "string", default = None)

    parser.add_option("--build", "-b", dest = "build",
        help = "Genome build to use. Note: must be in the TALON database.",
        type = "string")

    parser.add_option("--datasets", "-d",  dest = "datasets_file",
        help = """Optional: A file indicating which datasets should be
                  included (one dataset name per line). Default is to include
                  all datasets.""",
        metavar = "FILE", type = "string", default = None)

    parser.add_option("--o", dest = "outprefix", help = "Prefix for output file",
        metavar = "FILE", type = "string")


    (options, args) = parser.parse_args()
    return options

def create_outname(options):
    """ Creates filename for the output abundance that reflects the input options that
        were used. """

    outname = options.outprefix + "_talon_abundance"
    if options.whitelist != None:
        outname = "_".join([ outname, "filtered" ])

    outname += ".tsv"
    return outname

# def fetch_dataset_list(dataset_file, database):
#     """ Gets a list of all datasets in the database """
#
#     conn = sqlite3.connect(database)
#     cursor = conn.cursor()
#     all_db_datasets = qutils.fetch_all_datasets(cursor)
#     conn.close()
#
#     if dataset_file == None:
#
#         return all_db_datasets
#
#     else:
#         datasets = []
#         with open(dataset_file, 'r') as f:
#             for line in f:
#                 dataset = line.strip()
#                 if dataset not in all_db_datasets:
#                     raise ValueError("Dataset name '%s' not found in database" \
#                                       % (dataset))
#                 datasets.append(dataset)
#
#         return datasets

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
           7) number of exons in transcript

        Returns a list of tuples (one tuple per transcript)
    """

    # datasets = fetch_dataset_list(database)
    abundance = create_abundance_dict(database, datasets)

    col_query = """SELECT
	               t.gene_ID,
	               t.transcript_ID,
                       ga_id.value AS annot_gene_id,
                       ta_id.value AS annot_transcript_id,
	               ga_name.value AS annot_gene_name,
	               ta_name.value AS annot_transcript_name,
                       t.n_exons"""

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
                """ % (annot, annot, annot, annot)

    full_query = "\n".join([col_query, name_status_query, whitelist_string])

    try:
        abundance_tuples = (cursor.execute(full_query)).fetchall()
        colnames = list(abundance_tuples[0].keys()) + list(datasets)
    except Exception as e:
        print(e)
        raise RuntimeError("Something went wrong with the database query")

    conn.close()

    # Now iterate over the query results to incorporate the abundance information
    final_abundance = []
    for entry in abundance_tuples:
        transcript_ID = entry["transcript_ID"]

        # limit only to expressed transcripts
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

def write_abundance_file(abundances, col_names, prefix, n_places, datasets,
                         novelty_types, transcript_lengths, outfile):
    """ Writes abundances and metadata to an output file """

    o = open(outfile, 'w')

    novelty_type_cols = ["gene_novelty", "transcript_novelty", "ISM_subtype"]

    first_dataset_index = len(col_names) - len(datasets)
    first_colnames = col_names[0:first_dataset_index]
    dataset_colnames = col_names[first_dataset_index:]
    all_colnames = first_colnames + ["length"] + novelty_type_cols + dataset_colnames
    o.write("\t".join(all_colnames) + "\n")

    abundance_list = [list(elem) for elem in abundances]

    # Find indices of columns that may need 'None' replaced
    gene_ID_index = all_colnames.index("gene_ID")
    transcript_ID_index = all_colnames.index("transcript_ID")
    annot_gene_ID_index = all_colnames.index("annot_gene_id")
    annot_transcript_ID_index = all_colnames.index("annot_transcript_id")
    gene_name_index = all_colnames.index("annot_gene_name")
    transcript_name_index = all_colnames.index("annot_transcript_name")
    dataset_indices = [i for i,s in enumerate(all_colnames) if s in set(datasets)]

    # Iterate over abundances, fixing Nones, and write to file
    for transcript in abundances:
        curr_novelty = get_gene_and_transcript_novelty_types(transcript[gene_ID_index],
                                                             transcript[transcript_ID_index],
                                                             novelty_types)
        transcript = list(transcript)
        transcript = transcript[0:first_dataset_index] + \
                     [transcript_lengths[transcript[transcript_ID_index]]] + \
                     [ curr_novelty[x] for x in novelty_type_cols] + \
                     transcript[first_dataset_index:]

        alt_gene_name, alt_transcript_name = talon.construct_names(transcript[gene_ID_index], \
                                                             transcript[transcript_ID_index], \
                                                             prefix, n_places)

        if transcript[annot_gene_ID_index] == None:
            transcript[annot_gene_ID_index] = alt_gene_name

        if transcript[gene_name_index] == None:
            transcript[gene_name_index] = alt_gene_name

        if transcript[annot_transcript_ID_index] == None:
            transcript[annot_transcript_ID_index] = alt_transcript_name

        if transcript[transcript_name_index] == None:
            transcript[transcript_name_index] = alt_transcript_name

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

    # Look for gene type
    if gene_ID in novelty_type.antisense_genes:
        curr_novel["gene_novelty"] = "Antisense"
    elif gene_ID in novelty_type.intergenic_genes:
        curr_novel["gene_novelty"] = "Intergenic"
    elif gene_ID in novelty_type.known_genes:
        curr_novel["gene_novelty"] = "Known"
    else:
        print("Warning: Could not locate novelty type for gene %s" % gene_ID)

    # Look for transcript type
    if transcript_ID in novelty_type.ISM_transcripts:
        curr_novel["transcript_novelty"] = "ISM"
    elif transcript_ID in novelty_type.NIC_transcripts:
        curr_novel["transcript_novelty"] = "NIC"
    elif transcript_ID in novelty_type.NNC_transcripts:
        curr_novel["transcript_novelty"] = "NNC"
    elif transcript_ID in novelty_type.antisense_transcripts:
        curr_novel["transcript_novelty"] = "Antisense"
    elif transcript_ID in novelty_type.intergenic_transcripts:
        curr_novel["transcript_novelty"] = "Intergenic"
    elif transcript_ID in novelty_type.genomic_transcripts:
        curr_novel["transcript_novelty"] = "Genomic"
    elif transcript_ID in novelty_type.known_transcripts:
        curr_novel["transcript_novelty"] = "Known"
    else:
        print("Warning: Could not locate novelty type for transcript %s" % transcript_ID)

    # Look for ISM subtype
    if transcript_ID in novelty_type.ISM_prefix and \
       transcript_ID in novelty_type.ISM_suffix:
        curr_novel["ISM_subtype"] = "Both"
    elif transcript_ID in novelty_type.ISM_prefix:
        curr_novel["ISM_subtype"] = "Prefix"
    elif transcript_ID in novelty_type.ISM_suffix:
        curr_novel["ISM_subtype"] = "Suffix"
    else:
        curr_novel["ISM_subtype"] = "None"

    return curr_novel

# def check_annot_validity(annot, database):
#     """ Make sure that the user has entered a correct annotation name """
#
#     conn = sqlite3.connect(database)
#     cursor = conn.cursor()
#
#     cursor.execute("SELECT DISTINCT annot_name FROM gene_annotations")
#     annotations = [str(x[0]) for x in cursor.fetchall()]
#     conn.close()
#
#     if "TALON" in annotations:
#         annotations.remove("TALON")
#
#     if annot == None:
#         message = "Please provide a valid annotation name. " + \
#                   "In this database, your options are: " + \
#                   ", ".join(annotations)
#         raise ValueError(message)
#
#     if annot not in annotations:
#         message = "Annotation name '" + annot + \
#                   "' not found in this database. Try one of the following: " + \
#                   ", ".join(annotations)
#         raise ValueError(message)
#
#     return
#
# def check_build_validity(build, database):
#     """ Make sure that the user has entered a correct build name """
#
#     conn = sqlite3.connect(database)
#     cursor = conn.cursor()
#
#     cursor.execute("SELECT name FROM genome_build")
#     builds = [str(x[0]) for x in cursor.fetchall()]
#     conn.close()
#
#     if build == None:
#         message = "Please provide a valid genome build name. " + \
#                   "In this database, your options are: " + \
#                   ", ".join(builds)
#         raise ValueError(message)
#
#     if build not in builds:
#         message = "Build name '" + build + \
#                   "' not found in this database. Try one of the following: " + \
#                   ", ".join(builds)
#         raise ValueError(message)
#
#     return

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

# def fetch_naming_prefix(database):
#     """ Get naming prefix from the database run_info table """
#     conn = sqlite3.connect(database)
#     conn.row_factory = sqlite3.Row
#     cursor = conn.cursor()
#     cursor.execute("SELECT value FROM run_info WHERE item = 'idprefix'")
#     prefix = cursor.fetchone()[0]
#
#     conn.close()
#     return prefix
#
# def fetch_n_places(database):
#     """ Get length of name field from the database run_info table """
#     conn = sqlite3.connect(database)
#     conn.row_factory = sqlite3.Row
#     cursor = conn.cursor()
#     cursor.execute("SELECT value FROM run_info WHERE item = 'n_places'")
#     n_places = cursor.fetchone()[0]
#
#     conn.close()
#     return int(n_places)

# def get_transcript_lengths(database, build):
#     """ Read the transcripts from the database. Then compute the lengths.
#         Store in a dictionary """
#
#     transcript_lengths = {}
#
#     conn = sqlite3.connect(database)
#     conn.row_factory = sqlite3.Row
#     cursor = conn.cursor()
#
#     # Get the exon lengths
#     exon_lens = lu.get_all_exon_lengths(cursor, build)
#
#     cursor.execute("SELECT * FROM transcripts")
#     for transcript_row in cursor.fetchall():
#         transcript_ID = transcript_row['transcript_ID']
#         length = lu.get_transcript_length(transcript_row, exon_lens)
#         transcript_lengths[transcript_ID] = length
#
#     conn.close()
#     return transcript_lengths


def main():
    options = getOptions()
    database = options.database
    annot = options.annot
    build = options.build

    whitelist_file = options.whitelist
    dataset_file = options.datasets_file
    outfile = create_outname(options)

    # Make sure that the input database exists!
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)

    autils.check_annot_validity(annot, database)
    autils.check_build_validity(build, database)

    # Determine which transcripts to include
    whitelist = putils.handle_filtering(database,
                                        annot,
                                        False,
                                        whitelist_file,
                                        dataset_file)

    # create transcript whitelist
    transcript_whitelist = []
    for key,group in itertools.groupby(whitelist,operator.itemgetter(0)):
        for id_tuple in list(group):
            transcript_whitelist.append(str(id_tuple[1]))

    # Get transcript length dict
    transcript_lengths = autils.get_transcript_lengths(database, build)

    # Create the abundance file
    datasets = autils.fetch_dataset_list(dataset_file, database)
    novelty_type = make_novelty_type_struct(database, datasets)
    abundances, colnames = fetch_abundances(database, datasets, annot, transcript_whitelist)
    prefix = autils.fetch_naming_prefix(database)
    n_places = autils.fetch_n_places(database)
    write_abundance_file(abundances, colnames, prefix, n_places, datasets, novelty_type, transcript_lengths, outfile)

if __name__ == '__main__':
    main()
