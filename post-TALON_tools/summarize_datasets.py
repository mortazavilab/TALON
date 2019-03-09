import argparse
import sys
import sqlite3
import os
from pathlib import Path
script_path = os.path.abspath(__file__)
main_path = "/".join(script_path.split("/")[0:-2])
sys.path.append(main_path)
import query_utils as qutils

def get_args():
    """ Fetches the arguments for the program """

    program_desc = """Generates a tab-delimited file of gene and transcript 
                      counts for each dataset in the database (broken down
                      by category)."""
    parser = argparse.ArgumentParser(description=program_desc)
    parser.add_argument('--db', dest = 'database', metavar='FILE,', type = str,
        help='TALON database')
    parser.add_argument("--verbose", 
                        help = "Verbose mode: print out the counts in terminal",
                        action="store_true")
    parser.add_argument("--o", dest = "outprefix", 
                        help = "Prefix for output file", type = str)

    args = parser.parse_args()
    return args

def write_counts_file(cursor, outprefix, verbose = False):
    """ Create a log file with the following columns:
            - dataset name
            - Number of reads annotated
            - Number of known genes detected (total)
            - Number of novel genes detected (total)
            - Number of known transcripts detected (total)
            - Number of novel transcripts detected (total)
            Breakdowns by category
            - Number of antisense genes detected
            - Number of intergenic genes detected
            - Number of known transcripts
            - Number of FSM transcripts detected (perfect + with novelty)
            - Number of total ISM transcripts detected
            - Number of suffix ISMs detected
            - Number of antisense transcripts detected
            - Number of genomic transcripts detected
    """
    o = open(outprefix + "_talon_summary.tsv", 'w')
    columns = [ "dataset", "reads_annotated", "known_genes", "antisense_genes",
                "other_novel_genes", "known_transcripts", "novel_transcripts",
                 "ISMs", "prefix_ISMs", "suffix_ISMs", "NICs", "NNCs",
                "antisense_transcripts", "genomic_transcripts" ]
    o.write("\t".join(columns) + "\n")

    # Get dataset names
    cursor.execute(""" SELECT dataset_name FROM dataset """)
    datasets = [ str(x[0]) for x in cursor.fetchall() ]
    for dataset in datasets:
        # Get number of reads in the dataset
        reads = qutils.count_observed_reads(cursor, dataset)

        # Get the number of known genes detected
        known_genes = qutils.count_known_genes_detected(cursor, dataset)

        # Get the number of novel genes detected
        novel_genes = qutils.count_novel_genes_detected(cursor, dataset)

        # Get the number of known transcripts detected
        known_transcripts = len(qutils.fetch_all_known_transcripts_detected(cursor, dataset))

        # Get the number of novel transcripts
        novel_transcripts = len(qutils.fetch_novel_transcripts(cursor, dataset))

        # Get antisense genes
        antisense_genes = len(qutils.fetch_antisense_genes(cursor, dataset))

        # Get intergenic genes
        intergenic_genes = len(qutils.fetch_intergenic_novel_genes(cursor, dataset))

        # Get ISM transcripts
        ISMs = len(qutils.fetch_all_ISM_transcripts(cursor, dataset))
        prefix_ISMs = len(qutils.fetch_prefix_ISM_transcripts(cursor, dataset))
        suffix_ISMs = len(qutils.fetch_suffix_ISM_transcripts(cursor, dataset))

        # Get NIC transcripts
        NICs = len(qutils.fetch_NIC_transcripts(cursor, dataset))

        # Get NNC transcripts
        NNCs = len(qutils.fetch_NNC_transcripts(cursor, dataset))

        # Get antisense transcripts
        antisense_transcripts = len(qutils.fetch_antisense_transcripts(cursor, dataset))

        # Get genomic novel transcripts
        genomic_transcripts = len(qutils.fetch_genomic_transcripts(cursor, dataset))

        outputs = [ dataset, reads, known_genes, antisense_genes,
                intergenic_genes, known_transcripts, ISMs, prefix_ISMs,
                suffix_ISMs, NICs, NNCs, antisense_transcripts, 
                genomic_transcripts ]

        if verbose == True:
            print("---------------%s---------------" % dataset)
            print("Number of reads: %d" % reads)
            print("Known genes: %d" % known_genes)
            print("Novel genes: %d" % novel_genes)
            print("----Antisense genes: %d" % antisense_genes)
            print("----Other novel genes: %d" % intergenic_genes)
            print("Known transcripts: %d" % known_transcripts)
            print("Novel transcripts: %d" % novel_transcripts)
            print("----ISM transcripts: %d" % ISMs)
            print("--------Prefix ISM transcripts: %d" % prefix_ISMs)
            print("--------Suffix ISM transcripts: %d" % suffix_ISMs)
            print("----NIC transcripts: %d" % NICs)
            print("----NNC transcripts: %d" % NNCs)
            print("----antisense transcripts: %d" % antisense_transcripts)
            print("----genomic transcripts: %d" % genomic_transcripts)

        o.write("\t".join([str(x) for x in outputs]) + "\n")

    o.close()

def main():
    options = get_args()

    # Make sure that the input database exists!
    database = options.database
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    write_counts_file(cursor, options.outprefix, options.verbose) 
    conn.close()
    

if __name__ == '__main__':
    main()
