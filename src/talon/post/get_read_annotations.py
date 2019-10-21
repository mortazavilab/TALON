# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# get_read_annotations.py is a utility that queries a TALON 
# database in order to get read-specific annotation information.

import argparse
import sqlite3
import os
from .. import query_utils as qutils
from string import Template

def get_args():
    """ Fetches the arguments for the program """

    program_desc = ("This utility queries a TALON database in order to get "
                    "read-specific annotation information.")
    parser = argparse.ArgumentParser(description=program_desc)

    parser.add_argument('--db', dest = 'database', metavar='FILE,', type = str,
        help='TALON database')
    parser.add_argument('--build', dest = 'build', metavar='STRING,', type = str,
        help='Genome build (i.e. hg38) to use. Must be in the database.')
    parser.add_argument('--datasets', dest = 'datasets', metavar='STRING,', type = str,
        help=('Optional: Comma-delimited list of datasets to include. Default '
              'behavior is to include all datasets in the database.'),
        default = "all")
    parser.add_argument("--o", dest = "outprefix", help = "Prefix for output files",
        type = str)

    args = parser.parse_args()
    return args

def fetch_reads(database, build, tmp_file = None, datasets = "all"):
    """ Performs database query to fetch location and gene/transcript assignment
        info for each long read in the specified datasets. 
        If tmp_file is set to None (default), then the function will return
        the query results in a list of lists. If an alternate value is provided, 
        then the results will be written to a tmp file of that name."""

    if datasets != "all":
        # Format as a string for query
        dataset_str = qutils.format_for_IN(datasets)
        dataset_str = " AND dataset IN " + dataset_str
    else:
        dataset_str = ""

    with sqlite3.connect(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        query = """ SELECT os.read_name, 
                           os.dataset, 
                           loc1.genome_build,
                           os.gene_ID as gene_ID,
                           os.transcript_ID as transcript_ID,
                           loc1.chromosome as chrom,
                           loc1.position as start_vertex_pos,
                           loc2.position as end_vertex_pos,
                           genes.strand,
                           transcripts.n_exons,
                           os.read_length,
                           os.start_delta as TSS_diff,
                           os.end_delta as TTS_diff
                    FROM observed as os
                    LEFT JOIN location as loc1 ON 
                        loc1.location_ID = os.start_vertex 
                    LEFT JOIN location as loc2 ON 
                        loc2.location_ID = os.end_vertex
                    LEFT JOIN genes ON genes.gene_ID = os.gene_ID
                    LEFT JOIN transcripts ON 
                        transcripts.transcript_ID = os.transcript_ID
                    WHERE loc1.genome_build = '$build' 
                    AND loc2.genome_build = '$build' """
        query = Template(query + dataset_str)
        try:
            cursor.execute(query.substitute({"build": build}))
        except Exception as e:
            print(e)
            raise RuntimeError("Problem with reads database query")

        if tmp_file != None:
            o = open(tmp_file, 'w') 
        else:
            reads = []

        colnames = []

        count = 0
        for entry in cursor:
            strand = entry["strand"]
            if entry["TSS_diff"] == None:
                TSS_diff = 0
            else:
                TSS_diff = entry["TSS_diff"]

            if entry["TTS_diff"] == None:
                TTS_diff = 0
            else:
                TTS_diff = entry["TTS_diff"]

            if strand == "+":
                read_start = entry["start_vertex_pos"] + TSS_diff
                read_end = entry["end_vertex_pos"] + TTS_diff
            elif strand == "-":
                read_start = entry["start_vertex_pos"] - TSS_diff
                read_end = entry["end_vertex_pos"] - TTS_diff
            else:
                raise ValueError("Unrecognized strand value: " + str(strand))
            
            # Create entry for output
            out_read = (entry["read_name"], entry["dataset"],
                        entry["genome_build"], entry["gene_ID"],
                        entry["transcript_ID"], entry["chrom"], 
                        read_start, read_end, strand, entry["n_exons"],
                        entry["read_length"])
           
            if tmp_file != None:
                o.write("\t".join([ str(x) for x in out_read ]) + "\n")
            else:
                reads.append(out_read)
            count += 1

    # Return results or close file
    if count == 0:
        raise ValueError(("No reads detected. Make sure your dataset names are " 
                          "correct."))

    if tmp_file != None:
        o.close()
    else:
        return reads

def get_gene_novelty(database):
    """ Given a database, get the novelty status of each gene. """

    gene_novelty = {}
    with sqlite3.connect(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Fetch known genes
        cursor.execute("""SELECT ID FROM gene_annotations
                              WHERE attribute = "gene_status"
                              AND value = "KNOWN";""")
        for entry in cursor:
            gene_novelty[entry[0]] = "Known"

        # Fetch antisense genes
        cursor.execute("""SELECT ID FROM gene_annotations
                              WHERE attribute = "antisense_gene"
                              AND value = "TRUE";""")
        for entry in cursor:
            gene_novelty[entry[0]] = "Antisense"

        # Fetch intergenic genes
        cursor.execute("""SELECT ID FROM gene_annotations
                              WHERE attribute = "intergenic_novel"
                              AND value = "TRUE";""")
        for entry in cursor:
            gene_novelty[entry[0]] = "Intergenic"

    return gene_novelty

def get_transcript_novelty(database):
    """ Given a database, get the novelty status of each transcript. """
  
    transcript_novelty = {}
    with sqlite3.connect(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Fetch known transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "transcript_status"
                              AND value = "KNOWN";""")
        for entry in cursor:
            transcript_novelty[entry[0]] = "Known"
    
        # Fetch ISM transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "ISM_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            transcript_novelty[entry[0]] = "ISM"
    
        # Fetch NIC transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "NIC_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            transcript_novelty[entry[0]] = "NIC"
    
        # Fetch NNC transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "NNC_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            transcript_novelty[entry[0]] = "NNC"
    
        # Fetch antisense transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "antisense_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            transcript_novelty[entry[0]] = "Antisense"
    
        # Fetch intergenic transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "intergenic_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            transcript_novelty[entry[0]] = "Intergenic"
    
        # Fetch genomic transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "genomic_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            transcript_novelty[entry[0]] = "Genomic"
    
    return transcript_novelty

def get_ISM_novelty(database):
    """ Given a database, get the ISM subtype of each ISM transcript. """
    
    all_ISMs = set()
    prefix_ISMs = set()
    suffix_ISMs = set()
    ISM_novelty = {}
    
    with sqlite3.connect(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Fetch ISM transcripts
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "ISM_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            all_ISMs.add(entry[0])

        # Fetch Prefix ISMs
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "ISM-prefix_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            prefix_ISMs.add(entry[0])

        # Fetch Suffix ISMs
        cursor.execute("""SELECT ID FROM transcript_annotations
                              WHERE attribute = "ISM-suffix_transcript"
                              AND value = "TRUE";""")
        for entry in cursor:
            suffix_ISMs.add(entry[0])

    # Look for ISM subtype
    for transcript_ID in all_ISMs:
        if transcript_ID in prefix_ISMs and transcript_ID in suffix_ISMs:
            ISM_novelty[transcript_ID] = "Both"
        elif transcript_ID in prefix_ISMs:
            ISM_novelty[transcript_ID] = "Prefix"
        elif transcript_ID in suffix_ISMs:
            ISM_novelty[transcript_ID] = "Suffix"
        else:
            ISM_novelty[transcript_ID] = "None"

    return ISM_novelty

def get_gene_annotations(database): 
    """ Create a dictionary linking each TALON gene ID to its human-readable
        name and accession ID """

    gene_name = {}
    gene_ID = {}

    with sqlite3.connect(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        cursor.execute("""SELECT ID, ga.value FROM gene_annotations as ga 
                          WHERE attribute = "gene_name";""")
        for entry in cursor:
            gene_name[entry["ID"]] = entry["value"]

        cursor.execute("""SELECT ID, ga.value FROM gene_annotations as ga
                          WHERE attribute = "gene_id";""") 
        for entry in cursor:
            gene_ID[entry["ID"]] = entry["value"]

    return gene_name, gene_ID

def get_transcript_annotations(database):
    """ Create a dictionary linking each TALON transcript ID to its human-readable
        name and accession ID """

    transcript_name = {}
    transcript_ID = {}

    with sqlite3.connect(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        cursor.execute("""SELECT ID, ta.value FROM transcript_annotations as ta
                          WHERE attribute = "transcript_name";""")
        for entry in cursor:
            transcript_name[entry["ID"]] = entry["value"]

        cursor.execute("""SELECT ID, ta.value FROM transcript_annotations as ta
                          WHERE attribute = "transcript_id";""")
        for entry in cursor:
            transcript_ID[entry["ID"]] = entry["value"]

    return transcript_name, transcript_ID

def make_read_annot_file(database, build, outprefix, datasets = "all"):
    """ Creates an output file with the following columns:
            1. read_name
            2. dataset
            3. genome_build
            4. chrom
            5. read_start
            6. read_end
            7. strand
            8. n_exons
            9. read_length
            10. gene_ID (TALON)
            11. transcript_ID (TALON)
            12. annot_gene_id
            13. annot_transcript_id
            14. annot_gene_name
            15. annot_transcript_name
            16. gene_novelty
            17. transcript_novelty
            18. ISM_subtype

        By default, reads from all datasets in the database are included, but 
        this can be modified by supplying a list/tuple of dataset names to the 
        datasets parameter.
    """
    tmp_read_file = outprefix + "_reads.tmp"
    fetch_reads(database, build, tmp_file = tmp_read_file, datasets = datasets)

    # Make annotation dicts
    gene_names, gene_IDs = get_gene_annotations(database)
    transcript_names, transcript_IDs = get_transcript_annotations(database) 

    # Make novelty dicts
    gene_novelty = get_gene_novelty(database)
    transcript_novelty = get_transcript_novelty(database)
    ISM_novelty = get_ISM_novelty(database) 

    fname = outprefix + "_talon_read_annot.tsv"
    o = open(fname, 'w')
    colnames = [ "read_name", "dataset", "genome_build", "chrom", 
                 "read_start", "read_end", "strand", "n_exons", "read_length",
                 "gene_ID", "transcript_ID", "annot_gene_id", "annot_transcript_id",
                 "annot_gene_name", "annot_transcript_name", "gene_novelty", 
                 "transcript_novelty", "ISM_subtype"]
    o.write("\t".join(colnames) + "\n")

    with open(tmp_read_file, 'r') as f:
        for read_entry in f:
            read_name, dataset, genome_build, gene_ID, \
            transcript_ID, chrom, read_start, read_end, \
            strand, n_exons, read_length = read_entry.split("\t")

            gene_ID = int(gene_ID)
            transcript_ID = int(transcript_ID)
            # Get novelty info
            curr_gene_novelty = gene_novelty[gene_ID]
            curr_transcript_novelty = transcript_novelty[transcript_ID]

            if curr_transcript_novelty == "ISM":
                curr_ISM_novelty = ISM_novelty[transcript_ID]
            else:
                curr_ISM_novelty = "None"    
             
            # Get annotation info
            annot_gene_id = gene_IDs[gene_ID]
            annot_gene_name = gene_names[gene_ID]
            annot_transcript_id = transcript_IDs[transcript_ID]
            annot_transcript_name = transcript_names[transcript_ID]

            gene_ID = str(gene_ID)
            transcript_ID = str(transcript_ID)
            o.write("\t".join([read_name, dataset, genome_build, chrom,
                               read_start, read_end, strand, n_exons, read_length,
                               annot_gene_id, annot_transcript_id, 
                               annot_gene_name, annot_transcript_name, 
                               curr_gene_novelty, curr_transcript_novelty, 
                               curr_ISM_novelty]) + "\n")

    o.close()
    os.system("rm " + tmp_read_file)

def check_build_validity(build, database):
    """ Make sure that the user has entered a correct build name """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    cursor.execute("SELECT name FROM genome_build")
    builds = [str(x[0]) for x in cursor.fetchall()]
    conn.close()

    if build == None:
        message = "Please provide a valid genome build name. " + \
                  "In this database, your options are: " + \
                  ", ".join(builds)
        raise ValueError(message)

    if build not in builds:
        message = "Build name '" + build + \
                  "' not found in this database. Try one of the following: " + \
                  ", ".join(builds)
        raise ValueError(message)

    return

def main():
    options = get_args()
    database = options.database
    build = options.build
    datasets = options.datasets
    outprefix = options.outprefix

    check_build_validity(build, database)
    # Make sure that the input database exists!
    if not Path(database).exists():
        raise ValueError("Database file '%s' does not exist!" % database)

    if datasets != "all":
        datasets = datasets.split(",")   
   
    make_read_annot_file(database, build, outprefix, datasets = datasets)
    

if __name__ == '__main__':
    main()
