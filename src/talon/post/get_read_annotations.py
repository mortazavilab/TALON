# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# get_read_annotations.py is a utility that queries a TALON 
# database in order to get read-specific annotation information.

import sqlite3
from .. import query_utils as qutils
from string import Template

def fetch_reads(database, build, tmp_file = None, datasets = "all"):
    """ Performs database query to fetch location and gene/transcript assignment
        info for each long read in the specified datasets. 
        If tmp_file is set to None (default), then the function will return
        the query results in a list of lists. If an alternate value is provided, 
        then the results will be written to a tmp file of that name, and the 
        function will return a list of all the transcript IDs that were 
        observed."""

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
            transcripts = set()
        else:
            reads = []

        colnames = []

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
            print(out_read)
           
            if tmp_file != None:
                o.write("\t".join([ str(x) for x in out_read ]) + "\n")
                transcripts.add(entry["transcript_ID"])
            else:
                reads.append(out_read)

    # Return results or close file
    if tmp_file != None:
        o.close()
        return transcripts
    else:
        return reads

def make_read_annot_file(database, build, annot_name, outprefix, datasets = "all"):
    """ Creates an output file with the following columns:
            1. 

        The annot_name input indicates which source should be used for annotation
        attributes such as gene names etc. By default, reads from all datasets
        in the database are included, but this can be modified by supplying a
        list/tuple of dataset names to the datasets parameter.
    """
    tmp_read_file = outprefix + "_reads.tmp"
    transcript_IDs = fetch_reads(database, build, tmp_file = tmp_read_file, 
                                 datasets = datasets)


def main():
    pass

if __name__ == '__main__':
    main()
