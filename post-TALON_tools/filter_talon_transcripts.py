# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# filter_talon_transcripts.py is a utility that filters the transcripts inside
# a TALON database to produce a transcript whitelist. This list can then be 
# used by downstream analysis tools to determine which transcripts and other
# features should be reported (for example in a GTF file).

from optparse import OptionParser
import sqlite3
import warnings

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
    for pairing in pairing_list:  
        pairing_string = "(" + ','.join(['"' + x + '"' for x in pairing]) + ")"
        print "Processing group: " + pairing_string + "..."
        
        if len(pairing) <= 1 and novel_filtered == True:
            print "Warning: Only one dataset in group. This means that no " + \
                   "novel transcripts will pass the filter for this group." 

        # Query that joins transcript table to gene/transcript status, and also 
        # counts the number of datasets each transcript was detected in
        query = """
                SELECT 
                   t.gene_ID,
                   t.transcript_ID,
                   t.path,
                   ga.value,
                   ta.value,
                   COUNT(*)
               FROM transcripts t
               LEFT JOIN gene_annotations ga ON t.gene_ID = ga.ID
               LEFT JOIN transcript_annotations ta ON t.transcript_ID = ta.ID 
               INNER JOIN abundance ON t.transcript_ID = abundance.transcript_ID
               WHERE (ga.annot_name = %s OR ga.annot_name = "talon_run")
                   AND ga.attribute = "gene_status" 
                   AND (ta.annot_name = %s  OR ta.annot_name = "talon_run") 
                   AND ta.attribute = "transcript_status"
                   AND abundance.dataset IN %s
               GROUP BY t.transcript_ID;
            """
        annot_str = '"' + annot + '"'
        cursor.execute(query % (annot_str, annot_str, pairing_string))
        transcripts = cursor.fetchall()

        # Iterate over transcripts
        for transcript in transcripts:
            gene_ID = str(transcript[0])
            transcript_ID = str(transcript[1])
            path = transcript[2] 
            gene_status = transcript[3] 
            transcript_status = transcript[4]
            n_datasets = transcript[5]
           
            # Decide whether to add transcript to whitelist or not 
            if transcript_status == "KNOWN" or novel_filtered == False:
                transcript_whitelist.add((gene_ID, transcript_ID))
            else:
                n_exons = (len(path.split(",")) + 1)/2
                if n_datasets >= 2:
                    if novel_multiexon_reqmt == False:
                        transcript_whitelist.add((gene_ID, transcript_ID))
                    elif n_exons > 1:
                        transcript_whitelist.add((gene_ID, transcript_ID))
                

    # Disconnect from database
    conn.close()

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

def main():
    options = getOptions()

    if options.pairings_file != None:
        pairings = process_pairings(options.pairings_file)
        whitelist = filter_talon_transcripts(options.database, options.annot,
                                             dataset_pairings = pairings)
    else:
        whitelist = filter_talon_transcripts(options.database, options.annot)

    # Write transcript IDs to file
    o = open(options.outfile, 'w')
    print "Writing whitelisted gene-transcript pairs to " + options.outfile + "..."
    for transcript in whitelist:
        o.write(transcript[0] + "," + transcript[1] + "\n")
    o.close()
    

if __name__ == '__main__':
    main()
