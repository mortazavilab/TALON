# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

from exon import *
from exontree import *
from gene import *
from genetree import *
from optparse import OptionParser
from sam_transcript import *
from transcript import *
from transcript_match_tracker import *
import build_talon_annotation as ba
import sqlite3
import warnings

def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "infile_list", 
        help = "Comma-delimited list of input SAM files",
        metavar = "FILE", type = "string")
    parser.add_option("--d", dest = "dataset_label_list",
        help = "Comma-delimited list of dataset names (same order as sam)",
        type = "string")
    parser.add_option("--annot", "-a", dest = "annot",
        help = "TALON database. Created using build_talon_annotation.py",
        metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outfile",
        help = "Outfile name",
        metavar = "FILE", type = "string")
    parser.add_option("--bioRepMode", dest ="bioRepMode", action='store_true',
        help="In this mode, two or more biological replicate "+ \
        "sam files are run together. Novel transcripts discovered " + \
        "in at least two replicates are given a special label in " + \
        "the TALON database to allow for stricter filtering of "+ \
        "novel events.")
    
    (options, args) = parser.parse_args()
    return options

def read_annotation(annot):
    """ Imports data from the provided TALON database into gene, transcript, and
        exon objects. Also imports the number of novel discoveries from the 
        database so as to properly name discoveries in this run.
    """
    counter = {}
    gene_tree = GeneTree()
    transcripts = {}
    exon_tree = ExonTree()

    conn = sqlite3.connect(annot)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    # Get novel counters
    c.execute('SELECT * FROM counters')
    c.execute('SELECT "novel" FROM "counters" WHERE "category" = "genes"')
    counter["genes"] = int(c.fetchone()[0])
    c.execute('SELECT "novel" FROM "counters" WHERE "category" = "transcripts"')
    counter["transcripts"] = int(c.fetchone()[0])
    c.execute('SELECT "novel" FROM "counters" WHERE "category" = "exons"')
    counter["exons"] = int(c.fetchone()[0])

    # Require genes to be >= 200 bp in total length 
    c.execute('SELECT * FROM genes')
    for gene_row in c:
        if int(gene_row["length"]) >= 200:
            gene = get_gene_from_db(gene_row)
            gene_tree.add_gene(gene)

    c.execute('SELECT * FROM transcripts')
    exons_to_add = set()
    for transcript_row in c:
        # Only accept a transcript if we included the gene it belongs to
        if transcript_row['gene_id'] in gene_tree.genes:
            transcript = get_transcript_from_db(transcript_row)
            transcripts[transcript.identifier] = transcript

            # Add transcript to gene
            curr_gene_id = transcript_row['gene_id']
            curr_gene = gene_tree.genes[curr_gene_id]
            curr_gene.add_transcript(transcript)
        
            # Add exons that belong to transcript
            for exon_id in transcript_row['exon_ids'].split(","):
                exons_to_add.add(exon_id)

    # Add exons that belong to transcripts that we added    
    c.execute('SELECT * FROM exons')
    for exon_row in c:
        if exon_row['identifier'] in exons_to_add:
            exon = get_exon_from_db(exon_row)

            # Add exon to the transcripts it belongs to
            updated_ids = set()
            for transcript_id in exon_row['transcript_ids'].split(","):
                if transcript_id in transcripts:
                    curr_transcript = transcripts[transcript_id]
                    curr_transcript.add_exon(exon)
                    updated_ids.add(transcript_id)
           
            # Keep the exon if at least one transcript we're using has it
            if len(updated_ids) > 0:
                exon.transcript_ids = updated_ids
                exon_tree.add_exon(exon)
           
    conn.close()
    
    return gene_tree, transcripts, exon_tree, counter

def process_sam_file(sam_file):
    """ Reads transcripts from a SAM file

        Args:
            sam_file: Path to the SAM file

        Returns:
            sam_transcripts: List of sam_transcript objects
    """

    sam_transcripts = []

    with open(sam_file) as sam:
        for line in sam:
            line = line.strip()

            # Ignore header
            if line.startswith("@"):
                continue

            sam = line.split("\t")

            # Only use uniquely mapped transcripts for now
            if sam[1] not in ["0", "16"]:
                continue

            # Only use reads that are >= 200 bp long
            if len(sam[9]) < 200:
                continue
            
            sam_transcript = get_sam_transcript(sam)
            sam_transcripts.append(sam_transcript)

    return sam_transcripts

def identify_sam_transcripts(sam_transcripts, gene_tree, transcripts, 
                             exon_tree, counter, dataset, outprefix):
    """ Assign each sam transcript an annotated or a novel transcript identity


        Returns:
            Modified versions of gene_tree, transcripts, exon_tree, counter to
                which novel objects have been added
            abundance_dict: dictionary mapping transcript IDs to the number of 
                times each was observed in the sam data 
    """
    abundance_dict = {}
    novel_gene_ids = {}
    novel_transcript_ids = {}
    novel_exon_ids = {}

    out_sam = open(outprefix + "_talon.sam", 'w')
    out_txt = open(outprefix + "_talon.tsv", 'w')
    
    out_txt.write("\t".join(["dataset", "read_ID", "chromosome", "start", "end", 
                       "strand", "gene_id", "gene_name", "transcript_id", 
                       "transcript_name", "annotation_status", "match_type", 
                       "diff_5", "diff_3"]) + "\n")

    for sam_transcript in sam_transcripts:

        # Initialize
        gene_id = "NA"
        gene_name = "NA"
        transcript_id = "NA"
        transcript_name = "NA"
        annotation = "NA"
        chromosome = sam_transcript.chromosome
        start = sam_transcript.start
        end = sam_transcript.end
        strand = sam_transcript.strand

        # Look for full and partial matches
        best_match, match_type, exon_matches, diff = find_transcript_match(
                                                     sam_transcript, 
                                                     transcripts, exon_tree)
        # Full transcript match
        if match_type == "full":
            annot_transcript = best_match
            gene_id = annot_transcript.gene_id
            gene_name = gene_tree.genes[gene_id].name
            diff_5 = str(diff[0])
            diff_3 = str(diff[1])

        # If there is no full match, a novel transcript must be created
        else:
            # Find out which gene the novel transcript comes from. 
            if match_type == "partial":
                gene_id = best_match.gene_id
                gene_name = gene_tree.genes[gene_id].name

            # Search for gene using transcript coordinates                
            else:
                match_type = "none"
                chromosome = sam_transcript.chromosome
                start = sam_transcript.start
                end = sam_transcript.end
                strand = sam_transcript.strand
                #gene_candidates = gene_tree.get_genes_in_range(chromosome, start,
                #                                               end, strand)
                gene_obj, counter = create_novel_gene(chromosome, start, end,
                                                     strand, counter)
                gene_id = gene_obj.identifier
                novel_gene_ids[gene_id] = dataset
                gene_tree.add_gene(gene_obj)

            # Create the novel transcript and add to transcript dict
            # and novel tracker
            annot_transcript, counter = create_novel_transcript(chromosome, 
                                            start, end, strand, gene_id, counter) 

            novel_transcript_ids[annot_transcript.identifier] = dataset

            # Add exons to the novel transcript
            for sam_exon, exon_match in zip(sam_transcript.exons, exon_matches):
                exon_start = sam_exon.start
                exon_end = sam_exon.end               

                if exon_match == None:
                    exon_obj, counter = create_novel_exon(chromosome, exon_start, 
                                    exon_end, strand, counter)
                    novel_exon_ids[exon_obj.identifier] = dataset

                else: 
                    exon_obj = exon_tree.exons[exon_match]
                    if exon_start != exon_obj.start or exon_end != exon_obj.end:
                        exon_obj, counter = create_novel_exon(chromosome, 
                                              exon_start, exon_end, strand, counter)
                        novel_exon_ids[exon_obj.identifier] = dataset
 
                # Add exon to transcript and vice versa. Update data structure
                annot_transcript.add_exon(exon_obj)
                exon_obj.transcript_ids.add(annot_transcript.identifier)
                exon_tree.add_exon(exon_obj)

            # Set null 3' and 5' differences, other attributes
            diff_5 = "NA"
            diff_3 = "NA"
            #annotation = "novel"

            # Add transcript to data structure
            transcripts[annot_transcript.identifier] = annot_transcript
             
        # Update abundance dict
        transcript_id = annot_transcript.identifier
  
        if transcript_id in abundance_dict:
            abundance_dict[transcript_id] += 1
        else:
            abundance_dict[transcript_id] = 1

        if "talon" in gene_id or "talon" in transcript_id:
            annotation = "novel"
        else:
            annotation = "known"
        # Output assignment to file
        gene_id = annot_transcript.gene_id
        transcript_id = annot_transcript.identifier
        transcript_name = annot_transcript.name
        out_txt.write("\t".join([dataset, sam_transcript.identifier, \
                     sam_transcript.chromosome, str(sam_transcript.start), \
                     str(sam_transcript.end), sam_transcript.strand,
                     gene_id, gene_name, transcript_id, transcript_name, \
                     annotation, match_type, diff_5, diff_3]) + "\n")         
        
    out_txt.close()
    out_sam.close()   
                                
    return gene_tree, transcripts, exon_tree, counter, abundance_dict, \
           novel_gene_ids, novel_transcript_ids, novel_exon_ids


def find_transcript_match(query_transcript, transcripts, exon_tree):
    """ Performs search for matches to the query transcript, one exon at a time.

        Args:
            query_transcript: Transcript object to be matched

            transcripts: Dictionary of transcript_id -> transcript object.
            This is the catalog of transcripts that we've seen before.

            exon_tree: ExonTree structure storing the exons that we have seen
            before.            

        Returns:
            transcript_match: ___ if the query is matched to a known transcript 
                in full or in part. None otherwise.
            diff: [5' diff, 3' diff] from full transcript match, [None, None]
                  for partial or no match
            match_type: "full", "partial", or "none"

    """

    # Find annotation matches where possible for all exons
    transcript_match = "none"
    match_type = "none"
    diff = [None, None]
    tracker = MatchTracker(query_transcript)   
    tracker.match_all_exons(exon_tree) 

    # Find transcript matches
    tracker.compute_match_sets(transcripts)

    # If there is more than one full transcript match, select and return the
    # best one
    if len(tracker.full_matches) > 0:
        transcript_match, diff = tracker.get_best_full_match(transcripts)
        match_type = "full"  
        exon_matches = transcript_match.exons     

    # If there are no full matches, look for partial matches
    else:
        if len(tracker.partial_matches) > 0:
            transcript_match = tracker.get_best_partial_match(transcripts)
            match_type = "partial"

        exon_matches = tracker.get_best_exon_matches()

    return transcript_match, match_type, exon_matches, diff
        

def update_database(database, genetree, transcripts, exontree, counter, novel_gene_ids, 
                    novel_transcript_ids, novel_exon_ids):
    """ Add novel genes, transcripts, and exons to the supplied database. Also
        update existing entries.
    """
    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()

    n_genes = str(len(novel_gene_ids))
    n_transcripts = str(len(novel_transcript_ids))
    n_exons = str(len(novel_exon_ids))

    print "Adding " + n_genes + " novel genes to database..."
    for gene_id in novel_gene_ids.keys():
        dset = novel_gene_ids[gene_id]
        ba.add_gene_entry(c, dset, genetree.genes[gene_id])

    print "Adding " + n_transcripts + " novel transcripts to database..."
    for transcript_id in novel_transcript_ids.keys():
        dset = novel_transcript_ids[transcript_id]
        ba.add_transcript_entry(c, dset, transcripts[transcript_id])
 
    print "Adding " + n_exons + " novel exons to database..."
    for exon_id in novel_exon_ids.keys():
        dset = novel_exon_ids[exon_id]
        ba.add_exon_entry(c, dset, exontree.exons[exon_id])

    # Update the database counter
    update_g = 'UPDATE "counters" SET "novel" = ? WHERE "category" = "genes"'
    c.execute(update_g,[counter['genes']])

    update_t = 'UPDATE "counters" SET "novel" = ? WHERE "category" = "transcripts"'
    c.execute(update_t,[counter['transcripts']])

    update_e = 'UPDATE "counters" SET "novel" = ? WHERE "category" = "exons"'
    c.execute(update_e,[counter['exons']])

   
    conn.commit()
    conn.close() 
    return

def update_abundance_table(database, dset, abundance):
    """ Adds a column to the provided database that contains the abundance
        counts for each transcript observed in the dataset. """
    pass

def checkArgs(options):
    """ Makes sure that the options specified by the user are compatible with 
        each other """
    infile_list = options.infile_list
    dataset_list = options.dataset_label_list
    annot = options.annot
    bioRepMode = options.bioRepMode
    out = options.outfile

    # Number of data names must match number of datasets
    if len(dataset_list.split(",")) != len(infile_list.split(",")):
        raise ValueError('Number of provided data names must match number ' + \
                         'of datasets')

    # In Bio Rep Mode, more than one dataset is required
    if bioRepMode == True:
        if len(infile_list.split(",")) < 2:
            raise ValueError('At least two biological replicate datasets ' + \
                           'must be provided to run Biological Replicate Mode')  

    return

def main():
    options = getOptions()
    infile_list = options.infile_list
    dataset_list = options.dataset_label_list
    annot = options.annot
    bioRepMode = options.bioRepMode
    out = options.outfile

    # Check validity of input options
    checkArgs(options)

    # Process the annotations
    print "Processing annotation...................."
    gt, tscripts, et, ct = read_annotation(annot)

    # Process the SAM files
    print "Processing SAM file......................"
    sam_files = infile_list.split(",")
    dataset_list = dataset_list.split(",")
    for sam, d_name in zip(sam_files, dataset_list):
        print "Identifying transcripts in " + d_name + "..............."
        sam_tscripts = process_sam_file(sam)
        # TODO: novel dicts get overwritten when more than one dataset provided
        gt, tscripts, et, ct, abundance, novel_gene_ids, novel_transcript_ids, \
            novel_exon_ids= identify_sam_transcripts(sam_tscripts, 
                                          gt, tscripts, et, ct, d_name, out)
        update_abundance_table(d_name, abundance)

    update_database(annot, gt, tscripts, et, ct, novel_gene_ids, novel_transcript_ids,
                    novel_exon_ids) 

if __name__ == '__main__':
    main()
