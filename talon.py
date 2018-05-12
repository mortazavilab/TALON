# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

from gene import *
from exontree import *
from genetree import *
from optparse import OptionParser
from sam_transcript import *
from transcript import *
from transcript_match_tracker import *
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
    
    (options, args) = parser.parse_args()
    return options

def read_annotation(annot):
    """ Imports data from the provided TALON database into gene, transcript, and
        exon objects. Also imports the number of novel discoveries from the 
        database so as to properly name discoveries in this run.

        Args:
            annot

        Returns: 

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
    c.execute('SELECT * FROM genes WHERE "length" >= 200')
    for gene_row in c:
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
        
            for exon_id in transcript_row['exon_ids']:
                exons_to_add.add(exon_id)
    
    # Add exons that belong to transcripts that we added    
    c.execute('SELECT * FROM exons')
    for exon_row in c:
        if exon_row['identifier'] in exons_to_add:
            exon = get_exon_from_db(exon_row)
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

    return sam_transcripts

def identify_sam_transcripts(sam_transcripts, genes, transcripts, exon_tree):

    for sam_transcript in sam_transcripts:
        #print sam_transcript.sam_id
        match, diff, transcripts, genes, exon_tree = look_for_transcript_matches(sam_transcript, transcripts, genes, exon_tree)

        if match != None:
            #print sam_transcript.sam_id + " " +  match.identifier + " " + match.gene_id
            gene_id = match.gene_id
            gene_name = genes[match.gene_id].name
            transcript_id = match.identifier
            transcript_name = match.name
            if "novel" in transcript_id:
                annotation = "novel"
            else:
                annotation = "known"
        else:
            gene_id = "NA"
            gene_name = "NA"
            transcript_id = "NA"
            transcript_name = "NA"
            annotation = "novel"
        if None not in diff:
            diff_5 = str(diff[0])
            diff_3 = str(diff[1])
        else:
            diff_5 = "NA"
            diff_3 = "NA"
        transcripts_processed += 1
        #o.write("\t".join([dataset, sam_transcript.sam_id, \
        #             sam_transcript.chromosome, str(sam_transcript.start), \
        #             str(sam_transcript.end), sam_transcript.strand, 
        #             gene_id, gene_name, transcript_id, transcript_name, \
        #             annotation, diff_5, diff_3]) + "\n") 
    return genes, transcripts, exon_tree

def look_for_transcript_matches(query_transcript, transcripts, genes, exon_tree):
    """ Performs search for matches to the query transcript, one exon at a time.

        Args:
            query_transcript: Transcript object to be matched

            transcripts: Dictionary of transcript_id -> transcript object.
            This is the catalog of transcripts that we've seen before.

            exon_tree: ExonTree structure storing the exons that we have seen
            before.            

        Returns:

    """
    exons = query_transcript.exons
    chromosome = query_transcript.chromosome
    strand = query_transcript.strand
    tracker = MatchTracker(len(exons)/2)   
 
    exon_num = 0
    for i in range(0, len(exons), 2):
        start = exons[i]
        end = exons[i+1]

        cutoff_5, cutoff_3 = set_cutoffs_permissibleEnds(i, \
                             len(exons), strand)
        matches, diffs = get_loose_exon_matches(chromosome, start, end, \
                                                strand, exon_tree, cutoff_5, \
                                                cutoff_3)

        # Record exon matches and get transcripts that contain the exon  matches
        t_matches = set()
        for match, diff in zip(matches, diffs):
            tracker.add_exon_match(exon_num, chromosome, start, end, strand, \
                                   match, diff[0], diff[1])
            t_matches = t_matches.union(match.transcript_ids)
        tracker.transcript_matches.append(t_matches)
        exon_num += 1

    tracker.compute_match_sets(transcripts)
    if len(tracker.full_matches) > 0:
        transcript_match, diff = tracker.get_best_full_match(query_transcript, transcripts)
    else:
        # Create a novel transcript. Before the object can be created, the 
        # transcript must be assigned to a gene, or a novel gene created.
        diff = [None, None]
        transcript_match, transcripts, genes, exon_tree = create_novel_transcript(query_transcript, tracker, transcripts, genes, exon_tree)    
    return transcript_match, diff, transcripts, genes, exon_tree
        
def create_novel_transcript(sam_transcript, tracker, transcripts, genes, exon_tree):
    """ Creates a novel transcript from a sam_transcript object that can be
        added to the annotation. To do this, it first assigns the transcript to
        a gene so that an appropriate ID can be assigned. 

        Args:
            sam_transcript: sam_transcript object to be added

            tracker: MatchTracker object pertaining to the sam_transcript. Used
            to check whether the transcript has a partial match to any gene

            transcripts: a dictionary of transcript_id -> transcript object

            genes: a dictionary of gene_id -> gene object
    """
    chromosome = sam_transcript.chromosome 
    t_start = sam_transcript.start
    t_end = sam_transcript.end
    strand = sam_transcript.strand
    novel_transcript = None

    # Assign to partial match gene if possible
    if len(tracker.partial_matches) > 0:
        part_match = transcripts[tracker.partial_matches[0]]
        gene_match = genes[part_match.gene_id]

        # Create new IDs
        gene_match.novel += 1
        gene_id = gene_match.identifier
        gene_name = gene_match.name
        new_id = gene_id + "_novel_transcript." + str(gene_match.novel)
        new_name = new_id

        # Create novel transcript:
        novel_transcript = Transcript(new_id, new_name, chromosome, \
                                      t_start, t_end, strand, gene_id)
        transcripts[new_id] = novel_transcript       
 
        # Deal with exons
        exons = sam_transcript.exons
        exon_num = 0
        for i in range(0, len(exons), 2):
            start = exons[i]
            end = exons[i+1]
            perfect_match = None
            curr_exon_matches = tracker.exon_matches[exon_num]
            for entry in curr_exon_matches:
                if entry.tot_diff == 0:
                    # The exon entry is an exact match to the sam exon
                    # Add the current transcript ID to the exon
                    exon_id = entry.obj_id
                    exon_obj = exon_tree.exons[exon_id]
                    perfect_match = exon_obj
                    exon_t_set = exon_obj.transcript_ids
                    exon_t_set.add(new_id)
                    break
            if perfect_match == None:
                # Create a novel exon and add the transcript to it
                exon_tree.add_novel_exon(chromosome, start, end, strand, gene_id, new_id)
            transcripts[new_id].add_exon(start, end)
        
            exon_num += 1
    
    return novel_transcript, transcripts, genes, exon_tree
        

def set_cutoffs_permissibleEnds(exon_index, n_exons, strand):
    """ Selects 3' and 5' difference cutoffs depending on the exon index and 
        strand. For internal exons, both cutoffs should be set to zero, but
        more flexibility is allowed on the 3' and 5' ends.

        Args:
            exon_index: index we are at in the exon vector
            tot_exons: total number of exons in the transcript
            strand: strand that the exon is on. 
    """
    cutoff_5 = 0
    cutoff_3 = 0

    if exon_index == 0:
        if strand == "+":
            cutoff_5 = 100000
        else:
            cutoff_3 = 100000
    if exon_index == n_exons - 2:
        if strand == "+":
            cutoff_3 = 100000
        else:
            cutoff_5 = 100000
    return cutoff_5, cutoff_3

def main():
    options = getOptions()
    infile_list = options.infile_list
    dataset_list = options.dataset_label_list
    annot = options.annot
    out = options.outfile

    # Process the annotations
    genes, transcripts, exons, counter = read_annotation(annot)

    # Process the SAM files
    sam_files = infile_list.split(",")
    dataset_list = dataset_list.split(",")
    for sam, d_name in zip(sam_files, dataset_list):
        sam_transcripts = process_sam_file(sam)


    #o = open(out, 'w')
    #o.write("\t".join(["dataset", "read_ID", "chromosome", "start", "end", \
    #                   "strand", "gene_id", "gene_name", "transcript_id", \
    #                   "transcript_name", "annotation_status", "diff_5", \
    #                   "diff_3"]) + "\n")
    #sam_files = infile_list.split(",")
    #dataset_list = dataset_list.split(",")
    #for sam, dat in zip(sam_files, dataset_list):
    #    genes, transcripts, exons = process_sam_file(sam, dat, genes, transcripts, exons, o)
    #o.close()

if __name__ == '__main__':
    main()
