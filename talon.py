# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

from gene import *
from genetree import *
from optparse import OptionParser
from sam_transcript import *
from transcript import *
import warnings

def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "infile_list", 
        help = "Comma-delimited list of input SAM files",
        metavar = "FILE", type = "string")
    parser.add_option("--gtf", "-g", dest = "gtf_file",
        help = "GTF annotation containing genes, transcripts, and exons.",
        metavar = "FILE", type = "string")
    (options, args) = parser.parse_args()
    return options

def read_gtf_file(gtf_file):
    """ Reads gene, transcript, and exon information from a GTF file.

        Args:
            gtf_file: Path to the GTF file

        Returns:
            genes: A GeneTree object, which consists of a dictionary mapping 
            each chromosome to an interval tree data structure. Each interval 
            tree contains intervals corresponding to gene class objects. 
    """
    genes = {} # GeneTree()
    currGene = None
    currTranscript = None

    with open(gtf_file) as gtf:
        for line in gtf:
            line = line.strip()
            
            # Ignore header
            if line.startswith("#"):
                continue

            # Split into constitutive fields on tab
            tab_fields = line.split("\t")
            chrom = tab_fields[0]
            entry_type = tab_fields[2]            

            if entry_type == "gene": 
                gene = get_gene_from_gtf(tab_fields)
                genes[gene.identifier] = gene
            elif entry_type == "transcript":
                transcript = get_transcript_from_gtf(tab_fields)
                gene_id = transcript.gene_id
                if gene_id not in genes:
                    warnings.warn("Tried to add transcript " + \
                    transcript.identifier + " to a gene that doesn't " + \
                    "exist in dict (" + gene_id + "). " + \
                    "Check order of GTF file.", RuntimeWarning)
                else:
                    genes[gene_id].add_transcript(transcript)
            elif entry_type == "exon": 
                info = tab_fields[-1]
                transcript_id = (info.split("transcript_id ")[1]).split('"')[1]
                gene_id = (info.split("gene_id ")[1]).split('"')[1]

                if gene_id not in genes:
                    warnings.warn("Tried to add exon to a gene that doesn't"+ \
                    " exist in dict (" + gene_id + "). " + \
                    "Check order of GTF file.", RuntimeWarning)
                elif transcript_id not in genes[gene_id].transcripts:
                    warnings.warn("Tried to add exon to a transcript (" + \
                    transcript_id + ") that isn't in "+ \
                    " gene transcript set (" + gene_id + "). " + \
                    "Check order of GTF file.", RuntimeWarning) 
                else:
                    currTranscript = genes[gene_id].transcripts[transcript_id]
                    currTranscript.add_exon_from_gtf(tab_fields)

    # Now create the GeneTree structure
    gt = GeneTree()
    for gene_id in genes:
        gene = genes[gene_id]
        gt.add_gene(gene)

    return gt

def process_sam_file(sam_file, genes):
    """ Reads transcripts from a SAM file

        Args:
            sam_file: Path to the SAM file

        Returns: 
    """

    transcripts = {}
    known_detected = 0
    transcripts_processed = 0

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
   
            transcripts_processed += 1
            sam_transcript = get_sam_transcript(sam)
            chromosome = sam_transcript.chromosome
            start = sam_transcript.start
            end = sam_transcript.end
            strand = sam_transcript.strand
            gene = get_best_gene_match(chromosome, start, end, strand, genes)
            if gene != None:
                transcript_match = gene.lookup_transcript_permissive_both(sam_transcript, False)
                if transcript_match != None:
                    print sam_transcript.sam_id
                    known_detected += 1
                
                #    exit()

            #exit()
    #print transcripts_processed
    #print known_detected

def main():
    options = getOptions()
    infile_list = options.infile_list
    gtf_file = options.gtf_file

    # Process the GTF annotations
    genes = read_gtf_file(gtf_file)

    # Process the SAM files
    sam_files = infile_list.split(",")
    for sam in sam_files:
        process_sam_file(sam, genes)

if __name__ == '__main__':
    main()
