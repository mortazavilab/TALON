# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

from genetree import GeneTree
from optparse import OptionParser

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
            genes: A dictionary mapping each chromosome to an interval tree
            data structure. Each interval tree contains intervals corresponding 
            to gene class objects. 
    """
    genes = GeneTree()

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

            # Process genes
            if entry_type == "gene":
                add_gene_from_gtf(tab_fields, genes)

    genes.print_tree()
                
def add_gene_from_gtf(gene, gene_tree):
    """ Adds gene from a GTF file to an existing dictionary of chromosome 
        interval trees.

        Args: 
            gene: A list containing fields from a GTF file gene entry. 
            Example:
            ['chr1', 'HAVANA', 'gene', '11869', '14409', '.', '+', '.', 
             'gene_id "ENSG00000223972.5"; 
             gene_type "transcribed_unprocessed_pseudogene"; 
             gene_status "KNOWN"; gene_name "DDX11L1"; level 2; 
             havana_gene "OTTHUMG00000000961.2";'] 

            gene_tree: Object that stores genes as intervals. The provided gene
            will be added to this structure.
 
        Returns:
            gene_tree: Input gene_tree with one additional gene added to it.
    """
    chromosome = gene[0]
    description = gene[-1]
    gene_name = (description.split("gene_name ")[1]).split('"')[1]
    start = int(gene[3])
    end = int(gene[4])
    strand = gene[6]

    gene_tree.add_gene(gene_name, chromosome, start, end, strand)
    return gene_tree 

def main():
    options = getOptions()
    infile_list = options.infile_list
    gtf_file = options.gtf_file

    # Process the GTF annotations
    read_gtf_file(gtf_file)

if __name__ == '__main__':
    main()
