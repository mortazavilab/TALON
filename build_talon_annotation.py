# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program reads in a GTF-formatted transcript annotation (ie GENCODE) and
# structures it in the annotation format used by the TALON pipeline.

from optparse import OptionParser
from gene import *
from exontree import *
from transcript import *
from sam_transcript import *
import warnings


def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "gtf",
        help = "GTF annotation containing genes, transcripts, and exons.",
        metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix",
        help = "Outprefix for the annotation files",
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
    genes = {}
    transcripts = {}
    exons = ExonTree()
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
                    transcripts[transcript.identifier] = transcript
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
                    exon = create_exon_from_gtf(tab_fields)
                    exons.add_exon(exon, exon.identifier)


    return genes, transcripts, exons

def make_tracker(genes, transcripts, exons, outprefix, annot_name):
    """ Iterates over GTF-derived transcripts and extracts the following 
        information in order to create an annotation tracker entry:
            - Transcript ID
            - Transcript Name
            - Gene ID
            - Gene Name
            - Chromosome
            - Start
            - End
            - Strand 
            - Exon IDs
            - Exon starts (delimited list)
            - Exon ends (delimited list)

            These fields will have predetermined values
            - Annotated (0/1)
            - Dataset of origin
            - Datasets in which transcript was found so far (delimited list)
            
            Args:
                transcripts: a dictionary of Transcript objects obtained from
                the provided GTF file     
    """ 
    #genes_seen = {}

    o_transcript = open(outprefix + "_talon_transcript.annot", 'w')
    #o_gene = open(outprefix + "_talon_gene.annot")

    o_transcript.write("# version=1\n") # Annotation version
    o_transcript.write("# novel_genes=0\n") # Gene naming counter
    o_transcript.write("# novel_transcripts=0\n") # Transcript naming counter
    o_transcript.write("# novel_exons=0\n") # Exon naming counter
    o_transcript.write("\t".join(["transcript_id", "transcript_name", 
                                  "gene_id", "gene_name", "chromosome", 
                                  "transcript_start", "transcript_end", 
                                  "strand", "exon_ids", "exon_coords", 
                                  "annotated", "dataset_of_origin", 
                                  "datasets_containing_transcript"]) + "\n")
    transcript_entries = []
    for transcript_id in transcripts:
        # Filter out transcripts that are smaller than 200 bp
        # As a rule, we're mostly concerned with filtering out single-exon
        # transcripts such as miRNAs here, so we can just look at the start
        # and end.
        transcript = transcripts[transcript_id]
        if transcript.end - transcript.start < 200:
            continue
        transcript_name = transcript.name
        gene_id = transcript.gene_id
        gene_name = genes[gene_id].name
        chromosome = transcript.chromosome
        transcript_start = str(transcript.start)
        transcript_end = str(transcript.end)
        strand = transcript.strand

        # Fetch exon IDs
        exon_ids = []
        for i in range(0, len(transcript.exons), 2):
            exon_start = transcript.exons[i]
            exon_end = transcript.exons[i+1]
            exon = get_loose_exon_matches(chromosome, exon_start, exon_end, \
                                                strand, exons, 0, 0)[0][0]
            exon_ids.append(exon.identifier)
             
        exon_ids = ",".join(exon_ids)
        exon_coords = ",".join([ str(x) for x in transcript.exons])
        annotated = "1"
        dataset_of_origin = annot_name
        datasets_containing_transcript = "."

        outstr = "\t".join([ transcript_id, transcript_name, gene_id, \
                           gene_name, chromosome, transcript_start, \
                           transcript_end, strand, exon_ids, exon_coords, \
                           annotated, dataset_of_origin, \
                           datasets_containing_transcript ]) 
        o_transcript.write(outstr + "\n")
    o_transcript.close()
    return

def main():
    options = getOptions()
    gtf_file = options.gtf
    outprefix = options.outprefix
    annot_name = (gtf_file.split("/")[-1]).split(".gtf")[0]

    # Process the GTF annotations
    genes, transcripts, exons = read_gtf_file(gtf_file)
    make_tracker(genes, transcripts, exons, outprefix, annot_name)

if __name__ == '__main__':
    main()
