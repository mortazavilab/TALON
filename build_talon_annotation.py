# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program reads in a GTF-formatted transcript annotation (ie GENCODE) and
# structures it in the annotation format used by the TALON pipeline.

import sqlite3
from sqlite3 import Error
from optparse import OptionParser
from gene import *
from transcript import *
from sam_transcript import *
from exon import *
import warnings
import os

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

def create_database(path):
    """ Creates an SQLite database to store transcript, gene, and exon 
        information """
    
    if os.path.isfile(path):
        raise ValueError("Database with name '" + path + "' already exists!")

    try:
        conn = sqlite3.connect(path)
    except Error as e:
        print(e)
    finally:
        conn.commit()
        conn.close()

    return

def add_transcript_table(database):
    """ Add a table to the database to track transcripts. Attributes are:
        - Transcript ID
        - Transcript Name
        - Gene ID
        - Gene Name
        - Chromosome
        - Start
        - End
        - Strand
        - Length (sum of exon lengths)
        - Exon count
        - Exon IDs
        - Exon coords (comma-delimited list)
        - Annotated (0/1)
        - Dataset of origin
        - Datasets in which transcript was found (comma-delimited list)
    """

    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()

    # Add table and set primary key column, which will be the transcript ID
    table_name = "transcripts"
    c.execute('CREATE TABLE {tn} ({nf} {ft} PRIMARY KEY)'\
        .format(tn=table_name, nf="identifier", ft="TEXT"))

    # Add more columns (empty)
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="name", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="gene_id", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="gene_name", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="chromosome", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="start", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="end", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="strand", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="length", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="exon_count", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="exon_ids", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="exon_coords", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="annotated", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="dataset_of_origin", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="datasets_containing_transcript", ct="TEXT"))

    conn.commit()
    conn.close()
    return

def add_gene_table(database):
    """ Add a table to the database to track genes. Attributes are:
        - Gene ID
        - Gene Name
        - Chromosome
        - Start
        - End
        - Strand
        - Transcript IDs (comma-delimited list)
        - Annotated (0/1)
        - Dataset of origin
        - Datasets in which gene was found (comma-delimited list)
    """

    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()

    # Add table and set primary key column, which will be the transcript ID
    table_name = "genes"
    c.execute('CREATE TABLE {tn} ({nf} {ft} PRIMARY KEY)'\
        .format(tn=table_name, nf="identifier", ft="TEXT"))  
    
    # Add more columns (empty)
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="name", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="chromosome", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="start", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="end", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="strand", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="transcript_ids", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="annotated", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="dataset_of_origin", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="datasets_containing_gene", ct="TEXT"))

    conn.commit()
    conn.close()
    return

def add_exon_table(database):
    """ Add a table to the database to track genes. Attributes are:
        - Exon ID
        - Chromosome
        - Start
        - End
        - Strand
        - Length
        - Transcript IDs (comma-delimited list)
        - Annotated (0/1)
        - Dataset of origin
        - Datasets in which exon was found (comma-delimited list)
    """

    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()

    # Add table and set primary key column, which will be the transcript ID
    table_name = "exons"
    c.execute('CREATE TABLE {tn} ({nf} {ft} PRIMARY KEY)'\
        .format(tn=table_name, nf="identifier", ft="TEXT"))

    # Add more columns (empty)
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="chromosome", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="start", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="end", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="strand", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="transcript_ids", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="length", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="annotated", ct="INTEGER"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="dataset_of_origin", ct="TEXT"))
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="datasets_containing_gene", ct="TEXT"))

    conn.commit()
    conn.close()
    return

def add_counter_table(database):
    """ Add a table to the database to track novel events. Attributes are:
        - Category (gene, transcript, exon)
        - Novel (number of novel events of that category so far)
    """

    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()

    # Add table and set primary key column, which will be the transcript ID
    table_name = "counters"
    c.execute('CREATE TABLE {tn} ({nf} {ft} PRIMARY KEY)'\
        .format(tn=table_name, nf="category", ft="TEXT"))

    # Add novel column
    default_val = 0
    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
        .format(tn=table_name, cn="novel", ct="INTEGER", df=default_val))

    # Add rows for genes, transcripts, and exons
    c.execute("INSERT INTO {tn} ({idf}, {cn}) VALUES ('genes', 0)".\
        format(tn=table_name, idf="category", cn="novel"))
    c.execute("INSERT INTO {tn} ({idf}, {cn}) VALUES ('transcripts', 0)".\
        format(tn=table_name, idf="category", cn="novel"))
    c.execute("INSERT INTO {tn} ({idf}, {cn}) VALUES ('exons', 0)".\
        format(tn=table_name, idf="category", cn="novel"))

    conn.commit()
    conn.close()
    return

def read_gtf_file(gtf_file):
    """ Reads gene, transcript, and exon information from a GTF file.

        Args:
            gtf_file: Path to the GTF file

        Returns:
            genes: A dictionary mapping gene IDs to corresponding gene objects
            transcripts: A dictionary mapping gene IDs to corresponding 
                   transcript objects
            exons: A dictionary mapping exon IDs to corresponding exon objects
    """
    genes = {}
    transcripts = {}
    exons = {} 

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

            # Entry is a gene
            if entry_type == "gene":
                gene = get_gene_from_gtf(tab_fields)
                genes[gene.identifier] = gene

            # Entry is a transcript
            elif entry_type == "transcript":
                transcript = get_transcript_from_gtf(tab_fields)
                gene_id = transcript.gene_id
                if gene_id not in genes:
                    warnings.warn("Tried to add transcript " + \
                    transcript.identifier + " to a gene that doesn't " + \
                    "exist in dict (" + gene_id + "). " + \
                    "Skipping this entry.", RuntimeWarning)
                else:
                    genes[gene_id].add_transcript(transcript)
                    transcripts[transcript.identifier] = transcript

            # Entry is an exon
            elif entry_type == "exon":
                info = tab_fields[-1]
                transcript_id = (info.split("transcript_id ")[1]).split('"')[1]
                gene_id = (info.split("gene_id ")[1]).split('"')[1]

                if gene_id not in genes:
                    warnings.warn("Tried to add exon to a gene that doesn't"+ \
                    " exist in dict (" + gene_id + "). " + \
                    "Skipping this entry.", RuntimeWarning)
                elif transcript_id not in genes[gene_id].transcripts:
                    warnings.warn("Tried to add exon to a transcript (" + \
                    transcript_id + ") that isn't in "+ \
                    " gene transcript set (" + gene_id + "). " + \
                    "Skipping this entry.", RuntimeWarning)
                else:
                    exon_id = (info.split("exon_id ")[1]).split('"')[1]
                    if exon_id not in exons:
                        # Create new exon entry
                        exon = create_exon_from_gtf(tab_fields)
                        exons[exon_id] = exon
                    else:
                        # Update existing exon entry
                        exon = exons[exon_id]
                        exon.transcript_ids.add(transcript_id)
                    currTranscript = genes[gene_id].transcripts[transcript_id]
                    currTranscript.add_exon(exon)

    return genes, transcripts, exons

def populate_db(database, genes, transcripts, exons):
    """ Iterate over GTF-derived gene, transcript, and exon entries in order
        to add a record for each in the database
    """
    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()
    for gene_id in genes:
        gene = genes[gene_id]
        add_gene(c, gene)
    for transcript_id in transcripts:
        transcript = transcripts[transcript_id]
        add_transcript(c, transcript)
    
    conn.commit()
    conn.close()
    
    return

def add_gene(c, gene):
    """ Given a conn.curser (c) to a database and a gene object,
        this function adds an entry to the database's gene table """
    gene_id = str_wrap(gene.identifier)
    gene_name = str_wrap(gene.name)
    chrom = str_wrap(gene.chromosome)
    start = str(gene.start)
    end = str(gene.end)
    strand = str_wrap(gene.strand)
    annot = "1"
    
    table_name = "genes"
    cols = ", ".join(['{idf}','{nm}','{chrom}','{st}','{ed}','{std}','{an}'])
    cols = " (" + cols + ") "
    vals = ", ".join([gene_id, gene_name, chrom, start, end, strand, annot])
    vals = " (" + vals + ") " 
    command = "INSERT OR IGNORE INTO {tn}" + cols + "VALUES" + vals   
    c.execute(command.format(tn=table_name, idf="identifier", nm="name", \
              chrom="chromosome", st="start", ed="end", std="strand", \
              an="annotated"))
    return

def add_transcript(c, transcript):
    """ Given a curser (c) to a database and a transcript object,
        this function adds an entry to the database's transcript table.
        It also updates the gene table as appropriate.
     """
    # Information we can get from the transcript object
    t_id = transcript.identifier
    t_name = str_wrap(transcript.name)
    g_id = transcript.gene_id
    chrom = str_wrap(transcript.chromosome)
    start = str(transcript.start)
    end = str(transcript.end)
    strand = str_wrap(transcript.strand)
    length = transcript.get_length()
    n_exons = len(transcript.exons)
    exon_ids = str_wrap(",".join([x.identifier for x in transcript.exons]))
    exon_coords = str_wrap(",".join(str(x) for x in 
                               transcript.get_exon_coords()))
    annot = "1"

    # Get gene name and add the transcript ID to the gene
    gene_id = transcript.gene_id
    get_gene = 'SELECT "name","transcript_ids" FROM "genes" ' + \
               'WHERE "identifier" = ?'
    c.execute(get_gene,[g_id])
    gene_name, gene_transcripts = c.fetchone()
    gene_name = str(gene_name)
    if gene_transcripts == None:
        gene_transcripts = t_id
    else:
        if t_id not in gene_transcripts:
            gene_transcripts = str(gene_transcripts) + "," + t_id
    update_gene = 'UPDATE "genes" SET "transcript_ids" = ' + \
              str_wrap(gene_transcripts) + ' WHERE "identifier" = ?'
    c.execute(update_gene,[g_id])

    cols = " (" + ", ".join([str_wrap_double(x) for x in ["identifier", "name", 
                     "gene_id", "gene_name", "chromosome", "start","end",
                     "strand","length","exon_count","exon_ids","exon_coords", 
                     "annotated"]]) + ") "
    vals = [t_id, t_name, gene_id, gene_name, chrom, start, end, strand, 
            length, n_exons, exon_ids, exon_coords, annot]
    
    command = 'INSERT OR IGNORE INTO "transcripts"' + cols + "VALUES " + \
              '(?,?,?,?,?,?,?,?,?,?,?,?,?)'
    c.execute(command,vals)
    return
    
def str_wrap(string):
    """ Adds single quotes around the input string """
    return "'" + string + "'"

def str_wrap_double(string):
    """ Adds double quotes around the input string """
    return '"' + string + '"'

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

    # Make database
    db_name = outprefix + ".db"
    #create_database(db_name)
    #add_transcript_table(db_name)
    #add_gene_table(db_name)
    #add_exon_table(db_name)
    #add_counter_table(db_name)

    # Process the GTF annotations
    genes, transcripts, exons = read_gtf_file(gtf_file)
    populate_db(db_name, genes, transcripts, exons)
    #make_tracker(genes, transcripts, exons, outprefix, annot_name)

if __name__ == '__main__':
    main()
