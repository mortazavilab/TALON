# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program reads in a GTF-formatted transcript annotation (ie GENCODE) and
# creates a SQLite database with gene, transcript, and exon tables.
# This database is used by the TALON pipeline to maintain a registry of 
# known annotations as well as novel discoveries.

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
    parser.add_option("--a", dest = "annot_name",
        help = "Name of supplied annotation (will be used to label data)",
        type = "string")
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

def add_gene_table(database):
    """ Add a table to the database to track genes. Attributes are:
        - Gene ID
        - Gene Name
        - Chromosome
        - Start
        - End
        - Strand
        - Length (start to end of gene)
        - Annotated (0/1)
        - Dataset of origin
        - Datasets in which gene was found (comma-delimited list)
    """

    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()

    # Add table and set primary key column, which will be the gene ID
    c.execute('CREATE TABLE "genes" ("identifier" TEXT PRIMARY KEY)')

    # Add more columns (no default values)
    c.execute('ALTER TABLE "genes" ADD COLUMN "name" TEXT')
    c.execute('ALTER TABLE "genes" ADD COLUMN "chromosome" TEXT')
    c.execute('ALTER TABLE "genes" ADD COLUMN "start" INTEGER')
    c.execute('ALTER TABLE "genes" ADD COLUMN "end" INTEGER')
    c.execute('ALTER TABLE "genes" ADD COLUMN "strand" TEXT')
    c.execute('ALTER TABLE "genes" ADD COLUMN "length" TEXT')
    c.execute('ALTER TABLE "genes" ADD COLUMN "annotated" INTEGER')
    c.execute('ALTER TABLE "genes" ADD COLUMN "dataset_of_origin" TEXT')
    c.execute('ALTER TABLE "genes" ADD COLUMN "datasets_containing" TEXT')

    conn.commit()
    conn.close()
    return

def add_transcript_table(database):
    """ Add a table to the database to track transcripts. Attributes are:
        - Transcript ID
        - Transcript Name
        - Gene ID
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
    c.execute('CREATE TABLE "transcripts" ("identifier" TEXT PRIMARY KEY)')

    # Add more columns (empty)
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "name" TEXT')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "gene_id" TEXT')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "chromosome" TEXT')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "start" INTEGER')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "end" INTEGER')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "strand" TEXT')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "length" INTEGER')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "exon_count" INTEGER')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "exon_ids" TEXT')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "annotated" INTEGER')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "dataset_of_origin" TEXT')
    c.execute('ALTER TABLE "transcripts" ADD COLUMN "datasets_containing" TEXT')

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
        - Gene ID
        - Transcript IDs (comma-delimited list)
        - Annotated (0/1)
        - Dataset of origin
        - Datasets in which exon was found (comma-delimited list)
    """

    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()

    # Add table and set primary key column, which will be the exon ID
    c.execute('CREATE TABLE "exons" ("identifier" TEXT PRIMARY KEY)')

    # Add more columns (empty)
    c.execute('ALTER TABLE "exons" ADD COLUMN "chromosome" TEXT')
    c.execute('ALTER TABLE "exons" ADD COLUMN "start" INTEGER')
    c.execute('ALTER TABLE "exons" ADD COLUMN "end" INTEGER')
    c.execute('ALTER TABLE "exons" ADD COLUMN "strand" TEXT')
    c.execute('ALTER TABLE "exons" ADD COLUMN "length" INTEGER')
    c.execute('ALTER TABLE "exons" ADD COLUMN "gene_id" TEXT')
    c.execute('ALTER TABLE "exons" ADD COLUMN "transcript_ids" TEXT')
    c.execute('ALTER TABLE "exons" ADD COLUMN "annotated" INTEGER')
    c.execute('ALTER TABLE "exons" ADD COLUMN "dataset_of_origin" TEXT')
    c.execute('ALTER TABLE "exons" ADD COLUMN "datasets_containing" TEXT')

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

    # Add table and set primary key column
    table_name = "counters"
    c.execute('CREATE TABLE "counters" ("category" TEXT PRIMARY KEY)')

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

def add_abundance_table(database):
    """ Adds a table that will store transcript counts for each tracked 
        dataset. Primary key is transcript ID, and other columns represent
        the datasets. """

    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()
 
    # Add table and set primary key column
    c.execute('CREATE TABLE "abundance" ("transcript_id" TEXT PRIMARY KEY)')

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

def populate_db(database, annot_name, genes, transcripts, exons):
    """ Iterate over GTF-derived gene, transcript, and exon entries in order
        to add a record for each in the database
    """
    # Connecting to the database file
    conn = sqlite3.connect(database)
    c = conn.cursor()
    for gene_id in genes:
        gene = genes[gene_id]
        add_gene_entry(c, annot_name, gene)
    for transcript_id in transcripts:
        transcript = transcripts[transcript_id]
        add_transcript_entry(c, annot_name, transcript)
    for exon_id in exons:
        exon = exons[exon_id]
        add_exon_entry(c, annot_name, exon)

    conn.commit()
    conn.close()
    
    return

def add_gene_entry(c, annot_name, gene):
    """ Given a conn.curser (c) to a database and a gene object,
        this function adds an entry to the database's gene table """

    gene_id = gene.identifier
    gene_name = gene.name
    chrom = gene.chromosome
    start = gene.start
    end = gene.end
    strand = gene.strand
    length = end - start + 1
    annot = 1
    
    cols = " (" + ", ".join([str_wrap_double(x) for x in ["identifier", "name",
           "chromosome", "start","end","strand","length","annotated",
           "dataset_of_origin"]]) + ") "
    vals = [gene_id, gene_name, chrom, start, end, strand, length, annot,
            annot_name]

    command = 'INSERT OR IGNORE INTO "genes"' + cols + "VALUES " + \
              '(?,?,?,?,?,?,?,?,?)'
    c.execute(command,vals)
    return

def add_transcript_entry(c, annot_name, transcript):
    """ Given a curser (c) to a database and a transcript object,
        this function adds an entry to the database's transcript table.
        It also updates the gene table as appropriate.
     """

    t_id = transcript.identifier
    t_name = transcript.name
    g_id = transcript.gene_id
    chrom = transcript.chromosome
    start = transcript.start
    end = transcript.end
    strand = transcript.strand
    length = transcript.get_length()
    n_exons = len(transcript.exons)
    exon_ids = ",".join([x.identifier for x in transcript.exons])
    annot = 1

    # Get gene name and add the transcript ID to the gene
    #gene_id = transcript.gene_id
    #get_gene = 'SELECT "name","transcript_ids" FROM "genes" ' + \
    #           'WHERE "identifier" = ?'
    #c.execute(get_gene,[g_id])
    #gene_name, gene_transcripts = c.fetchone()
    #gene_name = str(gene_name)
    #if gene_transcripts == None:
    #    gene_transcripts = t_id
    #else:
    #    if t_id not in gene_transcripts:
    #        gene_transcripts = str(gene_transcripts) + "," + t_id
    #update_gene = 'UPDATE "genes" SET "transcript_ids" = ' + \
    #          str_wrap(gene_transcripts) + ' WHERE "identifier" = ?'
    #c.execute(update_gene,[g_id])

    cols = " (" + ", ".join([str_wrap_double(x) for x in ["identifier", "name", 
                     "gene_id", "chromosome", "start","end",
                     "strand","length","exon_count","exon_ids", 
                     "annotated", "dataset_of_origin"]]) + ") "
    vals = [t_id, t_name, g_id, chrom, start, end, strand, 
            length, n_exons, exon_ids, annot, annot_name]
    
    command = 'INSERT OR IGNORE INTO "transcripts"' + cols + "VALUES " + \
              '(?,?,?,?,?,?,?,?,?,?,?,?)'
    c.execute(command,vals)

    # Add transcript_id to abundance table
    c.execute('INSERT OR IGNORE INTO "abundance" ("transcript_id") VALUES (?)',
               [t_id])
    return

def add_exon_entry(c, annot_name, exon):
    """ Given a curser (c) to a database and an exon object,
        this function adds an entry to the database's exon table.
    """
    # Information we can get from the transcript object
    exon_id = exon.identifier
    chrom = exon.chromosome
    start = exon.start
    end = exon.end
    strand = exon.strand
    length = exon.length
    gene_id = exon.gene_id
    transcript_ids = ",".join(list(exon.transcript_ids))
    annot = 1    

    cols = " (" + ", ".join([str_wrap_double(x) for x in ["identifier", 
                    "chromosome", "start","end","strand","length","gene_id",
                    "transcript_ids", "annotated", "dataset_of_origin"]]) + ") "
    vals = [exon_id, chrom, start, end, strand, length, gene_id, transcript_ids, 
            annot, annot_name]

    command = 'INSERT OR IGNORE INTO "exons"' + cols + "VALUES " + \
              '(?,?,?,?,?,?,?,?,?,?)'
    c.execute(command,vals)
    return

def str_wrap_double(string):
    """ Adds double quotes around the input string """
    return '"' + string + '"'

def main():
    options = getOptions()
    gtf_file = options.gtf
    outprefix = options.outprefix
    annot_name = options.annot_name

    # Make database
    db_name = outprefix + ".db"
    create_database(db_name)
    add_transcript_table(db_name)
    add_gene_table(db_name)
    add_exon_table(db_name)
    add_counter_table(db_name)
    add_abundance_table(db_name)

    # Process the GTF annotations
    genes, transcripts, exons = read_gtf_file(gtf_file)
    populate_db(db_name, annot_name, genes, transcripts, exons)

if __name__ == '__main__':
    main()
