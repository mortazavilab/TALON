import sqlite3
from talon import talon
from talon import initialize_talon_database as itd
import os

def make_minimal_db_for_filtering(db_file, reads, datasets, annotations):
    """ Create a TALON-style database with 3 tables:
        dataset, observed, and transcript_annotations. Populate it with 
        mock reads (observed) in order to test the filtering function. 
    """
    # Remove if db exists
    if os.path.isfile(db_file):
        os.system("rm %s" %(db_file))

    itd.create_database(db_file)
    itd.add_dataset_table(db_file)
    itd.add_observed_table(db_file)
    itd.add_annotation_table(db_file, "transcript_annotations", "transcripts",
                             "transcript_ID")

    conn = sqlite3.connect(db_file)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Add reads to observed table
    cols = " (" + ", ".join([talon.str_wrap_double(x) for x in
                       ["obs_ID", "gene_ID", "transcript_ID", "read_name",
                        "dataset", "start_vertex", "end_vertex",
                        "start_exon", "end_exon", "start_delta", "end_delta", 
                        "read_length", "fraction_As", "custom_label", 
                        "allelic_label"]]) + ") "
    command = 'INSERT INTO "observed"' + cols + \
                          "VALUES " + '(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    cursor.executemany(command, reads)

    # Add transcript novelty annotations
    cols = " (" + ", ".join([talon.str_wrap_double(x) for x in
                       ["ID", "annot_name", "source", "attribute", "value"]]) + ") "
    command = 'INSERT OR IGNORE INTO "' + 'transcript_annotations" ' + \
              cols + "VALUES " + '(?,?,?,?,?)'
    cursor.executemany(command, annotations)    

    # Add datasets
    cols = " (" + ", ".join([talon.str_wrap_double(x) for x in
               ["dataset_ID", "dataset_name", "sample", "platform"]]) + ") "
    command = 'INSERT INTO "dataset"' + cols + \
                  "VALUES " + '(?,?,?,?)'
    cursor.executemany(command, datasets)


    conn.commit() 
    conn.close()
