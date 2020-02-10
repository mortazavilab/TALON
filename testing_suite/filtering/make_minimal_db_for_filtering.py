import sqlite3
from talon import talon
from talon import initialize_talon_database as itd
import os

def init_mock_db(db_file):
    """ Initialize the minimal database needed to run the next tests.
        For testing inclusion of Knowns
         1. A known transcript (id = 1) (dataset 1) with low fracA
         2. A known transcript (id = 2) (dataset 2) with high fracA
        For testing selective exclusion of genomic transcripts and min_dataset
        and min_count criteria
         3. A genomic transcript (id = 3) (dataset 1) with low fracA
         4. A genomic transcript (id = 3) (dataset 2) with low fracA
         5. A genomic transcript (id = 3) (dataset 3) with low fracA
        For testing inclusion of Novel transcripts with high FracA
         6. An ISM transcript (id = 4) (dataset 4) with high fracA
         7. An ISM transcript (id = 4) (dataset 4) with high fracA
         8. An ISM transcript (id = 4) (dataset 5) with high fracA
    """

    # Add reads. Fields that are not relevant for this purpose are set to None
    known = [(1, 1, 1, "read_1", "dataset_1", None, None, None, None, None, None,
              None, 0.2, None, None),
             (2, 1, 2, "read_2", "dataset_2", None, None, None, None, None, None,
              None, 0.7, None, None)]
    genomic = [(3, 1, 3, "read_3", "dataset_1", None, None, None, None, None, None,
                None, 0.2, None, None),
               (4, 1, 3, "read_4", "dataset_2", None, None, None, None, None, None,
               None, 0.2, None, None),
               (5, 1, 3, "read_5", "dataset_3", None, None, None, None, None, None,
              None, 0.2, None, None)]
    ISM = [(6, 1, 4, "read_6", "dataset_4", None, None, None, None, None, None,
                None, 0.7, None, None),
               (7, 1, 4, "read_7", "dataset_4", None, None, None, None, None, None,
               None, 0.8, None, None),
               (8, 1, 4, "read_8", "dataset_5", None, None, None, None, None, None,
              None, 0.9, None, None)]


    reads = known + genomic + ISM

    # Datasets
    datasets = [(1, "dataset_1", "test", "test"),
                (2, "dataset_2", "test", "test"),
                (3, "dataset_3", "test", "test"),
                (4, "dataset_4", "test", "test"),
                (5, "dataset_5", "test", "test")]

    # Annotations
    annotations = [(1, "toy", "", "transcript_status", "KNOWN"),
                   (2, "toy", "", "transcript_status", "KNOWN"),
                   (3, "TALON", "", "transcript_status", "NOVEL"),
                   (3, "TALON", "", "genomic_transcript", "TRUE"),
                   (4, "TALON", "", "ISM_transcript", "TRUE")]

    make_minimal_db_for_filtering(db_file, reads, datasets, annotations)

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
