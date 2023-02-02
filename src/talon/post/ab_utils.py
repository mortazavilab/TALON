import sqlite3
import itertools
import operator
from optparse import OptionParser
from pathlib import Path
import scanpy
import numpy as np

from . import filter_talon_transcripts as filt
from .. import dstruct as dstruct
from .. import length_utils as lu
from . import post_utils as putils
from .. import query_utils as qutils
from .. import talon as talon

def check_annot_validity(annot, database):
    """ Make sure that the user has entered a correct annotation name """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    cursor.execute("SELECT DISTINCT annot_name FROM gene_annotations")
    annotations = [str(x[0]) for x in cursor.fetchall()]
    conn.close()

    if "TALON" in annotations:
        annotations.remove("TALON")

    if annot == None:
        message = "Please provide a valid annotation name. " + \
                  "In this database, your options are: " + \
                  ", ".join(annotations)
        raise ValueError(message)

    if annot not in annotations:
        message = "Annotation name '" + annot + \
                  "' not found in this database. Try one of the following: " + \
                  ", ".join(annotations)
        raise ValueError(message)

    return

def check_build_validity(build, database):
    """ Make sure that the user has entered a correct build name """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    cursor.execute("SELECT name FROM genome_build")
    builds = [str(x[0]) for x in cursor.fetchall()]
    conn.close()

    if build == None:
        message = "Please provide a valid genome build name. " + \
                  "In this database, your options are: " + \
                  ", ".join(builds)
        raise ValueError(message)

    if build not in builds:
        message = "Build name '" + build + \
                  "' not found in this database. Try one of the following: " + \
                  ", ".join(builds)
        raise ValueError(message)

    return

def fetch_naming_prefix(database):
    """ Get naming prefix from the database run_info table """
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    cursor.execute("SELECT value FROM run_info WHERE item = 'idprefix'")
    prefix = cursor.fetchone()[0]

    conn.close()
    return prefix

def fetch_n_places(database):
    """ Get length of name field from the database run_info table """
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    cursor.execute("SELECT value FROM run_info WHERE item = 'n_places'")
    n_places = cursor.fetchone()[0]

    conn.close()
    return int(n_places)

def get_transcript_lengths(database, build):
    """ Read the transcripts from the database. Then compute the lengths.
        Store in a dictionary """

    transcript_lengths = {}

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Get the exon lengths
    exon_lens = lu.get_all_exon_lengths(cursor, build)

    cursor.execute("SELECT * FROM transcripts")
    for transcript_row in cursor.fetchall():
        transcript_ID = transcript_row['transcript_ID']
        length = lu.get_transcript_length(transcript_row, exon_lens)
        transcript_lengths[transcript_ID] = length

    conn.close()
    return transcript_lengths

def fetch_dataset_list(dataset_file, database):
    """ Gets a list of all datasets in the database """

    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    all_db_datasets = qutils.fetch_all_datasets(cursor)
    conn.close()

    if dataset_file == None:

        return all_db_datasets

    else:
        datasets = []
        with open(dataset_file, 'r') as f:
            for line in f:
                dataset = line.strip()
                if dataset not in all_db_datasets:
                    raise ValueError("Dataset name '%s' not found in database" \
                                      % (dataset))
                datasets.append(dataset)

        return datasets
