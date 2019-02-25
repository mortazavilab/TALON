# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Queries for interacting with a TALON database

import sqlite3

def count_observed_reads(cursor, dataset):
    """ Count the number of observed reads for the provided dataset """

    reads = cursor.execute(""" SELECT COUNT(obs_ID)
                                   FROM observed WHERE dataset = ? """,
                               [dataset]).fetchone()[0]
    return reads

def count_known_genes_detected(cursor, dataset):
    """ Count the number of known genes detected in the dataset (no filter). """

    query = """SELECT COUNT(DISTINCT(ga.ID)) 
                 FROM abundance 
                 LEFT JOIN transcripts AS t 
                   ON abundance.transcript_ID = t.transcript_ID
                 LEFT JOIN gene_annotations as ga 
                   ON ga.ID = t.gene_ID
                AND dataset = ?
                AND (ga.attribute = 'gene_status' AND ga.value = 'KNOWN')"""
 
    cursor.execute(query, [dataset])
    known_genes = cursor.fetchone()[0]
    return known_genes


