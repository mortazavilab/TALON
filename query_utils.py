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

def fetch_all_known_genes_detected(cursor, dataset):
    """ Get the IDs of all known genes found in a particular dataset (no 
        filter with respect to the type of transcript detected). """

    query = """SELECT DISTINCT(gene_ID) FROM observed
                   LEFT JOIN gene_annotations AS ga ON ga.ID = observed.gene_ID
                   WHERE (ga.attribute = 'gene_status' AND ga.value = 'KNOWN')
                   AND observed.dataset = ?;"""
    cursor.execute(query, [dataset])
    known_genes = [x[0] for x in cursor.fetchall()]
    return known_genes

def count_known_genes_detected(cursor, dataset):
    """ Count the number of known genes detected in the dataset (no filter
        with respect to the type of transcript detected). """

    known_genes = fetch_all_known_genes_detected(cursor, dataset)
    return len(known_genes)

def count_novel_genes_detected(cursor, dataset):
    """ Count the number of novel genes detected in the dataset (no filter
        with respect to the type of transcript detected). """

    novel_genes = fetch_all_novel_genes_detected(cursor, dataset)
    return len(novel_genes)

def fetch_all_novel_genes_detected(cursor, dataset):
    """ Get the IDs of all novel genes found in a particular dataset (no
        filter with respect to the type of transcript detected). """

    query = """SELECT DISTINCT(gene_ID) FROM observed
                   LEFT JOIN gene_annotations AS ga ON ga.ID = observed.gene_ID
                   WHERE (ga.attribute = 'gene_status' AND ga.value = 'NOVEL')
                   AND observed.dataset = ?;"""
    cursor.execute(query, [dataset])
    novel_genes = [x[0] for x in cursor.fetchall()]
    return novel_genes

def fetch_all_known_transcripts_detected(cursor, dataset):
    """ Get the IDs of all transcripts annotated as known. Does not include 
        novel FSMs """

    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'transcript_status' AND ta.value = 'KNOWN')
                   AND observed.dataset = ?;"""
    cursor.execute(query, [dataset])
    known_transcripts = [x[0] for x in cursor.fetchall()]
    return known_transcripts

def fetch_FSM_novel_transcripts(cursor, dataset):
    """ Fetch IDs of novel FSMs observed in the current dataset """

    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'FSM_transcript' AND ta.value = 'TRUE')
                   AND observed.dataset = ?;"""
    cursor.execute(query, [dataset])
    FSM_transcripts = [x[0] for x in cursor.fetchall()]
    return FSM_transcripts



    
