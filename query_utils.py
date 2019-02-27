# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Queries for interacting with a TALON database

import sqlite3

def fetch_reproducible_intergenic(cursor, datasets):
    """ Return the gene and transcript ID of any intergenic transcripts that were
        found in at least two of the supplied datasets """

    datasets = format_for_IN(datasets)
    query = """SELECT gene_ID,
                      a.transcript_ID
               FROM abundance as a
               LEFT JOIN transcript_annotations as ta
                   ON ta.ID = a.transcript_ID
               LEFT JOIN transcripts
                   ON transcripts.transcript_ID = a.transcript_ID
               WHERE ta.attribute = 'intergenic_transcript'
               AND a.dataset IN """ + datasets + \
           """ GROUP BY a.transcript_ID
               HAVING count(*) > 1;"""

    cursor.execute(query)
    intergenic = [(x[0], x[1], "intergenic_transcript") for x in cursor.fetchall()]
    return intergenic

def fetch_reproducible_antisense(cursor, datasets):
    """ Return the gene and transcript ID of any antisense transcripts that were
        found in at least two of the supplied datasets """

    datasets = format_for_IN(datasets)
    query = """SELECT gene_ID,
                      a.transcript_ID
               FROM abundance as a
               LEFT JOIN transcript_annotations as ta
                   ON ta.ID = a.transcript_ID
               LEFT JOIN transcripts
                   ON transcripts.transcript_ID = a.transcript_ID
               WHERE ta.attribute = 'antisense_transcript'
               AND a.dataset IN """ + datasets + \
           """ GROUP BY a.transcript_ID
               HAVING count(*) > 1;"""

    cursor.execute(query)
    antisense = [(x[0], x[1], "antisense_transcript") for x in cursor.fetchall()]
    return antisense

def fetch_reproducible_NNCs(cursor, datasets):
    """ Return the gene and transcript ID of any NNC transcripts that were
        found in at least two of the supplied datasets """

    datasets = format_for_IN(datasets)
    query = """SELECT gene_ID, 
                      a.transcript_ID 
               FROM abundance as a
	       LEFT JOIN transcript_annotations as ta
	           ON ta.ID = a.transcript_ID
               LEFT JOIN transcripts
	           ON transcripts.transcript_ID = a.transcript_ID
	       WHERE ta.attribute = 'NNC_transcript'
	       AND a.dataset IN """ + datasets + \
           """ GROUP BY a.transcript_ID
               HAVING count(*) > 1;"""

    cursor.execute(query)
    NNC = [(x[0], x[1], "NNC_transcript") for x in cursor.fetchall()]
    return NNC

def fetch_known_transcripts_with_gene_label(cursor, datasets):
    """ Fetch known transcripts along with the gene they belong to """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT gene_ID,transcript_ID FROM observed
                   LEFT JOIN transcript_annotations AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'transcript_status' AND ta.value = 'KNOWN')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    known_transcripts = [(x[0], x[1], "FSM_transcript") for x in cursor.fetchall()]
    return known_transcripts

def fetch_NIC_transcripts_with_gene_label(cursor, datasets):
    """ Fetch NIC transcripts along with the gene they belong to """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT gene_ID,transcript_ID FROM observed
                   LEFT JOIN transcript_annotations AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'NIC_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    known_transcripts = [(x[0], x[1], "NIC_transcript") for x in cursor.fetchall()]
    return known_transcripts

def count_observed_reads(cursor, dataset):
    """ Count the number of observed reads for the provided dataset """

    reads = cursor.execute(""" SELECT COUNT(obs_ID)
                                   FROM observed WHERE dataset = ? """,
                               [dataset]).fetchone()[0]
    return reads

def fetch_all_known_genes_detected(cursor, datasets):
    """ Get the IDs of all known genes found in a particular dataset (no 
        filter with respect to the type of transcript detected). """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(gene_ID) FROM observed
                   LEFT JOIN gene_annotations AS ga ON ga.ID = observed.gene_ID
                   WHERE (ga.attribute = 'gene_status' AND ga.value = 'KNOWN')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
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

def fetch_all_novel_genes_detected(cursor, datasets):
    """ Get the IDs of all novel genes found in a particular dataset (no
        filter with respect to the type of transcript detected). """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(gene_ID) FROM observed
                   LEFT JOIN gene_annotations AS ga ON ga.ID = observed.gene_ID
                   WHERE (ga.attribute = 'gene_status' AND ga.value = 'NOVEL')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    novel_genes = [x[0] for x in cursor.fetchall()]
    return novel_genes

def fetch_all_known_transcripts_detected(cursor, datasets):
    """ Get the IDs of all transcripts annotated as known. Does not include 
        novel FSMs """
 
    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'transcript_status' AND ta.value = 'KNOWN')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
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

def fetch_novel_transcripts(cursor, dataset):
    """ Fetch IDs of novel transcripts observed in the current dataset """

    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'transcript_status' AND ta.value = 'NOVEL')
                   AND observed.dataset = ?;"""
    cursor.execute(query, [dataset])
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts

def fetch_antisense_genes(cursor, datasets):
    """ Fetch IDs of antisense genes observed in the dataset(s) """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(gene_ID) FROM observed
                   LEFT JOIN gene_annotations AS ga ON ga.ID = observed.gene_ID
                   WHERE (ga.attribute = 'antisense_gene')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    genes = [x[0] for x in cursor.fetchall()]
    return genes

def fetch_intergenic_novel_genes(cursor, datasets):
    """ Fetch IDs of novel genes denoted as intergenic """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(gene_ID) FROM observed
                   LEFT JOIN gene_annotations AS ga ON ga.ID = observed.gene_ID
                   WHERE (ga.attribute = 'intergenic_novel')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    genes = [x[0] for x in cursor.fetchall()]
    return genes

def fetch_all_ISM_transcripts(cursor, datasets):
    """ Fetch IDs of all ISM transcripts """
    
    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations 
                       AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'ISM_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts

def fetch_prefix_ISM_transcripts(cursor, datasets):
    """ Fetch IDs of all ISM prefix transcripts """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations
                       AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'ISM-prefix_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts

def fetch_suffix_ISM_transcripts(cursor, datasets):
    """ Fetch IDs of all ISM suffix transcripts """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations
                       AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'ISM-suffix_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts

def fetch_NIC_transcripts(cursor, datasets):
    """ Fetch IDs of all NIC transcripts """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations
                       AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'NIC_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts

def fetch_NNC_transcripts(cursor, datasets):
    """ Fetch IDs of all NNC transcripts """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations
                       AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'NNC_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts

def fetch_antisense_transcripts(cursor, datasets):
    """ Fetch IDs of all antisense transcripts """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations
                       AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'antisense_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts

def fetch_genomic_transcripts(cursor, datasets):
    """ Fetch IDs of all genomic transcripts """

    datasets = format_for_IN(datasets)
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations
                       AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'genomic_transcript')
                   AND observed.dataset IN """ + datasets
    cursor.execute(query)
    transcripts = [x[0] for x in cursor.fetchall()]
    return transcripts


#-------------------------------------------------------------------------------
def format_for_IN(l):
    """ Converts input to string that can be used for IN database query """
    
    if type(l) is not list:
        l = [str(l)]

    return "(" + ','.join(['"' + x + '"' for x in l]) + ")" 
