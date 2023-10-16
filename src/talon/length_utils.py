# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Queries for working with exon and transcript lengths


def get_all_exon_lengths(cursor, build):
    """Compute all exon lengths and store in a dict"""

    exon_lengths = {}
    cursor.execute(
        """ SELECT edge_ID,
                          loc1.position AS pos1,
                          loc2.position AS pos2,
                          abs(loc1.position - loc2.position) + 1 AS diff
                       FROM edge 
                       LEFT JOIN location AS loc1 ON edge.v1 = loc1.location_ID
                       LEFT JOIN location AS loc2 ON edge.v2 = loc2.location_ID
                       WHERE edge_type = 'exon' 
                       AND loc1.genome_build = '%s'
                       AND loc2.genome_build = '%s' """
        % (build, build)
    )

    for exon in cursor.fetchall():
        exon_ID = exon["edge_ID"]
        length = exon["diff"]
        exon_lengths[exon_ID] = length

    return exon_lengths


def get_transcript_length(transcript_row, exon_lengths):
    """Compute the length of the supplied transcript model based on its
    exons. Expected input format consists of a transcript row from a
    TALON database."""

    length = 0
    start_exon = transcript_row["start_exon"]
    end_exon = transcript_row["end_exon"]
    n_exons = transcript_row["n_exons"]

    if n_exons == 1:
        return exon_lengths[start_exon]
    else:
        jn_path = transcript_row["jn_path"].split(",")
        all_exons = [start_exon] + [int(x) for x in jn_path[1::2]] + [end_exon]

        for exon in all_exons:
            length += exon_lengths[exon]

        return length
