# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# Contains functions that query the database to initialize various data
# structures for the TALON run.
# ---------------------------------------------------------------------
# make_temp_novel_gene_table
# make_temp_monoexonic_transcript_table
# make_location_dict
# make_edge_dict
# make_transcript_dict
# make_vertex_2_gene_dict
# make_gene_start_and_end_dict

from string import Template
import pandas as pd

def make_temp_novel_gene_table(cursor, build, chrom = None, start = None,
                               end = None, tmp_tab = "temp_gene"):
    """ Attaches a temporary database with a table that has the following fields:
            - gene_ID
            - chromosome
            - start
            - end
            - strand
        The purpose is to track novel genes from this run in order to match
        transcripts to them when other forms of gene assignment have failed.
    """
    if any(val == None for val in [chrom, start, end]):
        command = Template(""" CREATE TEMPORARY TABLE IF NOT EXISTS $tmp_tab AS
                                   SELECT gene_ID,
                                     chromosome,
                                     start,
                                     end,
                                     strand
                                    FROM (SELECT g.gene_ID,
                                              loc.chromosome,
                                              MIN(loc.position) as start,
                                              MAX(loc.position) as end,
                                              g.strand
                                        FROM genes as g
                                        LEFT JOIN vertex as v ON g.gene_ID = v.gene_ID
                                        LEFT JOIN location as loc ON loc.location_ID = v.vertex_ID
                                        WHERE loc.genome_build = '$build'
                                        GROUP BY g.gene_ID); """)
    else:
        command = Template(""" CREATE TEMPORARY TABLE IF NOT EXISTS $tmp_tab AS
                                   SELECT gene_ID,
                                     chromosome,
                                     start,
                                     end,
                                     strand
                                    FROM (SELECT g.gene_ID,
                                              loc.chromosome,
                                              MIN(loc.position) as start,
                                              MAX(loc.position) as end,
                                              g.strand
                                        FROM genes as g
                                        LEFT JOIN vertex as v ON g.gene_ID = v.gene_ID
                                        LEFT JOIN location as loc ON loc.location_ID = v.vertex_ID
                                        WHERE loc.genome_build = '$build'
                                        GROUP BY g.gene_ID)
                                    WHERE chromosome = '$chrom'
                                        AND ((start <= $start AND end >= $end)
                                          OR (start >= $start AND end <= $end)
                                          OR (start >= $start AND start <= $end)
                                          OR (end >= $start AND end <= $end)); """)

    command = command.substitute({'tmp_tab':tmp_tab, 'build':build, 'chrom':chrom,
                                  'start':start, 'end':end})
    cursor.execute(command)

    return tmp_tab

def make_temp_transcript_table(cursor, build, chrom = None,
                                          start = None, end = None,
                                          tmp_tab = "temp_transcript"):
    """ Attaches a temporary database with a table that has the following fields:
            - gene_ID
            - transcript_ID
            - chromosome
            - start (min position)
            - end (max position)
            - strand
        The purpose is to allow location-based matching tiebreaking
        transcripts. """

    if any(val == None for val in [chrom, start, end]):
        command = Template(""" CREATE TEMPORARY TABLE IF NOT EXISTS $tmp_tab AS
                                   SELECT t.gene_ID,
                                      t.transcript_ID,
                                      loc1.chromosome,
                                      genes.strand,
                                      MIN(loc1.position, loc2.position) as min_pos,
                                      MAX(loc1.position, loc2.position) as max_pos
                                   FROM transcripts as t
                                   LEFT JOIN location as loc1
                                       ON loc1.location_ID = t.start_vertex
                                   LEFT JOIN location as loc2
                                       ON loc2.location_ID = t.end_vertex
                                   LEFT JOIN genes
                                       ON genes.gene_ID = t.gene_ID
                                   WHERE loc1.genome_build = '$build'
                                       AND loc2.genome_build = '$build' """)
    else:
        command = Template(""" CREATE TEMPORARY TABLE IF NOT EXISTS $tmp_tab AS
                                   SELECT t.gene_ID,
                                      t.transcript_ID,
                                      loc1.chromosome,
                                      genes.strand,
                                      t.start_exon as exon_ID,
                                      MIN(loc1.position, loc2.position) as min_pos,
                                      MAX(loc1.position, loc2.position) as max_pos
                                   FROM transcripts as t
                                   LEFT JOIN location as loc1
                                       ON loc1.location_ID = t.start_vertex
                                   LEFT JOIN location as loc2
                                       ON loc2.location_ID = t.end_vertex
                                   LEFT JOIN genes
                                       ON genes.gene_ID = t.gene_ID
                                   WHERE loc1.genome_build = '$build'
                                   AND loc2.genome_build = '$build'
                                   AND loc1.chromosome = '$chrom'
                                   AND ((min_pos <= $start AND max_pos >= $end)
                                       OR (min_pos >= $start AND max_pos <= $end)
                                       OR (min_pos >= $start AND min_pos <= $end)
                                       OR (max_pos >= $start AND max_pos <= $end))""")

    command = command.substitute({'build':build, 'chrom':chrom,
                                  'start':start, 'end':end,
                                  'tmp_tab':tmp_tab})
    cursor.execute(command)

    return tmp_tab

def make_temp_monoexonic_transcript_table(cursor, build, chrom = None,
                                          start = None, end = None,
                                          tmp_tab = "temp_monoexon"):
    """ Attaches a temporary database with a table that has the following fields:
            - gene_ID
            - transcript_ID
            - chromosome
            - start (min position)
            - end (max position)
            - strand
        The purpose is to allow location-based matching for monoexonic query
        transcripts. """

    if any(val == None for val in [chrom, start, end]):
        command = Template(""" CREATE TEMPORARY TABLE IF NOT EXISTS $tmp_tab AS
                                   SELECT t.gene_ID,
                                      t.transcript_ID,
                                      loc1.chromosome,
                                      loc1.position as start,
                                      loc2.position as end,
                                      genes.strand,
                                      t.start_vertex,
                                      t.end_vertex,
                                      t.start_exon as exon_ID,
                                      MIN(loc1.position, loc2.position) as min_pos,
                                      MAX(loc1.position, loc2.position) as max_pos
                                   FROM transcripts as t
                                   LEFT JOIN location as loc1
                                       ON loc1.location_ID = t.start_vertex
                                   LEFT JOIN location as loc2
                                       ON loc2.location_ID = t.end_vertex
                                   LEFT JOIN genes
                                       ON genes.gene_ID = t.gene_ID
                                   WHERE n_exons = 1
                                       AND loc1.genome_build = '$build'
                                       AND loc2.genome_build = '$build' """)
    else:
        command = Template(""" CREATE TEMPORARY TABLE IF NOT EXISTS $tmp_tab AS
                                   SELECT t.gene_ID,
                                      t.transcript_ID,
                                      loc1.chromosome,
                                      loc1.position as start,
                                      loc2.position as end,
                                      genes.strand,
                                      t.start_vertex,
                                      t.end_vertex,
                                      t.start_exon as exon_ID,
                                      MIN(loc1.position, loc2.position) as min_pos,
                                      MAX(loc1.position, loc2.position) as max_pos
                                   FROM transcripts as t
                                   LEFT JOIN location as loc1
                                       ON loc1.location_ID = t.start_vertex
                                   LEFT JOIN location as loc2
                                       ON loc2.location_ID = t.end_vertex
                                   LEFT JOIN genes
                                       ON genes.gene_ID = t.gene_ID
                                   WHERE n_exons = 1
                                   AND loc1.genome_build = '$build'
                                   AND loc2.genome_build = '$build'
                                   AND loc1.chromosome = '$chrom'
                                   AND ((min_pos <= $start AND max_pos >= $end)
                                       OR (min_pos >= $start AND max_pos <= $end)
                                       OR (min_pos >= $start AND min_pos <= $end)
                                       OR (max_pos >= $start AND max_pos <= $end))""")

    command = command.substitute({'build':build, 'chrom':chrom,
                                  'start':start, 'end':end,
                                  'tmp_tab':tmp_tab})
    cursor.execute(command)

    return tmp_tab

def make_location_dict(genome_build, cursor, chrom = None, start = None, end = None):
    """ Format of dict:
        chromosome -> dict(position -> SQLite3 row from location table)

        old:
            Key: chromosome, pos
            Value: SQLite3 row from location table
    """
    location_dict = {}

    if any(val == None for val in [chrom, start,end]):
        query = Template("""SELECT * FROM location WHERE genome_build = '$build' """)
    else:
        query = Template("""SELECT * FROM location
                            WHERE genome_build = '$build'
                            AND chromosome = '$chrom'
                            AND position >= $start
                            AND position <= $end""")
    query = query.substitute({'build':genome_build, 'chrom':chrom,
                              'start':start, 'end':end})
    cursor.execute(query)
    for location in cursor.fetchall():
        chromosome = location["chromosome"]
        position = location["position"]
        try:
            location_dict[chromosome][position] = location
        except:
            location_dict[chromosome] = {position: location}

    return location_dict

def make_edge_dict(cursor, build = None, chrom = None, start = None, end = None):
    """ Format of dict:
            Key: vertex1_vertex2_type
            Value: SQLite3 row from edge table
    """
    edge_dict = {}
    if any(val == None for val in [chrom, start, end, build]):
        query = """SELECT * FROM edge"""
    else:
        query = Template("""SELECT e.*
                            FROM edge AS e
                            LEFT JOIN location as loc1 ON e.v1 = loc1.location_ID
                            LEFT JOIN location as loc2 ON e.v2 = loc2.location_ID
                            WHERE loc1.genome_build = '$build' AND loc2.genome_build = '$build'
                                 AND loc1.chromosome = "$chrom"
                                 AND (loc1.position >= $start AND loc1.position <= $end)
                                 AND (loc2.position >= $start AND loc2.position <= $end);
                         """)
        query = query.substitute({'build':build, 'chrom':chrom,
                                  'start':start, 'end':end})
    cursor.execute(query)
    for edge in cursor.fetchall():
        vertex_1 = edge["v1"]
        vertex_2 = edge["v2"]
        edge_type = edge["edge_type"]
        key = (vertex_1, vertex_2, edge_type)
        edge_dict[key] = edge

    return edge_dict

def make_transcript_dict(cursor, build, chrom = None, start = None, end = None):
    """ Format of dict:
            Key: tuple consisting of edges in transcript path
            Value: SQLite3 row from transcript table
    """
    transcript_dict = {}
    if any(val == None for val in [chrom, start, end]):
         query = Template("""SELECT t.*,
                                loc1.chromosome as chromosome,
                                loc1.position as start_pos,
                                loc2.position as end_pos
                            FROM transcripts AS t
                                LEFT JOIN location as loc1 ON t.start_vertex = loc1.location_ID
                                LEFT JOIN location as loc2 ON t.end_vertex = loc2.location_ID
                                WHERE loc1.genome_build = '$build' AND loc2.genome_build = '$build';
                          """)

    else:
        query = Template("""SELECT t.*,
                                loc1.chromosome as chrom,
                                loc1.position as start_pos,
                                loc2.position as end_pos,
                                MIN(loc1.position, loc2.position) as min_pos,
                                MAX(loc1.position, loc2.position) as max_pos
                            FROM transcripts AS t
                                LEFT JOIN location as loc1 ON t.start_vertex = loc1.location_ID
                                LEFT JOIN location as loc2 ON t.end_vertex = loc2.location_ID
                                WHERE loc1.genome_build = '$build' AND loc2.genome_build = '$build'
                                         AND chrom == '$chrom'
                                         AND ((min_pos <= $start AND max_pos >= $end)
                                           OR (min_pos >= $start AND max_pos <= $end)
                                           OR (min_pos >= $start AND min_pos <= $end)
                                           OR (max_pos >= $start AND max_pos <= $end))""")

    query = query.substitute({'build':build, 'chrom':chrom,
                                  'start':start, 'end':end})
    cursor.execute(query)
    for transcript in cursor.fetchall():
        transcript_path = transcript["jn_path"]
        if transcript_path != None:
            transcript_path = transcript_path.split(",") + \
                              [transcript["start_exon"], transcript["end_exon"]]
            transcript_path = frozenset([ int(x) for x in transcript_path])
        else:
            transcript_path = frozenset([transcript["start_exon"]])
        transcript_dict[transcript_path] = transcript

    return transcript_dict

def make_vertex_2_gene_dict(cursor, build = None, chrom = None, start = None, end = None):
    """ Create a dictionary that maps vertices to the genes that they belong to.
    """
    vertex_2_gene = {}
    if any(val == None for val in [chrom, start, end, build]):
        query = """SELECT vertex_ID,
                          vertex.gene_ID,
                          strand
                       FROM vertex
                       LEFT JOIN genes ON vertex.gene_ID = genes.gene_ID"""
    else:
        query = Template("""SELECT vertex_ID,
                                   vertex.gene_ID,
                                   strand
                            FROM vertex
                                LEFT JOIN genes ON vertex.gene_ID = genes.gene_ID
                                LEFT JOIN location AS loc ON vertex.vertex_ID = loc.location_ID
                                WHERE loc.genome_build = '$build'
                                     AND loc.chromosome = '$chrom'
                                     AND (loc.position >= $start AND loc.position <= $end)
                         """)
        query = query.substitute({'build':build, 'chrom':chrom,
                                  'start':start, 'end':end})

    cursor.execute(query)
    for vertex_line in cursor.fetchall():
        vertex = vertex_line["vertex_ID"]
        gene = vertex_line["gene_ID"]
        strand = vertex_line["strand"]

        if vertex in vertex_2_gene:
            vertex_2_gene[vertex].add((gene, strand))
        else:
            vertex_2_gene[vertex] = set()
            vertex_2_gene[vertex].add((gene, strand))

    return vertex_2_gene

def make_gene_start_or_end_dict(cursor, build, mode, chrom = None, start = None, end = None):
    """ Select the starts (or ends) of known genes in the database and store
        in a dict.
        Format of dict:
            Key: gene ID from database
            Value: dict mapping positions to start vertices (or end vertices) of
                   KNOWN transcripts from that gene
    """
    if mode not in ["start", "end"]:
        raise ValueError(("Incorrect mode supplied to 'make_gene_start_or_end_dict'."
                          " Expected 'start' or 'end'."))

    output_dict = {}
    if any(val == None for val in [chrom, start,end]):
        query = """SELECT gene_ID,
                          %s_vertex as vertex,
                          loc1.position as %s
                   FROM transcripts
                   LEFT JOIN transcript_annotations as ta
                       ON ta.ID = transcripts.transcript_ID
                   LEFT JOIN location as loc1
                       ON transcripts.%s_vertex = loc1.location_ID
                   WHERE ta.attribute = 'transcript_status'
                         AND ta.value = 'KNOWN'
                         AND loc1.genome_build = '%s' """
        cursor.execute(query % (mode, mode, mode, build))

    else:
        query = Template("""SELECT  gene_ID,
                                    ${mode}_vertex as vertex,
                                    loc1.chromosome as chrom,
                                    loc1.position as $mode

                            FROM transcripts
                            LEFT JOIN transcript_annotations as ta
                                ON ta.ID = transcripts.transcript_ID
                            LEFT JOIN location as loc1
                                ON transcripts.${mode}_vertex = loc1.location_ID
                            WHERE ta.attribute = 'transcript_status'
                                  AND ta.value = 'KNOWN'
                                  AND loc1.genome_build = '$build'
                                  AND chrom = '$chrom'
                                  AND ($mode >= $start AND $mode <= $end)""")
        query = query.substitute({'build':build, 'chrom':chrom,
                                  'start':start, 'end':end, 'mode':mode})
        cursor.execute(query)

    for entry in cursor.fetchall():
        gene_ID = entry['gene_ID']
        vertex = entry['vertex']
        pos = entry[mode]

        try:
            output_dict[gene_ID][pos] = vertex
        except:
            output_dict[gene_ID] = {}
            output_dict[gene_ID][pos] = vertex

    return output_dict
