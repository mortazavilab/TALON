# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program takes transcripts from one or more samples (SAM format) and
# assigns them transcript and gene identifiers based on a GTF annotation.
# Novel transcripts are assigned new identifiers.

import edge as Edge
import edgetree as EdgeTree
import gene as Gene
import genetree as GeneTree
from intervaltree import *
from optparse import OptionParser
import pdb
import sam_transcript as SamTranscript
import sqlite3
import transcript as Transcript
import transcript_match_tracker as TMT
import vertex as Vertex
import warnings

def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "config_file", 
        help = "Dataset config file: dataset name, sample description, " + \
               "platform, sam file (comma-delimited)", type = "string")
    parser.add_option("--annot", "-a", dest = "annot",
        help = "TALON database. Created using build_talon_annotation.py",
        metavar = "FILE", type = "string")
    parser.add_option("--build", "-b", dest = "build",
        help = "Genome build to use. Note: must be in the TALON database.",
        type = "string")
    parser.add_option("--o", dest = "outfile",
        help = "Outfile name",
        metavar = "FILE", type = "string")
    parser.add_option("--encode", dest ="encode_mode", action='store_true',
                      help = "If this option is set, TALON will require novel \
                      transcripts to be corroborated by at least one other dataset \
                      in order to be included in the output file.")
    parser.add_option("--noUpdate", dest ="noUpdate", action='store_true',
                      help = "If this option is set, the database will not be updated.\
                             Typically this is not used.")

    
    (options, args) = parser.parse_args()
    return options

def read_annotation(annot, genome_build):
    """ Imports data from the provided TALON database into gene, transcript, and
        exon objects. Also imports the number of novel discoveries from the 
        database so as to properly name discoveries in this run.
    """

    conn = sqlite3.connect(annot)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Get counters
    counter = get_counters(cursor)
   
    # Read in the genes and add them to the gene tree data structure
    gene_tree = read_genes(cursor, genome_build)

    vertices = read_vertices(cursor, genome_build)
 
    # Read in all exons and introns; add them to inteval tree data structures
    exon_tree = read_edges(cursor, genome_build, "exon")
    intron_tree = read_edges(cursor, genome_build, "intron")
     
    # Read in transcripts and edges
    transcripts = read_transcripts(cursor, exon_tree, intron_tree)

    conn.close()
    
    return gene_tree, transcripts, exon_tree, intron_tree, vertices, counter

def get_counters(cursor):
    """ Fetches counters for the different database categories and returns
        them in a dictionary """

    counter = {}
    cursor.execute('SELECT "count" FROM "counters" WHERE "category" = "genes"')
    counter["genes"] = int(cursor.fetchone()[0])
    cursor.execute('SELECT "count" FROM "counters" WHERE "category" = "transcripts"')
    counter["transcripts"] = int(cursor.fetchone()[0])
    cursor.execute('SELECT "count" FROM "counters" WHERE "category" = "vertex"')
    counter["vertices"] = int(cursor.fetchone()[0])
    cursor.execute('SELECT "count" FROM "counters" WHERE "category" = "edge"')
    counter["edges"] = int(cursor.fetchone()[0])
    cursor.execute('SELECT "count" FROM "counters" WHERE "category" = "dataset"')
    counter["datasets"] = int(cursor.fetchone()[0])
    cursor.execute('SELECT "count" FROM "counters" WHERE "category" = "observed_transcript_start"')
    counter["observed_starts"] = int(cursor.fetchone()[0])
    cursor.execute('SELECT "count" FROM "counters" WHERE "category" = "observed_transcript_end"')
    counter["observed_ends"] = int(cursor.fetchone()[0])

    return counter

def read_genes(cursor, build): 
    """ Fetches genes from database, creates Gene objects for each, and adds
        them to a GeneTree data structure"""
    gene_tree = GeneTree.GeneTree()
    query_starts = """SELECT gene_ID, chromosome,MIN(position),strand FROM 
                   (SELECT * from genes 
                        LEFT JOIN vertex ON genes.gene_ID = vertex.gene_ID 
                   ) AS SUBQUERY_START_PLUS
                    LEFT JOIN location ON SUBQUERY_START_PLUS.vertex_ID = location.location_ID 
                        WHERE genome_build = """ + str_wrap_double(build) + \
                        """  GROUP BY gene_ID; """ 
    query_ends = """SELECT gene_ID, chromosome,MAX(position),strand FROM
                   (SELECT * from genes
                        LEFT JOIN vertex ON genes.gene_ID = vertex.gene_ID
                   ) AS SUBQUERY_START_PLUS
                    LEFT JOIN location ON SUBQUERY_START_PLUS.vertex_ID = location.location_ID
                        WHERE genome_build = """ + str_wrap_double(build) + \
                        """  GROUP BY gene_ID; """

    cursor.execute(query_starts)
    starts = cursor.fetchall()
    cursor.execute(query_ends)
    ends = cursor.fetchall()

    for start, end in zip(starts, ends):
        gene = Gene.get_gene_from_db(start, end)
        gene_tree.add_gene(gene)
    return gene_tree

def read_vertices(cursor, build):
    """ Extracts vertices and their positions in the provided genome build from 
        the database and organizes them in the following data structure:
    
        chromosome -> location -> list of vertices (tuples)

        Note: Vertex and gene IDs are integers here
    """
    
    vertex_dict = {}
    query = """SELECT vertex_ID, chromosome, position, strand, gene_ID FROM
            vertex LEFT JOIN location ON vertex.vertex_id = location.location_ID
                  WHERE genome_build = """ + str_wrap_double(build) 
    cursor.execute(query)
    db_vertices = cursor.fetchall()
 
    for v in db_vertices:
        curr_chrom = v['chromosome']
        # Add new chromosome if necessary
        if curr_chrom not in vertex_dict:
            vertex_dict[curr_chrom] = {}

        # Add location
        pos = int(v['position'])
        new_vertex = Vertex.Vertex(v["vertex_ID"], curr_chrom, pos, 
                                   v["strand"], v["gene_ID"])
        if pos in vertex_dict[curr_chrom]:
            vertex_dict[curr_chrom][pos].append(new_vertex)
        else:
            vertex_dict[curr_chrom][pos] = [new_vertex]
    return vertex_dict

def read_edges(cursor, build, edge_type):
    """ Fetches edges of the specified type and organizes them in an
        interval tree. 'v1' refers to the first vertex in the edge, whereas 
        'v2' refers to the second"""

    edge_tree = EdgeTree.EdgeTree()

    # Get all Vertex1 locations
    v1_query = """SELECT edge_ID, vertex_ID, chromosome, position, strand, gene_id FROM
                  (SELECT edge_ID,v1,chromosome,position,strand FROM edge
                   LEFT JOIN location ON edge.v1 = location.location_ID
                  WHERE genome_build = """ + str_wrap_double(build) + \
                  """ AND edge_type = """ + str_wrap_double(edge_type) + """) AS SUBQ
                  LEFT JOIN vertex ON SUBQ.v1 = vertex.vertex_ID;"""
    # get all Vertex2 locations
    v2_query = """SELECT edge_ID, vertex_ID, chromosome, position, strand, gene_id FROM
                  (SELECT edge_ID,v2,chromosome,position,strand FROM edge
                   LEFT JOIN location ON edge.v2 = location.location_ID
                  WHERE genome_build = """ + str_wrap_double(build) + \
                  """ AND edge_type = """ + str_wrap_double(edge_type) + """) AS SUBQ
                  LEFT JOIN vertex ON SUBQ.v2 = vertex.vertex_ID;"""

    cursor.execute(v1_query)
    v1s = cursor.fetchall()
    cursor.execute(v2_query)
    v2s = cursor.fetchall()

    for v1, v2 in zip(v1s, v2s):
        edge = Edge.get_edge_from_db(v1, v2)
        if not (edge.start == edge.end):
            edge_tree.add_edge(edge)

    return edge_tree


def read_transcripts(cursor, exon_tree, intron_tree):
    """ Fetches transcripts, creates a Transcript object for each, and adds 
        them to a dictionary"""

    transcripts = {}

    cursor.execute('SELECT * FROM transcripts')
    transcript_rows = cursor.fetchall()
    for t in transcript_rows:
        transcript = Transcript.get_transcript_from_db(t, exon_tree, intron_tree)
        if transcript != None:
            transcripts[transcript.identifier] = transcript

    return transcripts

def str_wrap_double(string):
    """ Adds double quotes around the input string """
    return '"' + string + '"'

def process_sam_file(sam_file, dataset):
    """ Reads transcripts from a SAM file
        Args:
            sam_file: Path to the SAM file
        Returns:
            sam_transcripts: List of sam_transcript objects
    """

    sam_transcripts = []

    with open(sam_file) as sam:
        for line in sam:
            line = line.strip()

            # Ignore header
            if line.startswith("@"):
                continue

            sam = line.split("\t")

            # Only use uniquely mapped transcripts for now
            if sam[1] not in ["0", "16"]:
                continue

            # Only use reads that are >= 300 bp long
            if len(sam[9]) < 300:
                continue
            
            try: 
                sam_transcript = SamTranscript.get_sam_transcript(sam, dataset)
                sam_transcripts.append(sam_transcript)
            except:
                print "An error occurred while processing sam transcript " + \
                      sam[0] + ". Will skip this transcript."
        
    return sam_transcripts

def identify_sam_transcripts(sam_transcripts, gene_tree, transcripts, exon_tree, 
                             intron_tree, vertices, counter, dataset, 
                             novel_ids, abundance):
    """ Assign each sam transcript an annotated or a novel transcript identity
        Returns:
            Modified versions of gene_tree, transcripts, exon_tree, counter to
                which novel objects have been added
            abundance_dict: dictionary mapping transcript IDs to the number of 
                times each was observed in the sam data 
    """

    for sam_transcript in sam_transcripts:
        print sam_transcript.identifier
        chromosome = sam_transcript.chromosome
        start = sam_transcript.start
        end = sam_transcript.end
        strand = sam_transcript.strand

        # Look for full and partial matches
        best_match, match_type, edge_matches, diff = find_transcript_match(
                                                     sam_transcript, transcripts, 
                                                     exon_tree, intron_tree)

        # Full transcript match
        if match_type == "full":
            sam_transcript.gene_ID = best_match.gene_id
            sam_transcript.transcript_ID = best_match.identifier
            annot_transcript = best_match 
            sam_transcript.novel = "Seen_previously"

        # If there is no full match, a novel transcript must be created
        else:
            sam_transcript.novel = "Novel"

            # Search for gene using transcript coordinates
            gene_match_id = Vertex.search_for_gene(sam_transcript, vertices)     

            if gene_match_id == None:
                match_type = "None"
                chromosome = sam_transcript.chromosome
                start = sam_transcript.start
                end = sam_transcript.end
                strand = sam_transcript.strand
             
                gene_obj = Gene.create_novel_gene(chromosome, start, end,
                                                  strand, counter)
                sam_transcript.gene_ID = gene_obj.identifier 
                gene_tuple = (gene_obj.identifier, dataset)
                novel_ids['genes'][gene_obj.identifier] = gene_tuple
                gene_tree.add_gene(gene_obj)
                
            else:
                gene_obj = gene_tree.genes[gene_match_id]
                sam_transcript.gene_ID = gene_obj.identifier
            # Create the novel transcript and add to transcript dict
            # and novel tracker
            annot_transcript = make_novel_transcript(sam_transcript, transcripts, edge_matches, 
                                                     exon_tree, intron_tree, vertices,
                                                     counter, novel_ids)

            novel_tuple = (annot_transcript.identifier, annot_transcript.gene_id,
                           annot_transcript.get_edge_path(), dataset)
            novel_ids["transcripts"][annot_transcript.identifier] = novel_tuple
            transcripts[annot_transcript.identifier] = annot_transcript

        get_transcript_start_and_end_diffs(sam_transcript, annot_transcript,
                                       vertices, dataset, counter, novel_ids)    

        # Add transcript observation to abundance dict
        if annot_transcript.identifier in abundance:
            try:
                abundance[annot_transcript.identifier][dataset] += 1
            except:
                abundance[annot_transcript.identifier][dataset] = 1
        else:
            abundance[annot_transcript.identifier] = {}
            abundance[annot_transcript.identifier][dataset] = 1
            
    return

def get_transcript_start_and_end_diffs(sam_transcript, annot_transcript, 
                                       vertices, dataset, counter, novel_ids):
    """Must already have a completed match for the query transcript"""

    # Get IDs
    obs_start_id = str(counter["observed_starts"] + 1)
    counter["observed_starts"] += 1
    obs_end_id = str(counter["observed_ends"] + 1)
    counter["observed_ends"] += 1

    # Get the 5' and 3' difference between the sam transcript and the 
    # annotated transcript assigned to it
    a = [sam_transcript.start, sam_transcript.end]
    b = [annot_transcript.start, annot_transcript.end]
    diff_5, diff_3 = TMT.get_difference(a, b, sam_transcript.strand)
    sam_transcript.diff_5 = diff_5
    sam_transcript.diff_3 = diff_3

    # Now, get the vertices that the differences apply to
    chromosome = str(annot_transcript.chromosome)
    strand = sam_transcript.strand
    gene_id = sam_transcript.gene_ID
 
    first_vertex = Vertex.fetch_vertex(vertices, chromosome, 
                                    annot_transcript.start, gene_id).identifier
    last_vertex = Vertex.fetch_vertex(vertices, chromosome,
                                    annot_transcript.end, gene_id).identifier

    if strand == "+":
        observed_start = (obs_start_id, annot_transcript.identifier,
                          first_vertex, dataset, diff_5)
        observed_end = (obs_end_id, annot_transcript.identifier,
                          last_vertex, dataset, diff_3) 
    elif strand == "-":
        observed_start = (obs_start_id, annot_transcript.identifier,
                          last_vertex, dataset, diff_5)
        observed_end = (obs_end_id, annot_transcript.identifier,
                          first_vertex, dataset, diff_3)

    novel_ids['observed_starts'][obs_start_id] = observed_start
    novel_ids['observed_ends'][obs_end_id] = observed_end

    return 

def make_novel_transcript(sam_transcript, transcripts, edge_matches, exon_tree, 
                          intron_tree, vertices, counter, novel_ids):

    chromosome = sam_transcript.chromosome
    strand = sam_transcript.strand
    gene_id = sam_transcript.gene_ID
    novel_transcript_id = str(counter['transcripts'] + 1)
    dataset = sam_transcript.dataset

    # First, get edges in the novel transcript
    novel_transcript_exons = []
    novel_transcript_introns = []
    i = 0
    for sam_edge, edge_match in zip(sam_transcript.get_all_edges(), edge_matches):
        edge_start = sam_edge.start
        edge_end = sam_edge.end

        # Even indices are exons
        if i % 2 == 0:
            edge_tree = exon_tree
            edge_type = "exon"
            novel_list = novel_transcript_exons
        else:
            edge_tree = intron_tree
            edge_type = "intron"
            novel_list = novel_transcript_introns

        # To ensure that all edges in a transcript come from the same gene,
        # it is necessary to create a new edge if the match comes from a 
        # different gene
        match_gene_id = None        

        if edge_match in edge_tree.edges:
            edge_obj = edge_tree.edges[edge_match]
            match_gene_id = edge_obj.gene_id

        if match_gene_id != gene_id or edge_match == None:
            edge_obj = Edge.create_novel_edge(chromosome, edge_start,
                                              edge_end, strand, gene_id,
                                        sam_transcript.transcript_ID, counter)
                   
            Vertex.try_vertex_update(edge_obj, vertices, novel_ids, counter)
            novel_ids["edges"][edge_obj.identifier] = (edge_obj.identifier, 
                                                       edge_obj.v1,
                                                       edge_obj.v2,
                                                       edge_type, 
                                                       dataset)
                                                              
            edge_obj.transcript_ids.add(novel_transcript_id)
            edge_tree.add_edge(edge_obj)

        else:
            # Add the transcript ID to the exon
            edge_obj.transcript_ids.add(novel_transcript_id)   
 
        # Add edge to intron or exon list.
        novel_list.append(edge_obj)
        i += 1

    start = novel_transcript_exons[0].start
    end = novel_transcript_exons[-1].end
    annot_transcript = Transcript.create_novel_transcript(chromosome, 
                                          start, end, strand, gene_id, counter,
                                          novel_transcript_exons, 
                                          novel_transcript_introns) 
    sam_transcript.transcript_ID = annot_transcript.identifier
    return annot_transcript                        


def find_transcript_match(query_transcript, transcripts, exon_tree, intron_tree):
    """ Performs search for matches to the query transcript, one exon at a time.
        Args:
            query_transcript: Transcript object to be matched
            transcripts: Dictionary of transcript_id -> transcript object.
            This is the catalog of transcripts that we've seen before.
            exon_tree: ExonTree structure storing the exons that we have seen
            before.            
        Returns:
            transcript_match:  if the query is matched to a known transcript 
                in full or in part. None otherwise.
            diff: [5' diff, 3' diff] from full transcript match, [None, None]
                  for partial or no match
            match_type: "full", "partial", or "none"
    """
    # Find annotation matches where possible for all exons and introns
    transcript_match = "none"
    match_type = "none"
    diff = [None, None]
    tracker = TMT.MatchTracker(query_transcript)   
    tracker.match_all_edges(exon_tree, intron_tree) 

    # Find transcript matches
    tracker.compute_match_sets(transcripts)

    # If there is more than one full transcript match, select and return the
    # best one
    if len(tracker.full_matches) > 0:
        transcript_match, diff = tracker.get_best_full_match(transcripts)
        match_type = "full"  
        edge_matches = transcript_match.get_all_edges()     

    # If there are no full matches, look for partial matches
    else:
        if len(tracker.partial_matches) > 0:
            transcript_match = tracker.get_best_partial_match(transcripts)
            match_type = "partial"

        edge_matches = tracker.get_best_edge_matches()
    
    return transcript_match, match_type, edge_matches, diff
        

def update_database(database, datasets, transcripts, counter, novel_ids, 
                    abundances, batch_size, genome_build):
    """ Add novel entries to the supplied database
    """
    # Connecting to the database file
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    n_genes = str(len(novel_ids["genes"]))
    n_transcripts = str(len(novel_ids["transcripts"]))

    print "Adding " + n_genes + " novel genes to database..."
    batch_add_genes(cursor, novel_ids, batch_size)
        
    print "Adding " + n_transcripts + " novel transcripts to database..."
    batch_add_transcripts(cursor, novel_ids, batch_size)

    print "Adding edges and vertices to database............"
    batch_add_edges(cursor, novel_ids, batch_size)
    batch_add_vertices_and_locations(cursor, novel_ids, genome_build, batch_size)

    print "Adding observed starts and ends to database............"
    batch_add_observed_starts(cursor, novel_ids, batch_size)
    batch_add_observed_ends(cursor, novel_ids, batch_size)

    print "Adding datasets to database............"
    add_datasets(cursor, novel_ids, counter)

    print "Updating counter............."
    update_counter(cursor, counter)

    # Update abundance table
    batch_add_abundance(cursor, abundances, batch_size)

    # Before actually committing the changes to the database, perform a 
    # set of validity checks as a safeguard
    check_database_integrity(cursor)

    conn.commit()
    conn.close() 
    return

def check_database_integrity(cursor):
    """ Perform some checks on the database. Run before committing changes"""

    print "Validating database........"    
    # For each category, check that the number of table entries matches the counter
    counter_query = "SELECT * FROM counters"
    cursor.execute(counter_query)
    counters = cursor.fetchall()

    for table_name, curr_counter in counters:
        curr_counter = int(curr_counter)
        query = "select COUNT(*) from " + table_name
        cursor.execute(query)
        actual_count = int(cursor.fetchone()[0])

        if actual_count != curr_counter:
            print "table_count: "  + str(actual_count)
            print "counter_value: " + str(curr_counter)
            raise ValueError("Database counter for '" + table_name + \
                  "' does not match the number of entries in the table." + \
                  " Discarding changes to database and exiting...")

    return

def update_counter(cursor, counter):
    # Update the database counter
    update_g = 'UPDATE "counters" SET "count" = ? WHERE "category" = "genes"'
    cursor.execute(update_g,[counter['genes']])

    update_t = 'UPDATE "counters" SET "count" = ? WHERE "category" = "transcripts"'
    cursor.execute(update_t,[counter['transcripts']])

    update_e = 'UPDATE "counters" SET "count" = ? WHERE "category" = "edge"'
    cursor.execute(update_e,[counter['edges']])

    update_v = 'UPDATE "counters" SET "count" = ? WHERE "category" = "vertex"'
    cursor.execute(update_v,[counter['vertices']])

    update_d = 'UPDATE "counters" SET "count" = ? WHERE "category" = "dataset"'
    cursor.execute(update_d,[counter['datasets']])

    update_os = 'UPDATE "counters" SET "count" = ? WHERE "category" = "observed_transcript_start"'
    cursor.execute(update_os,[counter['observed_starts']])

    update_oe = 'UPDATE "counters" SET "count" = ? WHERE "category" = "observed_transcript_end"'
    cursor.execute(update_oe,[counter['observed_ends']])
  
    return

def batch_add_genes(cursor, novel_ids, batch_size):

    novel_tuples = novel_ids['genes'].values()
    gene_entries = []
    gene_annotations = []
    for nt in novel_tuples:
        gene_entries.append((nt[0],))
        gene_annotations.append((nt[0], "talon_run", nt[-1],
                                       "gene_status", "NOVEL"))

    index = 0
    while index < len(gene_entries):
        try:
            gene_batch = gene_entries[index:index + batch_size]
            annot_batch = gene_annotations[index:index + batch_size]
        except:
            gene_batch = gene_entries[index:]
            annot_batch = gene_annotations[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in 
                   ["gene_id"]]) + ") "
            command = 'INSERT INTO "genes"' + cols + "VALUES " + '(?)'
            cursor.executemany(command, gene_batch)
        except Exception as e:
            print(e)

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["ID", "annot_name", "source", "attribute", "value"]]) + ") "
            command = 'INSERT INTO "gene_annotations"' + cols + "VALUES " + \
                      '(?,?,?,?,?)'
            cursor.executemany(command, annot_batch)

        except Exception as e:
            print(e)

    return

def batch_add_transcripts(cursor, novel_ids, batch_size):

    # Using the novel IDs, extract transcripts that need to be added
    # and fetch their gene IDs and path (sequence of edges)

    novel_tuples = novel_ids['transcripts'].values()
    transcript_entries = []
    transcript_annotations = []
    for nt in novel_tuples:
        transcript_entries.append(nt[0:3])
        transcript_annotations.append((nt[0], "talon_run", nt[3], 
                                       "transcript_status", "NOVEL"))
    index = 0
    while index < len(transcript_entries):
        try:
            transcript_batch = transcript_entries[index:index + batch_size]
            annot_batch = transcript_annotations[index:index + batch_size]
        except:
            transcript_batch = transcript_entries[index:]
            annot_batch = transcript_annotations[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in 
                   ["transcript_id", "gene_id", "path"]]) +\
                    ") "
            command = 'INSERT INTO "transcripts"' + cols + "VALUES " + '(?,?,?)'
            cursor.executemany(command, transcript_batch)
        except Exception as e:
            print(e) 

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["ID", "annot_name", "source", "attribute", "value"]]) + ") "
            command = 'INSERT INTO "transcript_annotations"' + cols + \
                      "VALUES " + '(?,?,?,?,?)'
            cursor.executemany(command, annot_batch)

        except Exception as e:
            print(e)

    return

def batch_add_edges(cursor, novel_ids, batch_size):

    edge_tuples = novel_ids['edges'].values()
    edge_entries = []
    exon_annotations = []

    for entry in edge_tuples:
        edge_entries.append(entry[0:4])
        if entry[-2] == "exon":
            exon_annotations.append((entry[0], "talon_run", entry[-1],
                                       "exon_status", "NOVEL"))

    index = 0
    while index < len(edge_entries):
        try:
            batch = edge_entries[index:index + batch_size]
        except:
            batch = edge_entries[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["edge_ID", "v1", "v2", "edge_type"]]) + ") "
            command = 'INSERT INTO "edge"' + cols + "VALUES " + '(?,?,?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)

    index = 0
    while index < len(exon_annotations):
        try:
            annot_batch = exon_annotations[index:index + batch_size]
        except:
            annot_batch = exon_annotations[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["ID", "annot_name", "source", "attribute", "value"]]) + ") "
            command = 'INSERT INTO "exon_annotations"' + cols + \
                      "VALUES " + '(?,?,?,?,?)'
            cursor.executemany(command, annot_batch)

        except Exception as e:
            print(e)

    return

def batch_add_vertices_and_locations(cursor, novel_ids, genome_build, batch_size):

    novel_tuples = novel_ids['vertices'].values()
    vertex_entries = []
    location_entries = []
    for nt in novel_tuples:
        vertex_entries.append(nt[0:2])
        location_entries.append((nt[0], genome_build, nt[2], nt[3], nt[4]))

    index = 0
    while index < len(vertex_entries):
        try:
            vertex_batch = vertex_entries[index:index + batch_size]
            location_batch = location_entries[index:index + batch_size]
        except:
            vertex_batch = vertex_entries[index:]
            location_batch = location_entries[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["vertex_ID", "gene_id"]]) + ") "
            command = 'INSERT INTO "vertex"' + cols + "VALUES " + '(?,?)'
            cursor.executemany(command, vertex_batch)

            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["location_id", "genome_build", "chromosome", "position",
                    "strand"]]) + ") "
            command = 'INSERT INTO "location"' + cols + "VALUES " + \
                      '(?,?,?,?,?)'
            cursor.executemany(command, location_batch)

        except Exception as e:
            print(e)

    return

def batch_add_observed_starts(cursor, novel_ids, batch_size):
    observed = novel_ids['observed_starts'].values()
    index = 0
    while index < len(observed):
        try:
            batch = observed[index:index + batch_size]
        except:
            batch = observed[index:]
        index += batch_size

        # Add to database
        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["obs_start_ID", "transcript_ID", "vertex_ID", "dataset",
                     "delta"]]) + ") "
            command = 'INSERT INTO "observed_transcript_start"' + cols + \
                      "VALUES " + '(?,?,?,?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)
    return

def batch_add_observed_ends(cursor, novel_ids, batch_size):
    observed = novel_ids['observed_ends'].values()
    index = 0
    while index < len(observed):
        try:
            batch = observed[index:index + batch_size]
        except:
            batch = observed[index:]
        index += batch_size

    # Add to database
        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in
                   ["obs_end_ID", "transcript_ID", "vertex_ID", "dataset",
                     "delta"]]) + ") "
            command = 'INSERT INTO "observed_transcript_end"' + cols + \
                      "VALUES " + '(?,?,?,?,?)'
            cursor.executemany(command, batch)

        except Exception as e:
            print(e)
    return


def add_datasets(cursor, novel_ids, counter):
    datasets = novel_ids['datasets'].values()

    try:
        cols = " (" + ", ".join([str_wrap_double(x) for x in
               ["dataset_ID", "dataset_name", "sample", "platform"]]) + ") "
        command = 'INSERT INTO "dataset"' + cols + \
                  "VALUES " + '(?,?,?,?)'
        cursor.executemany(command, datasets)

    except Exception as e:
        print(e)
    return
        

def batch_add_abundance(cursor, abundance_dict, batch_size):

    abundances = []   

    for transcript_id in abundance_dict.keys():
        dataset_abundances = abundance_dict[transcript_id]
        for dataset in dataset_abundances.keys():
            abundance_tuple = (transcript_id, dataset, dataset_abundances[dataset])
            abundances.append(abundance_tuple)

    index = 0
    while index < len(abundances):
        try:
            batch = abundances[index:index + batch_size]
        except:
            batch = abundances[index:]
        index += batch_size

        try:
            cols = " (" + ", ".join([str_wrap_double(x) for x in 
                   ["transcript_id", "dataset", "count"]]) + ") "
            command = 'INSERT INTO "abundance"' + cols + "VALUES " + '(?,?,?)'
            cursor.executemany(command, batch)
        except Exception as e:
            print(e)
    return


def compute_abundances(sam_transcripts, dataset):
    """ Iterates over sam transcripts and records the transcript assignment of 
        each in a dictionary. This dictionary is then converted to a list of 
        tuples that can later be added to the TALON database"""
   
    abundance_dict = {}    
    for sam_transcript in sam_transcripts:
        transcript_id = sam_transcript.transcript_ID
        try:
            abundance_dict[transcript_id] += 1
        except:
            abundance_dict[transcript_id] = 1

    abundances = [(x, dataset, abundance_dict[x]) for x in abundance_dict.keys()]
    return abundances

def write_outputs(sam_transcripts, outprefix):
    """ Create a tab-delimited summary file that lists the gene and transcript
        assignments of every input transcript. """

    #out_sam = open(outprefix + "_talon.sam", 'w')
    out_txt = open(outprefix + "_talon.tsv", 'w')

    out_txt.write("\t".join(["dataset", "read_ID", "chromosome", "start", "end",
                       "strand", "gene_id", "transcript_id",
                       "annotation_status", "diff_5", "diff_3"]) + "\n") 

    for transcript in sam_transcripts:
        gene_id = transcript.gene_ID
        gene_name = "NA"
        transcript_id = transcript.transcript_ID
        annotation = transcript.novel
        chromosome = transcript.chromosome
        start = str(transcript.start)
        end = str(transcript.end)
        strand = transcript.strand
        length = str(transcript.get_length()) 
        diff_5 = str(transcript.diff_5)
        diff_3 = str(transcript.diff_3)
        dataset = str(transcript.dataset)

        try:
            out_txt.write("\t".join([dataset, transcript.identifier, \
                     chromosome, start, end, strand, gene_id, transcript_id, \
                     annotation, length, diff_5, diff_3]) + "\n") 
        except Exception as e: 
            print(e)
            
    out_txt.close()
    return

def write_abundance_file(transcripts, abundances, novel_ids, annot, datasets, 
                         outprefix):
    """ Create a tab-delimited summary file that lists the number of times
        each transcript was observed in each dataset from the run"""

    dataset_names = [ x[0] for x in datasets ]
    out_txt = open(outprefix + "_talon_abundance.tsv", 'w')
    columns = ["gene_ID", 
               "transcript_ID",
               "gene_name",
               "transcript_name",
               "gene_annotated",
               "gene_discovery_dataset", 
               "transcript_annotated",
               "transcript_discovery_dataset",
               "transcript_match_type",
               "n_exons"] + dataset_names
    out_txt.write("\t".join(columns) + "\n")

    # Using the database, get the annotation status of each gene and transcript
    gene_status = get_annotation_status_dict(annot, "gene")
    transcript_status = get_annotation_status_dict(annot, "transcript")
    exon_status = get_annotation_status_dict(annot, "exon")

    # Using the database, get the human-readable name of each gene and transcript
    gene_names = get_readable_name_dict(annot, "gene")
    transcript_names = get_readable_name_dict(annot, "transcript")

    for transcript_ID in abundances.keys():
        gene_ID = transcripts[transcript_ID].gene_id
        n_exons = str(len(transcripts[transcript_ID].exons))
        gene_transcript_info = get_transcript_info(transcript_ID, gene_ID, 
                                              gene_status, transcript_status)

        # Unpack info
        gene_annotated = gene_transcript_info[0]
        gene_discovery_dataset = gene_transcript_info[1]
        transcript_annotated = gene_transcript_info[2]
        transcript_discovery_dataset = gene_transcript_info[3] 

        # Get the type of transcript match
        match_type = get_match_type(transcripts[transcript_ID], gene_annotated,
                                    transcript_annotated, exon_status)

        # Get the gene and transcript names
        try:
            gene_name = gene_names[gene_ID]
        except:
            gene_name = "NA"
        try:
            transcript_name = transcript_names[transcript_ID]
        except:
            transcript_name = "NA"

        counts_by_dataset = abundances[transcript_ID]
        all_counts = []
        for dataset in dataset_names:
            if dataset in counts_by_dataset:
                count = counts_by_dataset[dataset]
            else:
                count = 0
            all_counts.append(str(count))

        output = [gene_ID, transcript_ID, gene_name, transcript_name] + \
                  gene_transcript_info + [match_type, n_exons] + all_counts
       
        out_txt.write("\t".join(output) + "\n")

    out_txt.close()
    return         

def get_match_type(transcript, gene_annotated, transcript_annotated, exon_status):

    if transcript_annotated == "KNOWN":
        match_type = "full_match"
    else:
        if gene_annotated == "KNOWN":
            exon_known = []
            for exon in transcript.exons:
                if exon.identifier in exon_status:
                    exon_known.append(exon_status[exon.identifier][0])
                else:
                    exon_known.append("AMBIGUOUS")

            if all([ x == "KNOWN" for x in exon_known]):
                match_type = "known_exons"
            else:
                match_type = "novel_exons"
        elif gene_annotated == "NOVEL":
            match_type = "novel_gene"
        else:
            match_type = "NA"
    return match_type

def get_transcript_info(transcript_ID, gene_ID, gene_status, transcript_status):
    """ Look up annotation info about the provided transcript and gene IDs. 
        If not found, set values to NA """

    try:
        gene_annotated = gene_status[gene_ID][0]
        gene_discovery_dataset = gene_status[gene_ID][1]
    except:
        gene_annotated = "NA"
        gene_discovery_dataset = "NA"

    try:
        transcript_annotated = transcript_status[transcript_ID][0]
        transcript_discovery_dataset = transcript_status[transcript_ID][1]
    except:
        transcript_annotated = "NA"
        transcript_discovery_dataset = "NA"

    return [gene_annotated, gene_discovery_dataset, transcript_annotated, 
            transcript_discovery_dataset]

def get_readable_name_dict(database, cat_type):
    """ Query the provided database and return a dictionary as follows:
            ID of gene or transcript -> [ KNOWN/NOVEL, SOURCE]
        The 'cat_type' variable should be either 'gene' or 'transcript',
        indicating which table to query. """

    annotation_names = {}

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    query = "SELECT * FROM " + cat_type + "_annotations WHERE attribute = '" + cat_type + "_name'"
    cursor.execute(query)
    annot_names = cursor.fetchall()

    for result in annot_names:
        ID = str(result[0])
        name = result[-1]
        annotation_names[ID] = name

    conn.close()
    return annotation_names

def get_annotation_status_dict(database, cat_type):
    """ Query the provided database and return a dictionary as follows:
            ID of gene or transcript -> [ KNOWN/NOVEL, SOURCE]
        The 'cat_type' variable should be either 'gene' or 'transcript',
        indicating which table to query. """

    annotation_status = {}

    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    query = "SELECT * FROM " + cat_type + "_annotations WHERE attribute = '" + cat_type + "_status'"
    cursor.execute(query)
    annot_status = cursor.fetchall()
   
    for result in annot_status:
        ID = str(result[0])
        annotation = result[-1]
        source = result[2]

        if ID not in annotation_status:
            annotation_status[ID] = [annotation, source]
        else:
            if annotation_status[ID][0] == "KNOWN":
                continue
            else:
                annotation_status[ID] = [annotation, source]

    conn.close()
    return annotation_status  

def filter_outputs_for_encode(curr_run_abundances, novel_ids, database):
    """ This function only gets run if 'encode' mode is enabled. It filters
        novel genes, transcripts, edges, vertices, and observed starts/stops
        so that only those belonging to novel transcripts observed in more
        than one dataset get added to the final database. """

    filtered_abundances = {}

    filtered_novel_ids = {'datasets': {}, \
                          'genes': {}, \
                          'transcripts': {}, \
                          'edges': {}, \
                          'vertices': {}, \
                          'observed_starts': {}, \
                          'observed_ends': {}}

    verified_genes = {}
    verified_edges = {}
    verified_vertices = {}
    database_abundances = {}

    # Find out how many datasets every transcript has been seen in
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    query = "SELECT transcript_ID FROM abundance"
    cursor.execute(query)
    
    for transcript_id in cursor.fetchall():
        transcript_id = str(transcript_id[0])
        try:
            database_abundances[transcript_id] += 1
        except:
            database_abundances[transcript_id] = 1
    conn.close()

    # Get novel/known status of transcripts
    transcript_status = get_annotation_status_dict(database, "transcript")

    for transcript_id in curr_run_abundances.keys():
        try:
            transcript_is_known = transcript_status[transcript_id][0] == "KNOWN"
        except:
            transcript_is_known = False

        if transcript_is_known or database_abundances[transcript_id] >= 2:
            filtered_abundances[transcript_id] = curr_run_abundances[transcript_id]

            if not transcript_is_known and transcript_id in novel_ids['transcripts']:
               filtered_novel_ids['transcripts'][transcript_id] = novel_ids['transcripts'][transcript_id]
               gene = novel_ids['transcripts'][transcript_id][1]
               edges = novel_ids['transcripts'][transcript_id][2]
               if gene in novel_ids['genes']:
                   verified_genes[gene] = 1

               # Whitelist the novel edges from this transcript
               for edge_id in edges.split(","):
                   if edge_id in novel_ids['edges']:
                       verified_edges[edge_id] = 1
                       v1 = novel_ids['edges'][edge_id][1]
                       v2 = novel_ids['edges'][edge_id][2]
                       verified_vertices[v1] = 1
                       verified_vertices[v2] = 1

    for gene_id in novel_ids['genes']:
        if gene_id in verified_genes:
            filtered_novel_ids['genes'][gene_id] = novel_ids['genes'][gene_id]

    for edge_id in novel_ids['edges']:
        if edge_id in verified_edges:
            filtered_novel_ids['edges'][edge_id] = novel_ids['edges'][edge_id]
        
    for vertex_id in novel_ids['vertices']:
        if vertex_id in verified_vertices:
            filtered_novel_ids['vertices'][vertex_id] = novel_ids['vertices'][vertex_id]

    for start in novel_ids['observed_starts']:
        transcript = novel_ids['observed_starts'][start][1]
        transcript_is_whitelisted = transcript in filtered_novel_ids['transcripts']
        transcript_is_known = transcript not in novel_ids['transcripts']
        if transcript_is_whitelisted or transcript_is_known:
            filtered_novel_ids['observed_starts'][start] = novel_ids['observed_starts'][start]

    for end in novel_ids['observed_ends']:
        transcript = novel_ids['observed_ends'][end][1]
        transcript_is_whitelisted = transcript in filtered_novel_ids['transcripts']
        transcript_is_known = transcript_id not in novel_ids['transcripts']
        if transcript_is_whitelisted or transcript_is_known:
            filtered_novel_ids['observed_ends'][end] = novel_ids['observed_ends'][end]

    return filtered_abundances, filtered_novel_ids 
    


def checkArgs(options):
    """ Makes sure that the options specified by the user are compatible with 
        each other """
    
    # Check that interval tree is working correctly
    test_tree = IntervalTree()
    test_tree[1:5] = 1
    if len(test_tree[1:5]) == 0:
        raise RuntimeError("Failed IntervalTree test. Make sure that " + \
                            "you are running a Python version >= 2.7.13.")

    config_file = options.config_file
    annot = options.annot
    build = options.build
    out = options.outfile

    # Make sure that the genome build exists in the provided TALON database.
    conn = sqlite3.connect(annot)
    cursor = conn.cursor() 
    cursor.execute(""" SELECT DISTINCT name FROM genome_build """)
    annot_builds = cursor.fetchone()
    if build not in annot_builds:
        build_names = ", ".join(list(annot_builds))
        raise ValueError("Please specify a genome build that exists in the" + 
                          " database. The choices are: " + build_names)
    annot_builds = cursor.fetchall()

    # Make sure that the dataset is not already in the database, and
    # also make sure that each dataset name is unique 
    sam_files = []
    dataset_metadata = []

    cursor.execute(""" SELECT dataset_name FROM dataset """)
    existing_datasets = [ str(x[0]) for x in cursor.fetchall() ]

    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            if len(line) != 4:
                raise ValueError('Incorrect number of comma-separated fields'+ \
                                 ' in config file. There should be four ' + \
                                 '(dataset name, sample description, ' + \
                                 'platform, associated transcript sam file).')

            metadata = (line[0], line[1], line[2])
            dataname = metadata[0]
            if dataname in existing_datasets and options.noUpdate == False:
                print "Ignoring dataset with name '" + dataname + "' because"+\
                      " it is already in the database."
            else:
                dataset_metadata.append(metadata)
                if not line[3].endswith(".sam"):
                    raise ValueError('Last field in config file must be a .sam file')
                sam_files.append(line[3])

    conn.close()

    # If the 'encode' flag is set, then the config file is required to have
    # at least two datasets in it (biological replicates)
    #if options.encode_mode:
    #    if len(sam_files) < 2:
    #        raise ValueError("When running TALON with 'encode' mode enabled,"+ \
    #                         " you must provide two or more biological " + \
    #                         "replicates in the config file.")

    return sam_files, dataset_metadata

def main():
    options = getOptions()

    # Check validity of input options
    sam_files, dataset_list = checkArgs(options)

    if len(sam_files) == 0:
        print "No new SAM files included in input. Exiting..."
        exit()

    annot = options.annot
    build = options.build
    out = options.outfile

    # Process the annotations
    print "Processing annotation...................."
    gene_tree, annot_transcripts, exon_tree, intron_tree, vertices, counter \
                                                 = read_annotation(annot, build)
    
    # Process the SAM files
    print "Processing SAM file......................"
    novel_ids = {'datasets': {}, \
                 'genes': {}, \
                 'transcripts': {}, \
                 'edges': {}, \
                 'vertices': {}, \
                 'observed_starts': {}, \
                 'observed_ends': {}}
                 
    all_sam_transcripts = []
    abundances = {}
    
    # Identify input sam transcripts
    for sam, d_metadata in zip(sam_files, dataset_list):

        # Create new dataset entry for the database
        d_id = str(counter['datasets'] + 1)
        novel_tuple = (d_id, d_metadata[0], d_metadata[1], d_metadata[2])
        d_name = d_metadata[0]
        novel_ids['datasets'][d_id] = novel_tuple
        counter['datasets'] += 1

        print "Identifying transcripts in " + d_name + "..............."
        sam_transcripts = process_sam_file(sam, d_name)
        if len(sam_transcripts) == 0:
            print "Warning: no transcripts detected in file " + sam
        identify_sam_transcripts(sam_transcripts, gene_tree, annot_transcripts, 
                                 exon_tree, intron_tree, vertices, counter, 
                                 d_name, novel_ids, abundances)
        
        all_sam_transcripts += sam_transcripts 

    # Update database
    if options.noUpdate == None:
        print "Updating TALON database.................."
        batch_size = 10000
        update_database(annot, dataset_list, annot_transcripts, counter,
                        novel_ids, abundances, batch_size, build)

    print "Writing SAM and summary file outputs..............."
    if options.encode_mode:
        abundances, novel_ids = filter_outputs_for_encode(abundances, novel_ids, 
                                                          annot)

    write_outputs(all_sam_transcripts, out)
    write_abundance_file(annot_transcripts, abundances, novel_ids, annot, 
                         dataset_list, out)


if __name__ == '__main__':
    main()
