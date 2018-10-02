# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

import edge as Edge
import edgetree as EdgeTree
import sam_transcript as SamTranscript
import transcript as Transcript
import pdb

class MatchTracker(object):
    """ Stores information to track full and partial matches to a query 
        transcript.
        Attributes:
            n_edges: Number of edges i the query transcript. Full matches must 
            contain exactly this number of edges.
 
            edge_matches: List of sets. The set at index i corresponds to 
            the annotation edge IDs that matched to edge[i]
            transcript_matches: List of sets. The set at index i corresponds to
            the annotation transcript IDs that matched to edge_matches[i]
            full_matches: collection of transcript matches that satisfy the 
            conditions to fully match the query (i.e. all edges match)
            partial_matches: collection of transcript IDs that were partial
            matches to the query transcript (i.e. some edges match)
    """
    def __init__(self, query_transcript):
        self.query_transcript = query_transcript
        self.n_edges = len(query_transcript.exons) + len(query_transcript.introns)    
        self.edge_matches = [[] for i in range(self.n_edges)]
        self.transcript_matches = []

        self.full_matches = []
        self.partial_matches = []

    def match_all_edges(self, exon_tree, intron_tree):
        """ Iterates over provided edges and finds edge matches for each of 
            them where possible. These matches are stored in the Tracker object.
            As we go, keep track of which transcripts could match each edge as
            well. 
        """
        query_transcript = self.query_transcript
        chromosome = query_transcript.chromosome
        strand = query_transcript.strand
        
        all_edges = query_transcript.get_all_edges()
        n_edges = len(all_edges)
        for i in range(0, n_edges):
            if i % 2 == 0: # exon
                edge_tree = exon_tree
            else:
                edge_tree = intron_tree

            q_edge = all_edges[i]
            cutoff_5, cutoff_3 = set_cutoffs_permissiveEnds(i, n_edges, strand)
            matches = get_edge_matches(q_edge, edge_tree, cutoff_5, cutoff_3)

            # Iterate over matching edges we found and pull transcript ID set
            # associated with each to add to the transcript match set
            transcript_matches = set()
             
            for match in matches:
                match_edge_id = match.obj_id
                transcript_ids = edge_tree.edges[match_edge_id].transcript_ids
                transcript_matches = transcript_matches.union(transcript_ids)
            self.edge_matches[i] = matches
            self.transcript_matches.append(transcript_matches)

        return    

    def compute_match_sets(self, transcript_dict):
        """ Use the edge_matches field of the MatchTracker to figure out which
            transcript matches are full matches and which are partial.
            The transcript_dict is needed right now to check how many edges the
            match transcript has.
        """

        if len(self.transcript_matches) != self.n_edges:
            raise ValueError('Cannot run compute_match_sets until all ' + \
                             ' query edges have been processed.')

        tmp_full_matches = set.intersection(*(self.transcript_matches))
        full_matches = set()

        # Screen the full matches to make sure they have the correct number of
        # edges
        for match in tmp_full_matches:
            #if match not in transcript_dict:
                
            n_edges_match = len(transcript_dict[match].get_all_edges())
            if self.n_edges == n_edges_match:
                full_matches.add(match) 
        partial_matches = set.union(*(self.transcript_matches))^full_matches

        self.full_matches = list(full_matches)
        self.partial_matches = list(partial_matches)

        return

    def get_best_edge_matches(self):
        """ Iterates over each each edge and compares the edge matches in order
            to find the best one for each. This is done by computing the 
            differences at the 3' and 5' ends. It isn't necessary to use
            different rules by edge context at this point because the edge 
            matches were already selected under those rules. 
            Logic:
                1) If diff_3 = diff_5 = 0, that is the best.
                2) Next best is diff_3 = 0.
                3) Next best is diff_5 = 0.
                4) After that, just minimize tot_diff
            This function returns a list, with each element consisting of an
            edge object (the best match), or None if there was no match.
            """
        best_matches = []
        
        for i in range(0,self.n_edges):
            curr_matches = self.edge_matches[i]
            if len(curr_matches) == 0:
                best_matches.append(None)
                continue
            best_match = None
            best_diff_3 = 1000000
            best_diff_5 = 1000000
            best_tot_diff = 1000000

            match_tuples = []
            for match in curr_matches:
                diff_5 = match.diff_5
                diff_3 = match.diff_3
                tot_diff = abs(diff_3) + abs(diff_5)
                match_diffs = (match.obj_id, diff_5, diff_3, tot_diff)
                match_tuples.append(match_diffs)

            # First, sort by 3' dist and secondarily, by 5' dist
            match_tuples.sort(key=lambda x: (abs(x[2]), abs(x[1])))
            if match_tuples[0][2] == 0: 
                best_match = match_tuples[0][0]
                best_matches.append(best_match)
                continue
            
            # Try sorting by 5' first and 3' second
            match_tuples.sort(key=lambda x: (abs(x[1]), abs(x[2])))
            
            if match_tuples[0][1] == 0:
                best_match = match_tuples[0][0]
                best_matches.append(best_match)
                continue

            # That failing, sort by tot_dist
            match_tuples.sort(key=lambda x: abs(x[3]))
            best_match = match_tuples[0][0]
            best_matches.append(best_match)
                
        return best_matches     
               


    def get_best_full_match(self, transcripts):
        """ Iterates over the full matches in the tracker and determines 
            which one is the best fit to the query transcript. This is done
            by computing the differences at the 3' and 5' ends for each.
            Logic: 
                1) If diff_3 = diff_5 = 0, that is the best. 
                2) Next best is diff_3 = 0.
                3) Next best is diff_5 = 0.
                4) After that, just minimize tot_diff
                    
            If there are no full matches, it returns None.
            Args: 
                transcripts: dictionary mapping transcript_id -> transcript
                object. Necessary in order to get from the transcript ids
                stored in the match tracker to the objects themselves. 
        """
        query = self.query_transcript
        if len(self.full_matches) == 0:
            return None, ["NA", "NA"]
        best_match = None
        best_diff_3 = 1000000
        best_diff_5 = 1000000
        best_tot_diff = 1000000
        query_pos = [query.start, query.end]
        
        strand = query.strand
        for match_id in self.full_matches:
            match = transcripts[match_id]
            match_pos = [match.start, match.end]
            diff_5, diff_3 = get_difference(query_pos, match_pos, strand)
            tot_diff = abs(diff_3) + abs(diff_5)
            
            if diff_5 == diff_3 == 0:
                return match, [diff_5, diff_3]
            elif diff_3 == 0:
                best_match = match
                best_diff_3 = diff_3
                best_diff_5 = diff_5
                best_tot_diff = tot_diff
            elif diff_5 == 0:
                if best_diff_3 != 0:
                    best_match = match
                    best_diff_3 = diff_3
                    best_diff_5 = diff_5
                    best_tot_diff = tot_diff
            elif tot_diff < best_tot_diff:
                best_match = match
                best_diff_3 = diff_3
                best_diff_5 = diff_5
                best_tot_diff = tot_diff
        return best_match, [best_diff_5, best_diff_3]

    def get_best_partial_match(self, transcripts):
        """ Iterates over the partial matches in the tracker and determines
            which one is the best fit to the query transcript. This is done
            by first computing the number of matching edges, and using 
            differences at the 3' and 5' ends as a tiebreaker. The purpose of
            selecting a partial match is mainly to choose the best gene match 
            for the query so that we know which gene to assign the novel 
            transcript to. It is not intended to be the end all be all of 
            transcript similarity since it doesn't consider the configuration
            of the edge matches.
            Best match: max number of matching edges
            Tiebreaker Logic:
                1) If diff_3 = diff_5 = 0, that is the best.
                2) Next best is diff_3 = 0.
                3) Next best is diff_5 = 0.
                4) After that, just minimize tot_diff
            If there are no partial matches, it returns None.
            Args:
                transcripts: dictionary mapping transcript_id -> transcript
                object. Necessary in order to get from the transcript ids
                stored in the match tracker to the objects themselves.
        """
        query = self.query_transcript
        if len(self.partial_matches) == 0:
            return None
        best_match = None
        best_diff_3 = 1000000
        best_diff_5 = 1000000
        best_tot_diff = 1000000
        best_n_edges = 0
        query_pos = [query.start, query.end]

        strand = query.strand
        for match_id in self.partial_matches:

            # Count the number of edges that match between query & partial match
            match = transcripts[match_id]
            match_pos = [match.start, match.end]
            diff_5, diff_3 = get_difference(query_pos, match_pos, strand)
            tot_diff = abs(diff_3) + abs(diff_5)
            shared_edges = 0
            for i in range(0,self.n_edges):
                if match_id in self.transcript_matches[i]:
                    shared_edges += 1

            if shared_edges > best_n_edges:
                best_match = match
                best_n_edges = shared_edges
                best_diff_3 = diff_3
                best_diff_5 = diff_5

            elif shared_edges == best_n_edges:
                # Apply tiebreaker rules
                tiebreakers = [ shared_edges > best_n_edges,
                    (diff_5 == diff_3 == 0),
                    (diff_3 == 0 and best_diff_3 != 0),
                    (diff_5 == 0 and best_diff_5 != 0 and best_diff_3 != 0),
                    (best_diff_3 != 0 and best_diff_5 != 0 and tot_diff < best_tot_diff) ]
                if any(tiebreakers):
                    best_match = match
                    best_n_edges = shared_edges
                    best_diff_3 = diff_3
                    best_diff_5 = diff_5

        return best_match

class Match(object):
    """ Describes the relationship of a query interval to an annotated object
        such as an edge or transcript.
 
        Attributes:
            chromosome: Chromosome of query and match
            start: start position of query
            end: end position of query
            strand: strand of query and match
            obj_id: identifier associated with the annotated object
            diff_5: 5' end difference between query and transcript
            diff_3: 3' end difference between query and transcript
            tot_diff: total difference between query and transcript
    """
    def __init__(self, chromosome, start, end, strand, obj_id, diff_5, diff_3):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.obj_id = obj_id
        self.diff_5 = diff_5
        self.diff_3 = diff_3
        self.tot_diff = abs(self.diff_3) + abs(self.diff_5)

#--------------------- Functions ------------------------------------------

def set_cutoffs_permissiveEnds(edge_index, n_edges, strand):
    """ Selects 3' and 5' difference cutoffs depending on the edge index and
        strand. For internal edges, both cutoffs should be set to zero, but
        more flexibility is allowed on the 3' and 5' ends.
        Args:
            edge_index: index we are at in the edge vector
            tot_edges: total number of edges in the transcript
            strand: strand that the edge is on.
    """
    cutoff_5 = 0
    cutoff_3 = 0

    if edge_index == 0:
        if strand == "+":
            cutoff_5 = 100000
        else:
            cutoff_3 = 100000
    if edge_index == n_edges-1:
        if strand == "+":
            cutoff_3 = 100000
        else:
            cutoff_5 = 100000
    return cutoff_5, cutoff_3


def get_edge_matches(query_edge, edge_tree, cutoff_5, cutoff_3):
    """ Finds and returns all annotation edges that overlap the query location
        and which start and end within the 3' and 5' cutoff of the query.
    """
    matches = []
    # Get all annotation edges that overlap the query location
    chromosome = query_edge.chromosome 
    start = query_edge.start
    end = query_edge.end
    strand = query_edge.strand
    edge_matches = edge_tree.get_edges_in_range(chromosome, start, end, strand)
    
    # Compute the 5' and 3' differences
    for match in edge_matches:
        query_range = [start, end]
        match_range = [match.start, match.end]

        diff_5, diff_3 = get_difference(query_range, match_range, strand)

        # Enforce similarity cutoff
        if abs(diff_5) <= cutoff_5 and abs(diff_3) <= cutoff_3:

            # Create match object
            match_obj = Match(chromosome, start, end, strand, match.identifier, 
                          diff_5, diff_3)
            matches.append(match_obj)
            
    return matches
    

def get_difference(a, b, strand):
    """ Computes the 5' and 3' difference between two edge intervals.
        Example: a = [ 0, 10]  b = [ 2, 8 ] on + strand
            5' difference =  -2
            3' difference =  +2
        Args:
            a: First interval, formattted as a list
            b: Second interval, formatted as a list
    """
    if strand == "+":
        diff_5 = a[0] - b[0]
        diff_3 = a[1] - b[1]

    elif strand == "-":
        diff_5 = b[1] - a[1]
        diff_3 = b[0] - a[0]

    return diff_5, diff_3

def get_overlap(a, b):
    """ Computes the amount of overlap between two intervals.
        Returns 0 if there is no overlap. The function treats the start and
        ends of each interval as inclusive, meaning that if a = b = [10, 20],
        the overlap reported would be 11, not 10.
        Args:
            a: First interval, formattted as a list
            b: Second interval, formatted as a list
    """
    overlap = max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
    return overlap
