# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from exon import *
from exontree import *
from transcript import *
from sam_transcript import *

class MatchTracker(object):
    """ Stores information to track full and partial matches to a query 
        transcript.

        Attributes:
            n_exons: Number of exons in the query transcript. Full matches must 
            contain exactly this number of exons.
 
            exon_matches: List of sets. The set at index i corresponds to 
            the annotation exon IDs that matched to exon[i]

            transcript_matches: List of sets. The set at index i corresponds to
            the annotation transcript IDs that matched to exon_matches[i]

            full_matches: collection of transcript matches that satisfy the 
            conditions to fully match the query (i.e. all exons match)

            partial_matches: collection of transcript IDs that were partial
            matches to the query transcript (i.e. some exons match)

    """
    def __init__(self, query_transcript):
        self.query_transcript = query_transcript
        self.n_exons = query_transcript.n_exons    
        self.exon_matches = [[] for i in range(self.n_exons)]
        self.transcript_matches = []

        self.full_matches = []
        self.partial_matches = []

    def match_all_exons(self, exon_tree):
        """ Iterates over provided exons and finds exon matches for each of 
            them where possible. These matches are stored in the Tracker object.
            As we go, keep track of which transcripts could match each exon as
            well. 
        """
        query_transcript = self.query_transcript
        chromosome = query_transcript.chromosome
        strand = query_transcript.strand
        n_exons = query_transcript.n_exons 
         
        
        for i in range(0, n_exons):
            q_exon = query_transcript.exons[i]
            cutoff_5, cutoff_3 = set_cutoffs_permissiveEnds(i, n_exons, strand)
 
            matches = get_exon_matches(q_exon, exon_tree, cutoff_5, cutoff_3)

            # Iterate over matching exons we found and pull transcript ID set
            # associated with each to add to the transcript match set
            transcript_matches = set()
            for match in matches:
                match_exon_id = match.obj_id
                transcript_ids = exon_tree.exons[match_exon_id].transcript_ids
                transcript_matches = transcript_matches.union(transcript_ids)

            self.exon_matches[i] = matches
            self.transcript_matches.append(transcript_matches)

        return


    def compute_match_sets(self, transcript_dict):
        """ Use the exon_matches field of the MatchTracker to figure out which
            transcript matches are full matches and which are partial.

            The transcript_dict is needed right now to check how many exons the
            match transcript has.
        """

        if len(self.transcript_matches) != self.n_exons:
            raise ValueError('Cannot run compute_match_sets until all ' + \
                             ' query exons have been processed.')

        tmp_full_matches = set.intersection(*(self.transcript_matches))
        full_matches = set()

        # Screen the full matches to make sure they have the correct number of
        # exons
        for match in tmp_full_matches:
            n_exons_match = len(transcript_dict[match].exons)
            if self.n_exons == n_exons_match:
                full_matches.add(match) 
        partial_matches = set.union(*(self.transcript_matches))^full_matches

        self.full_matches = list(full_matches)
        self.partial_matches = list(partial_matches)

        return

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
            by first computing the number of matching exons, and using 
            differences at the 3' and 5' ends as a tiebreaker. The purpose of
            selecting a partial match is mainly to choose the best gene match 
            for the query so that we know which gene to assign the novel 
            transcript to. It is not intended to be the end all be all of 
            transcript similarity since it doesn't consider the configuration
            of the exon matches.

            Best match: max number of matching exons
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
        best_n_exons = 0
        query_pos = [query.start, query.end]

        strand = query.strand
        for match_id in self.partial_matches:

            # Count the number of exons that match between query & partial match
            match = transcripts[match_id]
            match_pos = [match.start, match.end]
            diff_5, diff_3 = get_difference(query_pos, match_pos, strand)
            tot_diff = abs(diff_3) + abs(diff_5)
            shared_exons = 0
            for i in range(0,self.n_exons):
                if match_id in self.transcript_matches[i]:
                    shared_exons += 1

            if shared_exons > best_n_exons:
                best_match = match
                best_n_exons = shared_exons
                best_diff_3 = diff_3
                best_diff_5 = diff_5

            elif shared_exons == best_n_exons:
                # Apply tiebreaker rules
                tiebreakers = [ shared_exons > best_n_exons,
                    (diff_5 == diff_3 == 0),
                    (diff_3 == 0 and best_diff_3 != 0),
                    (diff_5 == 0 and best_diff_5 != 0 and best_diff_3 != 0),
                    (best_diff_3 != 0 and best_diff_5 != 0 and tot_diff < best_tot_diff) ]
                if any(tiebreakers):
                    best_match = match
                    best_n_exons = shared_exons
                    best_diff_3 = diff_3
                    best_diff_5 = diff_5

        return best_match

class Match(object):
    """ Describes the relationship of a query interval to an annotated object
        such as an exon or transcript.
 
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

def set_cutoffs_permissiveEnds(exon_index, n_exons, strand):
    """ Selects 3' and 5' difference cutoffs depending on the exon index and
        strand. For internal exons, both cutoffs should be set to zero, but
        more flexibility is allowed on the 3' and 5' ends.

        Args:
            exon_index: index we are at in the exon vector
            tot_exons: total number of exons in the transcript
            strand: strand that the exon is on.
    """
    cutoff_5 = 0
    cutoff_3 = 0

    if exon_index == 0:
        if strand == "+":
            cutoff_5 = 100000
        else:
            cutoff_3 = 100000
    if exon_index == n_exons-1:
        if strand == "+":
            cutoff_3 = 100000
        else:
            cutoff_5 = 100000
    return cutoff_5, cutoff_3


def get_exon_matches(query_exon, exon_tree, cutoff_5, cutoff_3):
    """ Finds and returns all annotation exons that overlap the query location
        and which start and end within the 3' and 5' cutoff of the query.

        Args:
            chromosome: Chromosome that the query is located on
            (format "chr1")

            start: The start position of the query with respect to the
            forward strand

            end: The end position of the query with respect to the
            forward strand

            strand: "+" if the query is on the forward strand, and "-" if
            it is on the reverse strand

            exons: An ExonTree object with exons organized by chromosome in an
            interval tree data structure

            cutoff_5: Permitted difference at 5' end

            cutoff_3: Permitted difference at 3' end

        Returns:
            matches: list of exon Match objects to query
    """
    matches = []

    # Get all annotation exons that overlap the query location
    chromosome = query_exon.chromosome 
    start = query_exon.start
    end = query_exon.end
    strand = query_exon.strand
    exon_matches = exon_tree.get_exons_in_range(chromosome, start, end, strand)
    
    # Compute the 5' and 3' differences
    for match in exon_matches:
        query_range = [start, end]
        match_range = [match.start, match.end]

        diff_5, diff_3 = get_difference(query_range, match_range, strand)

        # Enforce similarity cutoff
        if abs(diff_5) <= cutoff_5 and abs(diff_3) <= cutoff_3:

            # Create match object
            match = Match(chromosome, start, end, strand, match.identifier, 
                          diff_5, diff_3)
            matches.append(match)

    return matches
    

def get_difference(a, b, strand):
    """ Computes the 5' and 3' difference between two exon intervals.

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

