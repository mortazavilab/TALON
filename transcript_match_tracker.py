# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

from exon import *
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
    #TODO: update description once I implement match class
    def __init__(self, n_exons):
        self.n_exons = n_exons
        self.exon_matches = [[] for i in range(n_exons)]
        self.transcript_matches = []

        self.full_matches = []
        self.partial_matches = []

    def add_exon_match(self, exon_index, chromosome, start, end, strand, \
                       match, diff_5, diff_3):
        """ Adds a record to the exon_match list that describes a match
            between the exon at that index and an exon in the annotation
        """

        obj_id = match.identifier
        
        new_exon_match = Match(chromosome, start, end, strand, obj_id, \
                                diff_5, diff_3)
        (self.exon_matches[exon_index]).append(new_exon_match)
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
        #for elem in self.exon_matches:
        #    print elem

        tmp_full_matches = set.intersection(*(self.transcript_matches))
        full_matches = set()

        # Screen the full matches to make sure they have the correct number of
        # exons
        for match in tmp_full_matches:
            n_exons_match = len(transcript_dict[match].exons)/2
            if self.n_exons == n_exons_match:
                full_matches.add(match) 
        partial_matches = set.union(*(self.transcript_matches))^full_matches

        self.full_matches = list(full_matches)
        self.partial_matches = list(partial_matches)

        return

class Match(object):
    """ Describes the relationship of a query interval to an annotated object
        such as an exon or transcript.
 
        Attributes:
            chromosome:
            start:
            end: 
            strand:
            obj_id: identifier associated with the annotated object
            diff_5:
            diff_3:
            tot_diff:

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

