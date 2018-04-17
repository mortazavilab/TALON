# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------

class MatchTracker(object):
    """ Stores information to track full and partial matches to a query 
        transcript.

        Attributes:
            n_exons: Number of exons in the query transcript. Full matches must 
            contain exactly this number of exons.
 
            exon_matches: List of sets. The set at index i corresponds to 
            the transcript IDs that contain one or more of matches(exon[i])

            full_matches: collection of transcript matches that satisfy the 
            conditions to fully match the query (i.e. all exons match)

            partial_matches: collection of transcript IDs that were partial
            matches to the query transcript (i.e. some exons match)

    """
    #TODO: update description once I implement match class
    def __init__(self, n_exons):
        self.n_exons = n_exons
        self.exon_matches = []

        self.full_matches = []
        self.partial_matches = []

    def compute_match_sets(self, transcript_dict):
        """ Use the exon_matches field of the MatchTracker to figure out which
            transcript matches are full matches and which are partial.

            The transcript_dict is needed right now to check how many exons the
            match transcript has.
        """

        if len(self.exon_matches) != self.n_exons:
            raise ValueError('Cannot run compute_match_sets until all ' + \
                             ' query exons have been processed.')
        #for elem in self.exon_matches:
        #    print elem

        tmp_full_matches = set.intersection(*(self.exon_matches))
        full_matches = set()

        # Screen the full matches to make sure they have the correct number of
        # exons
        for match in tmp_full_matches:
            n_exons_match = len(transcript_dict[match].exons)/2
            if self.n_exons == n_exons_match:
                full_matches.add(match) 
        partial_matches = set.union(*(self.exon_matches))^full_matches

        self.full_matches = list(full_matches)
        self.partial_matches = list(partial_matches)

        return


