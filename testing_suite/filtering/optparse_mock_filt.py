class OptParseMockFilt(object):
     """Class to mimic option parser for filter_talon_transcripts.py"""
     def __init__(self, database, annot, max_frac_A = 0.5,
                  min_count = 2, min_datasets = 2, allow_genomic = False,
                  outprefix = "talon_filter"):
         self.database = database
         self.annot = annot
         self.max_frac_A = max_frac_A
         self.min_count = min_count
         self.min_datasets = min_datasets
         self.outprefix = outprefix
         self.allow_genomic = allow_genomic
