class OptParseMock(object):
     """Class to mimic option parser for talon_label_reads.py"""
     def __init__(self, sam_file, genome_file, threads = 2, fracA_range_size = 10,
                  tmp_dir = "tmp_label_reads", delete_tmp = False,
                  outprefix = "talon_prelabels"):
         self.sam_file = sam_file
         self.genome_file = genome_file
         self.threads = threads
         self.fracA_range_size = fracA_range_size
         self.tmp_dir = tmp_dir
         self.delete_tmp = delete_tmp
         self.outprefix = outprefix
