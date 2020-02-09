from talon import talon_label_reads as tlr
import optparse_mock
import os

def test_main_fn():
    """ Run script all the way through from main and check the outfiles """
    
    sam_file = "talon_label_reads/test_inputs/plus_strand_read.sam"
    genome_file = "talon_label_reads/test_inputs/toy_genome.fa"
    tmp_dir = "scratch/test_main/tmp"
    outprefix = "scratch/test_main/test"
    options = optparse_mock.OptParseMock(sam_file, genome_file,
                           tmp_dir = tmp_dir, outprefix = outprefix)

    if os.path.exists("scratch/test_main"):
            os.system("rm -r %s" % ("scratch/test_main"))

    tlr.main(options=options)

    # Check that outfiles exist
    final_sam = "scratch/test_main/test_labeled.sam"
    final_log =  "scratch/test_main/test_read_labels.tsv"

    assert os.path.isfile(final_sam)
    assert os.path.isfile(final_log) 
