import pytest
import pandas as pd
import os
import subprocess
import sys
from talon.post import get_transcript_sjs as tsj
@pytest.mark.integration

class TestGetTranscriptSJs(object):

    def test_gtf_2_dfs(self):
        """ Create location, edge and transcript data frames from GTF file and 
            make sure they contain the correct entries """
        gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf" 
        ref_loc_df, ref_edge_df, ref_t_df = tsj.create_dfs_gtf(gtf_file)
        print(ref_loc_df)
        run_df_tests(ref_loc_df, ref_edge_df, ref_t_df)

    def test_db_to_dfs(self):
        """ Initialize a TALON database from the same GTF annotation. Then,
            create location, edge and transcript data frames from GTF file and
            make sure they contain the correct entries. """

        os.system("mkdir -p scratch/get_transcript_sjs")
        try:
            os.remove("scratch/get_transcript_sjs/talon.db")
        except:
            pass
        try:
            subprocess.check_output(
                ["talon_initialize_database",
                 "--f", "input_files/test_get_transcript_sjs_util/annot.gtf",
                 "--a",  "toy_annot",
                 "--l", "0",
                 "--g",  "toy_build", "--o", "scratch/get_transcript_sjs/talon"])
        except Exception as e:
            print(e)
            sys.exit("Database initialization failed on toy annotation")       
 
        database = "scratch/get_transcript_sjs/talon.db"
        ref_loc_df, ref_edge_df, ref_t_df = tsj.create_dfs_db(database)
        run_df_tests(ref_loc_df, ref_edge_df, ref_t_df)

    def test_subset_edges(self):
        """ Make sure that this function returns only exons or introns
            as indicated """
        edge_df = pd.DataFrame({'edge_type': ['exon', 'intron'],
                                'strand': ['+', '+'],
                                'v1': [ 0, 1 ], 'v2': [1, 2],
                                'edge_id': [ (0, 1), (1, 2) ],
                                'chrom': ['chr1', 'chr1'],
                                'start': [1, 100],
                                'stop': [100, 500]})

        intron_df = tsj.subset_edges(edge_df, mode='intron')
        exon_df = tsj.subset_edges(edge_df, mode='exon')
        assert len(intron_df) == len(exon_df) == 1
        assert list(intron_df.iloc[0]) == ['intron', '+', 1, 2, (1,2), 'chr1', 100, 500]
        assert list(exon_df.iloc[0]) == ['exon', '+', 0, 1, (0,1), 'chr1', 1, 100]

    def test_determine_sj_novelty_Known_intron(self):
        """ Test that chr1:100-500 gets classified as all known """
        gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
        ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'intron')

        query_gtf = "input_files/test_get_transcript_sjs_util/known.gtf"
        loc_df, edge_df, t_df = prep_gtf(query_gtf, 'intron')

        edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
        assert edge_df.iloc[0].start_known == True
        assert edge_df.iloc[0].stop_known == True
        assert edge_df.iloc[0].combination_known == True

    def test_determine_sj_novelty_Known_exon(self):
        """ Test that chr1:1-100 gets classified as all known """
        gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
        ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'exon')

        query_gtf = "input_files/test_get_transcript_sjs_util/known.gtf"
        loc_df, edge_df, t_df = prep_gtf(query_gtf, 'exon')

        edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
        assert edge_df.iloc[0].start_known == True
        assert edge_df.iloc[0].stop_known == True
        assert edge_df.iloc[0].combination_known == True
       
    def test_determine_sj_novelty_NIC_intron(self):
        """ Test that chr1:100-900 gets classified as having a known start and stop,
            but a novel combination """

        gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
        ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'intron')

        query_gtf = "input_files/test_get_transcript_sjs_util/intron_NIC.gtf"
        loc_df, edge_df, t_df = prep_gtf(query_gtf, 'intron')

        edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
        assert edge_df.iloc[0].start_known == True
        assert edge_df.iloc[0].stop_known == True
        assert edge_df.iloc[0].combination_known == False
        
    def test_determine_sj_novelty_NNC_intron_donor(self):
         """ Test that chr1:90-900 gets classified as having a known stop and
             novel start"""
 
         gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
         ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'intron')
 
         query_gtf = "input_files/test_get_transcript_sjs_util/intron_NNC_donor.gtf"
         loc_df, edge_df, t_df = prep_gtf(query_gtf, 'intron')
 
         edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
         assert edge_df.iloc[0].start_known == False
         assert edge_df.iloc[0].stop_known == True
         assert edge_df.iloc[0].combination_known == False

    def test_determine_sj_novelty_NNC_exon_end(self):
         """ Test that chr1:1-90 gets classified as having a known start and
             novel stop"""

         gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
         ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'exon')

         query_gtf = "input_files/test_get_transcript_sjs_util/intron_NNC_donor.gtf"
         loc_df, edge_df, t_df = prep_gtf(query_gtf, 'exon')

         edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
         exon = edge_df.loc[edge_df['start'] == 1].iloc[0]
         print(exon)
         assert exon.start_known == True
         assert exon.stop_known == False
         assert exon.combination_known == False    

    def test_determine_sj_novelty_NNC_intron_acceptor(self):
         """ Test that chr1:100-800 gets classified as having a known start and 
             novel stop"""

         gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
         ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'intron')

         query_gtf = "input_files/test_get_transcript_sjs_util/intron_NNC_acceptor.gtf"
         loc_df, edge_df, t_df = prep_gtf(query_gtf, 'intron')

         edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
         assert edge_df.iloc[0].start_known == True
         assert edge_df.iloc[0].stop_known == False
         assert edge_df.iloc[0].combination_known == False

    def test_determine_sj_novelty_NNC_exon_start(self):
         """ Test that chr1:800-1000 gets classified as having a known stop and
             novel start"""

         gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
         ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'exon')

         query_gtf = "input_files/test_get_transcript_sjs_util/intron_NNC_acceptor.gtf"
         loc_df, edge_df, t_df = prep_gtf(query_gtf, 'exon')

         edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
         exon = edge_df.loc[edge_df['start'] == 800].iloc[0]
         print(exon)
         assert exon.start_known == False
         assert exon.stop_known == True
         assert exon.combination_known == False

    def test_determine_sj_novelty_antisense(self):
         """ Test that chr1:600-1000 on - strand gets classified as all novel"""

         gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
         ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'intron')

         query_gtf = "input_files/test_get_transcript_sjs_util/intron_novel_antisense.gtf"
         loc_df, edge_df, t_df = prep_gtf(query_gtf, 'intron')

         edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
         
         assert edge_df.iloc[0].start_known == False
         assert edge_df.iloc[0].stop_known == False
         assert edge_df.iloc[0].combination_known == False

    def test_determine_exon_novelty_antisense(self):
         """ Test that chr1:1-1000 on - strand gets classified as all novel"""

         gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
         ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'exon')

         query_gtf = "input_files/test_get_transcript_sjs_util/antisense_exon.gtf"
         loc_df, edge_df, t_df = prep_gtf(query_gtf, 'exon')
         edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
         exon = edge_df.loc[edge_df['start'] == 100].iloc[0]

         assert edge_df.iloc[0].start_known == False
         assert edge_df.iloc[0].stop_known == False
         assert edge_df.iloc[0].combination_known == False

    def test_transcript_exon_assignment(self):
        """ Test that exon chr1:1-1000 (+) gets assigned only to transcripts
            1 and 2 """
        gtf_file = "input_files/test_get_transcript_sjs_util/annot.gtf"
        ref_loc_df, ref_edge_df, ref_t_df = prep_gtf(gtf_file, 'exon')

        query_gtf = "input_files/test_get_transcript_sjs_util/transcript_exon_assignment.gtf"
        loc_df, edge_df, t_df = prep_gtf(query_gtf, 'exon')
        edge_df = tsj.determine_sj_novelty(ref_edge_df, edge_df)
        edge_df =tsj.find_tids_from_sj(edge_df, t_df, mode = 'exon')
        print(edge_df)

 
def prep_gtf(gtf, mode):
    """ Wrapper for GTF processing steps used by get_transcript_sjs main """
    loc_df, edge_df, t_df = tsj.create_dfs_gtf(gtf)
    edge_df = tsj.add_coord_info(edge_df, loc_df)
    edge_df = tsj.subset_edges(edge_df, mode=mode)
    edge_df = tsj.format_edge_df(edge_df)

    return loc_df, edge_df, t_df
        
def run_df_tests(ref_loc_df, ref_edge_df, ref_t_df):
    """ Runs the location, edge, and transcript tests on dfs """
    # Make sure the correct positions made it into the location df
    expected_locs = [("chr1", 1), ("chr1", 100), ("chr1", 500), 
                     ("chr1", 600), ("chr1", 900), ("chr1", 1000), 
                     ("chr4", 1000), ("chr4", 4000), 
                     ("chr1", 1500), ("chr1", 2000)]
    assert len(ref_loc_df) == len(expected_locs)

    for item in expected_locs:
        chrom = item[0]
        loc = item[1]
        assert ((ref_loc_df.chrom == chrom) & (ref_loc_df.coord == loc)).any()

    # Make sure that the correct edges are in the edge df
    assert len(ref_edge_df.loc[ref_edge_df.edge_type == 'intron']) == 3
    assert len(ref_edge_df.loc[ref_edge_df.edge_type == 'exon']) == 6

    # Check that the transcript position paths are correct across the tables
    expected_paths = {"ENST01": [1, 100, 500, 600, 900, 1000],
                      "ENST07": [4000, 1000],
                      "ENST03": [2000, 1500, 1000, 900]}
    for tid in ["ENST01", "ENST07", "ENST03"]:
        v_path = list(ref_t_df.loc[ref_t_df.tid == tid].path)[0]
        pos_path = []
        for vertex in v_path :
            pos = int(ref_loc_df.loc[ref_loc_df.vertex_id == vertex].coord)
            pos_path.append(pos)
        assert pos_path == expected_paths[tid]
