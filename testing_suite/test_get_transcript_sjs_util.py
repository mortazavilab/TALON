import pytest
import pandas
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
        run_df_tests(ref_loc_df, ref_edge_df, ref_t_df)

    def test_db_to_dfs(self):
        """ Initialize a TALON database from the same GTF annotation. Then,
            create location, edge and transcript data frames from GTF file and
            make sure they contain the correct entries. """

        os.system("mkdir -p scratch/get_transcript_sjs")
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
        
def run_df_tests(ref_loc_df, ref_edge_df, ref_t_df):
    """ Runs the location, edge, and transcript tests on dfs """
    # Make sure the correct positions made it into the location df
    expected_locs = [("chr1", 1), ("chr1", 100), ("chr1", 500), ("chr1", 600),
                     ("chr1", 900), ("chr1", 1000), ("chr4", 1000),
                     ("chr4", 4000), ("chr1", 1500), ("chr1", 2000),
                     ("chr1", 900), ("chr1", 1000)]
    print(ref_loc_df)
    assert 1 == 2
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
