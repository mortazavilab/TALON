import pytest
import pandas as pd
import os
import subprocess
import sys
from talon.post import call_longest_ends as cle
@pytest.mark.integration

class TestLongestEnds(object):

    # make sure get_longest ends works
    # * tes, novel, all datasets
    # * tss, novel, all datasets
    # * tes, known, all datasets
    # * tss, known, all datasets
    def test_get_longest(self):
        """ Make sure that get_longest_ends works for various settings """
        df = get_test_annot()

        # tes, novel, all datasets
        print('TES, novel, all datasets...')
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='novel', datasets='all')
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_end.tolist()
        ctrl_tids = ['t_1', 't_3']
        ctrl_ends = [7, 20]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

        # tss, novel, all datasets
        print('TSS, novel, all datasets...')
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='novel', datasets='all')
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_start.tolist()
        ctrl_tids = ['t_1', 't_3']
        ctrl_ends = [0, 25]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

        # tes, known+novel, all datasets
        print('TES, known+novel, all datasets...')
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='all', datasets='all')
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_end.tolist()
        ctrl_tids = ['t_1', 't_2', 't_3', 't_4']
        ctrl_ends = [7, 15, 20, 29]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

        # tss, known+novel, all datasets
        print('TSS, known+novel, all datasets...')
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='all', datasets='all')
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_start.tolist()
        ctrl_tids = ['t_1', 't_2', 't_3', 't_4']
        ctrl_ends = [0, 9, 25, 35]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

        # tes, novel, all datasets
        print('TES, novel, dataset a')
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='novel', datasets=['a'])
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_end.tolist()
        ctrl_tids = ['t_1', 't_3']
        ctrl_ends = [5, 20]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

        # tss, novel, all datasets
        print('TSS, novel, dataset a')
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='novel', datasets=['a'])
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_start.tolist()
        ctrl_tids = ['t_1', 't_3']
        ctrl_ends = [0, 25]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

        # tes, known+novel, dataset a
        print('TES, known+novel, dataset a')
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='all', datasets=['a'])
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_end.tolist()
        ctrl_tids = ['t_1', 't_2', 't_3', 't_4']
        ctrl_ends = [5, 15, 20, 30]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

        # tss, known+novel, dataset a
        print('TSS, known+novel, all dataset a')
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='all', datasets=['a'])
        tids = end_df.transcript_ID.tolist()
        ends = end_df.read_start.tolist()
        ctrl_tids = ['t_1', 't_2', 't_3', 't_4']
        ctrl_ends = [0, 10, 25, 35]
        check_lists(tids, ctrl_tids)
        check_lists(ends, ctrl_ends)

    def test_replace_gtf_ends(self):
        """ Make sure that replace_gtf_end_coords works for various settings """
        df = get_test_annot()

        # tes, novel, all datasets
        print('TES, novel, all datasets...')
        gtf_df = get_test_gtf()
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='novel', datasets='all')
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tes')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                      11,11,11,12,11,11,14,
                      20,20,24,20,
                      24,26,30,26,24,26,24]
        ctrl_stops = [7,7,2,7,
                     16,13,12,13,16,13,16,
                     25,25,25,23,
                     34,34,34,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        # tss, novel, all datasets
        print('TSS, novel, all datasets...')
        gtf_df = get_test_gtf()
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='novel', datasets='all')
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tss')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                      11,11,11,12,11,11,14,
                      20,20,24,20,
                      24,26,30,26,24,26,24]
        ctrl_stops = [9,9,2,9,
                     16,13,12,13,16,13,16,
                     25,25,25,23,
                     34,34,34,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        # tes, known+novel, all datasets
        print('TES, known+novel, all datasets...')
        gtf_df = get_test_gtf()
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='all', datasets='all')
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tes')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                      11,11,11,12,11,11,14,
                      20,20,24,20,
                      24,29,30,29,24,26,24]
        ctrl_stops = [7,7,2,7,
                     16,15,12,15,16,13,16,
                     25,25,25,23,
                     34,34,34,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        ctrl_starts = [0,0,0,3,
                       11,11,11,12,11,11,14,
                       20,20,24,20,
                       24,26,30,26,24,26,24]
        ctrl_stops = [9,9,2,9,
                     16,13,12,13,16,13,16,
                     25,25,25,23,
                     34,34,34,29,30,30,25]

        # tss, known+novel, all datasets
        print('TSS, known+novel, all datasets...')
        gtf_df = get_test_gtf()
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='all', datasets='all')
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tss')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                       9,9,9,12,11,11,14,
                       20,20,24,20,
                       24,26,30,26,24,26,24]
        ctrl_stops = [9,9,2,9,
                     16,13,12,13,16,13,16,
                     25,25,25,23,
                     35,35,35,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        # tes, novel, dataset a
        print('TES, novel, dataset a')
        gtf_df = get_test_gtf()
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='novel', datasets=['a'])
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tes')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                       11,11,11,12,11,11,14,
                       20,20,24,20,
                       24,26,30,26,24,26,24]
        ctrl_stops = [5,5,2,5,
                     16,13,12,13,16,13,16,
                     25,25,25,23,
                     34,34,34,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        # tss, novel, dataset a
        print('TSS, novel, dataset a')
        gtf_df = get_test_gtf()
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='novel', datasets=['a'])
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tss')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                       11,11,11,12,11,11,14,
                       20,20,24,20,
                       24,26,30,26,24,26,24]
        ctrl_stops = [9,9,2,9,
                     16,13,12,13,16,13,16,
                     25,25,25,23,
                     34,34,34,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        # tes, known+novel, dataset a
        print('TES, known+novel, dataset a')
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='all', datasets=['a'])
        gtf_df = get_test_gtf()
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tes')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                       11,11,11,12,11,11,14,
                       20,20,24,20,
                       24,30,30,30,24,26,24]
        ctrl_stops = [5,5,2,5,
                     16,15,12,15,16,13,16,
                     25,25,25,23,
                     34,34,34,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        # tss, known+novel, dataset a
        print('TSS, known+novel, all dataset a')
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='all', datasets=['a'])
        gtf_df = get_test_gtf()
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tss')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                       10,10,10,12,11,11,14,
                       20,20,24,20,
                       24,26,30,26,24,26,24]
        ctrl_stops = [9,9,2,9,
                     16,13,12,13,16,13,16,
                     25,25,25,23,
                     35,35,35,29,30,30,25]
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)

        # tss+tes, known+novel, all datasets
        print('TSS and TES, known+novel, all datasets')
        end_df = cle.get_longest_ends(df, how='tss',
            novelty='all', datasets='all')
        gtf_df = get_test_gtf()
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tss')
        end_df = cle.get_longest_ends(df, how='tes',
            novelty='all', datasets='all')
        gtf_df = cle.replace_gtf_end_coords(gtf_df, end_df, how='tes')
        starts = gtf_df.start.tolist()
        stops = gtf_df.stop.tolist()
        ctrl_starts = [0,0,0,3,
                       9,9,9,12,11,11,14,
                       20,20,24,20,
                       24,29,30,29,24,26,24]
        ctrl_stops = [7,7,2,7,
                     16,15,12,15,16,13,16,
                     25,25,25,23,
                     35,35,35,29,30,30,25]
        print(gtf_df[['gene_id', 'transcript_id', 'entry_type', 'start', 'stop']])
        check_lists(starts, ctrl_starts)
        check_lists(stops, ctrl_stops)





def check_lists(test, ctrl):
    assert test == ctrl

def get_test_annot():
    df = pd.read_csv('input_files/longest_ends/test_annot.tsv', sep='\t')
    return df

def get_test_gtf():
    df = pd.read_csv('input_files/longest_ends/test_gtf.gtf', sep='\t')
    return df
