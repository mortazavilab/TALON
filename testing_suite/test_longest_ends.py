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



def check_lists(test, ctrl):
    assert test == ctrl

def get_test_annot():
    df = pd.read_csv('input_files/longest_ends/test_annot.tsv', sep='\t')
    return df
