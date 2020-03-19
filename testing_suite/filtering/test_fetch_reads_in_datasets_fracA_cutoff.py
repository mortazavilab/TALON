import os
from talon.post import filter_talon_transcripts as filt

def test_fetch_reads_all_dsets():
    """ Should return all reads except read_2, which has fraction_A = 0.7"""

    database = "scratch/filter/test.db"
    reads = filt.fetch_reads_in_datasets_fracA_cutoff(database, None, 0.5)

    assert set(list(reads.read_name)) == set(["read_1", "read_3", "read_4", "read_5"])

def test_fetch_reads_dataset2():
    """ Should return read_4 but not read_2 """

    database = "scratch/filter/test.db"
    reads = filt.fetch_reads_in_datasets_fracA_cutoff(database, ["dataset_2"], 0.5)

    assert list(reads.read_name) == ["read_4"]

def test_fetch_reads_null_fraction_As():
    """ In the provided database, all fraction_As values are None, so no reads 
        should be returned. """

    database = "scratch/toy_mod.db"
    reads = filt.fetch_reads_in_datasets_fracA_cutoff(database, None, 1)
    assert len(reads) == 0
