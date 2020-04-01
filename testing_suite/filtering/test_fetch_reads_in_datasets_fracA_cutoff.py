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

def test_fetch_reads_null_fraction_As(capfd):
    """ In the provided database, all fraction_As values are None, so no reads 
        should be returned. Also make sure a specific printed warning shows up"""

    database = "scratch/toy_mod.db"
    reads = filt.fetch_reads_in_datasets_fracA_cutoff(database, None, 1)
    assert len(reads) == 0

    # check for our printed warning
    out, _ = capfd.readouterr()
    print(out)
    assert "No reads passed maxFracA cutoff. Is this expected?" in out

def test_fetch_reads_unlabelled_fraction_As(capfd):
    """ In the provided database, all fraction_As values are unlabelled; """
    """ should print an unlabelled dataset warning."""

    database = "scratch/toy_mod.db"
    reads = filt.fetch_reads_in_datasets_fracA_cutoff(database, None, 0.5)
    assert len(reads) == 0

    # check for our printed warning
    out, _ = capfd.readouterr()
    assert "Reads in dataset toy appear to be unlabelled. Only known transcripts will pass the filter." in out
