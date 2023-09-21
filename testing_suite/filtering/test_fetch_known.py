import make_minimal_db_for_filtering as mmdb
import os
from talon.post import filter_talon_transcripts as filt

def test_get_known_transcripts_all_datasets():
    """ Make sure the get_known_transcripts function returns both known
        transcripts when datasets are not speicified """

    database = "scratch/filter/test.db"
    include_annot = False
    known = filt.get_known_transcripts(database, "toy",
                                       include_annot,
                                       datasets = None)
    assert list(known.gene_ID) == [1, 1]
    assert list(known.transcript_ID) == [1, 2]

def test_get_known_transcripts_dataset_1_include_annot():
    """ Make sure the get_known_transcripts function returns all known
        transcripts with the include_annot function """

    database = "scratch/filter/test.db"
    include_annot = False
    known = filt.get_known_transcripts(database, "toy",
                                       include_annot,
                                       datasets = None)
    assert list(known.gene_ID) == [1, 1]
    assert list(known.transcript_ID) == [1, 2]

def test_get_known_transcripts_specific_dataset():
    """ Now make sure the correct transcript is returned when the dataset is
        specified. """

    database = "scratch/filter/test.db"
    include_annot = False

    # Both datasets
    known = filt.get_known_transcripts(database, "toy",
                                       include_annot,
                                       datasets = ["dataset_1", "dataset_2"])
    assert list(known.gene_ID) == [1, 1]
    assert list(known.transcript_ID) == [1, 2]

    # Dataset 1
    known = filt.get_known_transcripts(database, "toy",
                                       include_annot,
                                       datasets = ["dataset_1"])
    assert list(known.iloc[0]) == [1, 1]
    assert len(known) == 1

    # Dataset 2
    known = filt.get_known_transcripts(database, "toy",
                                       include_annot,
                                       datasets = ["dataset_2"])
    assert list(known.iloc[0]) == [1, 2]
    assert len(known) == 1

    # Dataset 3
    known = filt.get_known_transcripts(database, "toy",
                                       include_annot,
                                       datasets = ["dataset_3"])
    assert len(known) == 0
