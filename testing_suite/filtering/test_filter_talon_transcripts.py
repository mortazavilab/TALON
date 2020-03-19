import os
import optparse_mock_filt as omf
from talon.post import filter_talon_transcripts as filt

def test_filter_lax():
    """ Working with our mock dataset, the very lax settings
          max_frac_A = 1; 
          min_count = 1; 
          min_datasets = 1 ; 
          allow_genomic = True;
        should give us all of the original transcripts:
             gene_ID  transcript_ID
                1           1
                1           2
                1           3
                1           4
    """

    database = "scratch/filter/test.db"
    annot = "toy"

    options = omf.OptParseMockFilt(database, annot, max_frac_A = 1,
                                   min_count = 1, min_datasets = 1, 
                                   allow_genomic = True) 

    datasets = ["dataset_1", "dataset_2", "dataset_3", "dataset_4", "dataset_5"]
    filtered = filt.filter_talon_transcripts(database, annot, datasets, options)
    print(filtered)

    assert len(filtered) == 4
    assert list(filtered.iloc[0]) == [1, 1]
    assert list(filtered.iloc[1]) == [1, 2]
    assert list(filtered.iloc[2]) == [1, 3]
    assert list(filtered.iloc[3]) == [1, 4]

def test_filter_keep_genomic():
    """ Working with our mock dataset, these settings
          max_frac_A = 0.5;
          min_count = 1;
          min_datasets = 2;
          allow_genomic = True;
        should give us the known transcripts (because they are known)
        and the genomic transcript, but not the ISM:
             gene_ID  transcript_ID
                1           1
                1           3
    """
    database = "scratch/filter/test.db"
    annot = "toy"

    options = omf.OptParseMockFilt(database, annot, max_frac_A = 0.5,
                                   min_count = 1, min_datasets = 2,
                                   allow_genomic = True)

    datasets = ["dataset_1", "dataset_2", "dataset_3", "dataset_4", "dataset_5"]
    filtered = filt.filter_talon_transcripts(database, annot, datasets, options)
    print(filtered)

    assert len(filtered) == 3
    assert list(filtered.iloc[0]) == [1, 1]
    assert list(filtered.iloc[1]) == [1, 2]
    assert list(filtered.iloc[2]) == [1, 3]

def test_filter_discard_genomic():
    """ Working with our mock dataset, these settings
          max_frac_A = 1;
          min_count = 2;
          min_datasets = 2;
          allow_genomic = False;
        should give us the known transcripts (because they are known)
        and ISM, but not the genomic:
             gene_ID  transcript_ID
                1           1
                1           4
    """
    database = "scratch/filter/test.db"
    annot = "toy"

    options = omf.OptParseMockFilt(database, annot, max_frac_A = 1,
                                   min_count = 2, min_datasets = 1,
                                   allow_genomic = False)

    datasets = ["dataset_1", "dataset_2", "dataset_3", "dataset_4", "dataset_5"]
    filtered = filt.filter_talon_transcripts(database, annot, datasets, options)
    print(filtered)

    assert len(filtered) == 3
    assert list(filtered.iloc[0]) == [1, 1]
    assert list(filtered.iloc[1]) == [1, 2]
    assert list(filtered.iloc[2]) == [1, 4]

