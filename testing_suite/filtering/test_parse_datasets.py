import os
import pytest
from talon.post import filter_talon_transcripts as filt

def test_parse_datasets_None():
    """ Should return None"""

    datasets = filt.parse_datasets(None, None)
    assert datasets == None

def test_parse_datasets_str():
    """ Should return ['dataset_1', 'dataset_3'] """
    database = "scratch/filter/test.db"

    assert filt.parse_datasets("dataset_1,dataset_3", database) == ['dataset_1', 'dataset_3']
   
def test_parse_datasets_file():
    """ Should return ['dataset_1', 'dataset_2'] """
    database = "scratch/filter/test.db"
    dataset_file = "filtering/input_files/dset_file.txt"

    assert filt.parse_datasets(dataset_file, database) == ['dataset_1', 'dataset_2'] 

def test_parse_datasets_crash():
    """ Should crash when provided a dataset name that is not in the database """
    database = "scratch/filter/test.db"
    dataset_str = "dataset_1,crash"
    
    with pytest.raises(ValueError) as e:
        assert filt.parse_datasets(dataset_str, database)
    print(e)
    assert str(e.value) == ("Problem parsing datasets. The following names are "
                            "not in the database: 'crash'. \nValid dataset names: "
                            "'dataset_1, dataset_2, dataset_3'")
