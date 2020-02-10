from talon.post import filter_talon_transcripts as filt
import pandas as pd

def test_get_novelty_df():
    """ Should return a data frame like this:
        transcript_ID    transcript_novelty
              1               Known
              2               Known
              3               Genomic
              4               ISM
    """

    database = "scratch/filter/test.db"
    novelty = filt.get_novelty_df(database)
 
    assert list(novelty.iloc[0]) == [ 1, "Known"]
    assert list(novelty.iloc[1]) == [ 2, "Known"]
    assert list(novelty.iloc[2]) == [ 4, "ISM"]
    assert list(novelty.iloc[3]) == [ 3, "Genomic"]
   
def test_merge_reads_with_novelty():
    """ Should return a data frame like this:
        read_name  ... transcript_ID      transcript_novelty
          read_1   ...          1               Known
          read_2   ...          2               Known
          read_3   ...          3               Genomic
          read_4   ...          3               Genomic
          read_5   ...          3               Genomic
    """

    database = "scratch/filter/test.db"
    reads = pd.DataFrame([['read_1', 1], ['read_2', 2], ['read_3', 3],
                          ['read_4', 3], ['read_5', 3]],
                          columns = ["read_name", "transcript_ID"])
  
    novelty = filt.get_novelty_df(database)
    merged = filt.merge_reads_with_novelty(reads, novelty) 
    print(merged)

    assert len(merged) == 5
    assert list(merged.iloc[0]) == [ 'read_1', 1, "Known"]
    assert list(merged.iloc[1]) == [ 'read_2', 2, "Known"]
    assert list(merged.iloc[2]) == [ 'read_3', 3, "Genomic"]
    assert list(merged.iloc[3]) == [ 'read_4', 3, "Genomic"]
    assert list(merged.iloc[4]) == [ 'read_5', 3, "Genomic"]


