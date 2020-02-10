from talon.post import filter_talon_transcripts as filt

def test_get_novelty_df():
    """ Should return a data frame like this:
        transcript_ID    transcript_novelty
              1               Known
              2               Known
              3               Genomic
    """

    database = "scratch/filter/test.db"
    novelty = filt.get_novelty_df(database)
 
    assert list(novelty.iloc[0]) == [ 1, "Known"]
    assert list(novelty.iloc[1]) == [ 2, "Known"]
    assert list(novelty.iloc[2]) == [ 3, "Genomic"]
    
