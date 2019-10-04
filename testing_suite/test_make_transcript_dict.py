import pytest
from talon import talon
import sqlite3
@pytest.mark.unit

class TestTranscriptDict(object):
    def test_all_transcripts(self):
        """ Get all transcripts in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        t_dict = talon.make_transcript_dict(cursor, build)
     
        conn.close()

        assert set(t_dict.keys()) == set([ frozenset([1,2,3,4,5]), frozenset([6,7,8]),
                                           frozenset([9,10,11]),
                                           frozenset([12,13,14,15,16]), 
                                           frozenset([17,18,19,20,21,22,23,24,25,26,27]),
                                           frozenset([21,22,23,28,27,29,30]),
                                           frozenset([31]) ])

    def test_interval_transcripts(self):
        """ Get all transcripts in the database within the specified interval """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        t_dict = talon.make_transcript_dict(cursor, build, chrom = "chr1",
                                          start = 1, end = 1000)

        conn.close()
        assert set(t_dict.keys()) == set([ frozenset([1,2,3,4,5]), 
                                           frozenset([6,7,8])])


def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor 
