import pytest
from .helper_fns import get_db_cursor
from talon import talon, init_refs
import sqlite3

@pytest.mark.unit

class TestLocationDict(object):
    def test_all_locations(self):
        """ Get all locations in the database """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        l_dict = init_refs.make_location_dict(build, cursor)
     
        conn.close()

        chroms = l_dict.keys()
        
        # Get positions
        pos = []
        for chrom in chroms:
            pos += (list(l_dict[chrom].keys()))
        
        assert set(chroms) == set(["chr1", "chr2", "chr3", "chr4"])
        assert set(pos) == set([1, 100, 500, 600, 900, 1000, 1500, 2000, 5000, 
                                5550, 6000, 6500, 200, 400, 800, 1200, 1400, 
                                1600, 1800, 2200, 2400, 2600, 4000])

    def test_interval_locations(self):
        """ Get all locations in the database within the specified interval """

        conn, cursor = get_db_cursor()
        build = "toy_build"

        l_dict = init_refs.make_location_dict(build, cursor, chrom = "chr1",
                                          start = 1, end = 1000)

        conn.close()

        chroms = l_dict.keys()
        # Get positions
        pos = []
        for chrom in chroms:
            pos += (list(l_dict[chrom].keys()))

        assert set(chroms) == set(["chr1"])
        assert set(pos) == set([1, 100, 500, 600, 900, 1000])

