import timeit
import sqlite3
import sys
sys.path.append("..")
import talonQ as talon

def run_vertex_queries(location_dict, n):
    """ Run vertex query n times"""
    for i in range(0,n):
        talon.search_for_vertex_at_pos("chr1", 1, location_dict)

def run_timetest(cursor, n):
    fn = lambda: run_vertex_queries(cursor, n)
    out = "Time for %d query: %f seconds"
    print(out % (n, timeit.timeit(fn, number = 1)))

def main():
    #conn = sqlite3.connect("../qtests/scratch/toy.db")
    conn = sqlite3.connect("../../Temp_TALON_database_experiments/unmodified_full_gencode_v24_12-20-18.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    location_dict = talon.make_location_dict("hg38", cursor)
    fn = lambda: run_vertex_queries(location_dict, n)

    run_timetest(cursor, 1)
    run_timetest(cursor, 10)
    run_timetest(cursor, 100)
    run_timetest(cursor, 1000)
    run_timetest(cursor, 10000)
    run_timetest(cursor, 100000)
    

    conn.close()

if __name__ == '__main__':
    main()




