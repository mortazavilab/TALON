import timeit
import sqlite3
import sys
sys.path.append("..")
import talonQ as talon

def run_vertex_queries(cursor, n):
    """ Run vertex query n times"""
    for i in range(0,n):
        talon.search_for_vertex_at_pos("toy_build", "chr1", 1, cursor)

def run_timetest(cursor, n):
    fn = lambda: run_vertex_queries(cursor, n)
    out = "Time for %d query: %f seconds"
    print(out % (n, timeit.timeit(fn, number = 1)))

def main():
    conn = sqlite3.connect("../qtests/scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    fn = lambda: run_vertex_queries(cursor, n)

    run_timetest(cursor, 1)
    run_timetest(cursor, 10)
    run_timetest(cursor, 100)
    run_timetest(cursor, 1000)
    run_timetest(cursor, 10000)
    run_timetest(cursor, 100000)
    

    conn.close()

if __name__ == '__main__':
    main()




