import sqlite3

def fetch_correct_ID(name, feature_type, cursor):
    """ Given the name of a transcript or gene, find its TALON ID in the
        database """

    query = """ SELECT ID FROM %s_annotations
                    WHERE attribute = '%s_name'
                    AND value = '%s' """

    cursor.execute(query % (feature_type, feature_type, name))
    feature_ID = cursor.fetchone()["ID"]

    return(feature_ID)

def get_db_cursor():
    conn = sqlite3.connect("scratch/toy.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    return conn, cursor
