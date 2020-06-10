import sqlite3


def query_db(db, selected, table, condition):
    con = sqlite3.connect(db)
    cur = con.cursor()
    cur.execute(f'SELECT {selected} FROM {table} WHERE {condition};')
    return cur.fetchall()
    
