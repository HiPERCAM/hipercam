#!/usr/bin/env python

"""
Tries to identify problematic runs in the log database files.
"""

import sys
import re
import sqlite3
import os
import pandas as pd

import hipercam as hcam

if __name__ == '__main__':

    # Get database files.

    # Search in standard directory
    dbases_dir = os.path.join(
        os.environ.get(
            'HIPERCAM_ENV',
            os.path.join(os.environ["HOME"],'.hipercam')
        ),
        'dbases'
    )

    instruments = ('ultracam','ultraspec','hipercam')

    results = []
    for instrument in instruments:

        # connect to database
        dbase = os.path.join(dbases_dir, instrument + '.db')
        if not os.path.exists(dbase):
            print(f'Database = {dbase} does not exist')
            continue

        conn = sqlite3.connect(f"file:{dbase}?mode=ro",uri=True)

        # build query string
        query = f'SELECT * FROM {instrument}\n'
        query += """WHERE mjd_start IS NULL"""

        print(f'\nQuerying database {dbase} with SQL string:\n\n{query}')
        res = pd.read_sql_query(query, conn)
        if len(res):
            print(res)
            results.append(res)
        else:
            print('   no runs found')

        # close connection
        conn.close()

    if False:
        total = pd.concat(results,sort=False)
        total.to_csv(output)
