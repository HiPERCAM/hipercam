#!/usr/bin/env python

"""Filters the ULTRASPEC database file (ultraspec.db) for any use
of the avalanche mode. Dumps results to a csv file which can be
loaded into excel or whatever.

Written by T.R.Marsh

"""

import os
import sqlite3
import pandas as pd
import numpy as np

import pandas as pd
from hipercam.utils import format_hlogger_table

if __name__ == '__main__':

    # connect to database
    dbase = os.path.join(
        os.environ.get(
            'HIPERCAM_ENV',
            os.path.join(os.environ["HOME"],'.hipercam')
        ),
        'dbases','ultraspec.db'
    )
    conn = sqlite3.connect(dbase)

    # Query string where all the criteria are specified. NB
    # "hipercam" is the name of the table within hipercam.db

    query = """SELECT * FROM ultraspec
WHERE output == 'A'
"""

    # Run query
    res = pd.read_sql_query(query, conn)

    # close connection
    conn.close()

    if len(res):
        # Found something; save to disk
        print(f'Found {len(res)} matching runs')
        res.to_csv('avalanche.csv')
        print('Saved results to avalanche.csv')
    else:
        print('No matching runs found')
