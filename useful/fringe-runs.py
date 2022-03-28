#!/usr/bin/env python

"""Filters the HiPERCAM database file (hipercam.db) for runs suited to
making fringe maps. The database file can be found on the Warwick
hipercam log webpages.  Dumps results to a csv file which can be
loaded into excel or whatever. It's quick to run. As of 2022-02-08,
this returns just 54 runs making the jobs of finding runs suitable for
making fringes relatively easy. The file dumped is called
'fringe-runs.csv'.

This is an example of a script that selects runs given a set of
criteria that could easily be adapted for other purposes.

Selection criteria:

  1) >= 5 frames in the run
  2) Full frame readout
  3) No binning
  4) Nodding / dithering enabled
  5) Slide out or undefined
  6) Run type == data
  7) Cadence > 10 secs
  8) No binning
  9) Sun < -18 at start and end
 10) Moon phase < 0.3 or the Moon is below the horizon.

This returns relatively few runs; manual work after that..

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
        'dbases','hipercam.db'
    )
    conn = sqlite3.connect(dbase)

    # Query string where all the criteria are specified. NB
    # "hipercam" is the name of the table within hipercam.db

    query = """SELECT * FROM hipercam
WHERE nframe > 5 AND read_mode = 'FULL' AND
nodding = 'Y' AND (fpslide > 1000 OR fpslide = -99) AND
run_type = 'data' and cadence > 10 and binning = '1x1' AND
sun_alt_start < -18 AND sun_alt_end < -18 AND
(moon_phase < 0.3 OR (moon_alt_start < 0 AND moon_alt_end < 0))
"""

    # Runs query
    res = pd.read_sql_query(query, conn)

    # close connection
    conn.close()

    if len(res):
        # Found something; save to disk
        print(f'Found {len(res)} matching runs')
        res.to_csv('fringe-runs.csv')
    else:
        print('No matching runs found')
