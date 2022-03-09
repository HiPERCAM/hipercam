import sys
import re
import sqlite3
import os
import keyring
import getpass
import subprocess

import numpy as np
import pandas as pd

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline
from hipercam.utils import target_lookup

__all__ = [
    "calsearch",
]

#######################################################################
#
# calsearch -- search for calibration runs matching a set of input runs
#
#######################################################################

def calsearch(args=None):
    description = \
    """``calsearch runs output``

    Given a csv file from |logsearch| (possibly with rows edited, but
    all the same columns), this searches for matching calibration
    files. This is to aid data export.

    It works by searching through the database .db files which are
    accessed by password in the same way vias access to your keyring
    as described in |logsearch|.

    Arguments::

       runs : str
          csv input file

       output : str
          Name of CSV file to store the results.

    """

    command, args = utils.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("runs", Cline.LOCAL, Cline.PROMPT)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        runs = cl.get_value(
            "runs", "input csv file of runs",
            cline.Fname('results', '.csv')
        )
        output = cl.get_value(
            "output", "name of spreadsheet of results ['none' to ignore]",
            cline.Fname('results', '.csv', cline.Fname.NEW), ignore="none"
        )

    runs_df = pd.read_csv(runs)
    
    # Get database files.

    # First create directory for them if need be
    dbases_dir = os.path.join(
        os.environ.get(
            'HIPERCAM_ENV',
            os.path.join(os.environ["HOME"],'.hipercam')
        ),
        'dbases'
    )
    os.makedirs(dbases_dir, 0o700, True)

    # Then download them. Passwords will be prompted and, if the
    # subsequent download is successful, will be stored in the
    # system keyring

    server = 'https://cygnus.astro.warwick.ac.uk/phsaap/'

    dbases = []
    for dbase in ('ultracam', 'ultraspec', 'hipercam'):

        pword = keyring.get_password("Data logs", dbase)
        prompted = False
        if pword is None:
            pword = getpass.getpass(f'{dbase} logs password: ')
            prompted = True

        # accumulate list of files and equivalent table names
        fname = os.path.join(dbases_dir, f'{dbase}.db')
        dbases.append((fname, dbase))

        if pword != "":
            # use 'curl' to download. Check timestamp to see if
            # file is updated.
            start_time = os.path.getmtime(fname)
            args = [
                'curl','-u', f'{dbase}:{pword}','-o',fname,
                '-z',fname,f'{server}/{dbase}/logs/{dbase}.db'
            ]
            result = subprocess.run(
                args, capture_output=True, universal_newlines=True
            )
            if result.returncode and not os.path.exists(fname):
                raise hcam.HipercamError(
                    f'Failed to download {dbase}.db. Return from curl:'
                    + 'stdout={result.stdout}, stderr={result.stderr}'
                )
            elif result.returncode:
                print(
                    f'Failed to download {dbase}.db. Will use old'
                    'local copy although it may be out of date'
                )
            elif prompted:
                # successful, will store password in the keyring
                keyring.set_password("Data logs", dbase, pword)
                print(f'Downloaded {dbase}.db')
                print(f' stored password for {dbase} in keyring')

            end_time = os.path.getmtime(fname)
            if end_time == start_time:
                print(f' {dbase}.db unchanged on server')
            else:
                print(f' {dbase}.db updated from server')
        else:
            print(f'No attempt to update {fname}')

    print()

    results = []
    for dbase, dtable in dbases:

        # identify runs for the instrument
        tab = runs_df[runs_df['instrument'] == dtable]

        # connect to database
        conn = sqlite3.connect(dbase)
        
        for index, row in tab.iterrows():
            print(f"Searching for runs matching: run {row['obs_run']}, night {row['night']}, run {row['run_no']}")
        
            # build query string
            query = f'SELECT * FROM {dtable}\n'
            query += f'WHERE (obs_run == "{row["obs_run"]}")'
            
            #conn.create_function("REGEXP", 2, regexp)
            #query += f'WHERE (REGEXP("{regex}",target) AND total > {tmin})'

            res = pd.read_sql_query(query, conn)
            if len(res):
                print(res)
            else:
                print('   no runs found')

        # close connection
        conn.close()

    total = pd.concat(results,sort=False)
    total.to_csv(output)

