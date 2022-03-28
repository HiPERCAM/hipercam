import sys
import re
import sqlite3
import os
import keyring
import getpass
import subprocess
import tempfile

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
    with the same columns), this searches for matching calibration
    files. This is to aid data export. It first searches for flat
    fields, then for the combined list of data files and flats, it
    searches for bias frames. The flats do not try to determine if
    the filter matches in the case of ultracam because the filters
    can be unreliable.

    It works by searching through the database .db files which are
    accessed by password in the same way vias access to your keyring
    as described in |logsearch|.

    Arguments::

       runs : str
          csv input file

       diff : int
          maximum time difference in days to allow between a frame and
          a matching calibration file. Mainly here to prevent
          excessive numbers of matches. The "night" column is
          used. Thus diff=0 only allows calibrations from the same
          night to be considered, but diff=1 allows any preceding or
          following nights too.

       output : str
          Name of CSV file to store the results. The results include
          the original runs along with matching flats and biases and
          biases for the flats as well. This is readable by oocalc for
          instance (UTF-8, semi-colon separators disabled).  """

    command, args = utils.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("runs", Cline.LOCAL, Cline.PROMPT)
        cl.register("diff", Cline.LOCAL, Cline.PROMPT)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        runs = cl.get_value(
            "runs", "input csv file of runs",
            cline.Fname('results', '.csv')
        )
        diff = cl.get_value(
            "diff", "maximum time difference for matching calibrations (days)", 10, 0
        )
        output = cl.get_value(
            "output", "name of spreadsheet of results ['none' to ignore]",
            cline.Fname('results', '.csv', cline.Fname.NEW), ignore="none"
        )

    # Read the runs into pandas
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
            if os.path.exists(fname):
                start_time = os.path.getmtime(fname)
            else:
                start_time =''
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

    # write runs to be checked to junk file. this is because
    # i can't get in memory option to work for some reason
    dbname = 'zzz_junk.db'
    cnx = sqlite3.connect(dbname)
    runs_df.to_sql(name='tab', con=cnx, if_exists='replace')
    cnx.commit()
    cnx.close()

    # Search for flat field frames
    results = []
    for dbase, dtable in dbases:

        # connect to big database
        conn = sqlite3.connect(f"file:{dbase}?mode=ro", uri=True)

        # Add database / table runs.tab representing the runs we wish
        # to search over
        cursor = conn.cursor()
        cursor.execute(f'ATTACH "{dbname}" as runs')

        # Build query string to locate matching bias frames.
        #
        # Designed to:
        #
        # 1) only return entries from big table
        # 2) only from matching observing runs
        # 3) should not be the same run
        # 4) should match read speed and binning
        # 5) have more than 10 frames
        # 6) have some indication by name, type or comment that it is a bias.
        # 7) main table referred to as m, runs we are looking for biases for as t.
        # 8) main table referred to as m, runs we are looking for biases for as t.

        query = f'SELECT DISTINCT m.* FROM main.{dtable} AS m\n'
        query += f"""INNER JOIN runs.tab AS t
ON (m.obs_run = t.obs_run AND m.instrument = t.instrument)
WHERE (m.night != t.night OR m.run_no != t.run_no)
AND m.binning = '1x1' AND m.nframe > 10
AND (m.target LIKE '%flat%' OR m.run_type = 'flat' OR m.comment LIKE '%flat%')
AND ABS(JULIANDAY(m.night)-JULIANDAY(t.night)) <= {diff}
AND ((m.sun_alt_start IS NULL OR (m.sun_alt_start > -15 AND m.sun_alt_start < 0))
OR (m.sun_alt_end IS NULL OR (m.sun_alt_end > -15 AND m.sun_alt_end < 0)))
"""
        if dtable == 'ultraspec':
            query += f

        print(f'Searching for flats in "{dbase}" with SQL code:\n\n{query}\n')

        res = pd.read_sql_query(query, conn)
        if len(res):
            print(res)
            results.append(res)

        # close connection
        conn.close()

    # Add the "flats" to the selected runs.
    dbname = 'zzz_junk.db'
    cnx = sqlite3.connect(dbname)
    for res in results:
        res.to_sql(name='tab', con=cnx, if_exists='append')
    cnx.commit()
    cnx.close()

    # Now search the whole lot for bias frames
    results = [runs_df,]
    for dbase, dtable in dbases:

        # connect to big database
        conn = sqlite3.connect(f"file:{dbase}?mode=ro", uri=True)

        # Add database / table runs.tab representing the runs we wish
        # to search over
        cursor = conn.cursor()
        cursor.execute(f'ATTACH "{dbname}" as runs')

        # Build query string to locate matching bias frames.
        #
        # Designed to:
        #
        # 1) only return entries from big table
        # 2) only from matching observing runs
        # 3) should not be the same run
        # 4) should match read speed and binning
        # 5) have more than 10 frames
        # 6) have some indication by name, type or comment that it is a bias.
        # 7) main table referred to as m, runs we are looking for biases for as t.
        # 8) main table referred to as m, runs we are looking for biases for as t.

        query = f'SELECT DISTINCT m.* FROM main.{dtable} AS m\n'
        query += f"""INNER JOIN runs.tab AS t
ON (m.obs_run = t.obs_run AND m.instrument = t.instrument)
WHERE (m.night != t.night OR m.run_no != t.run_no)
AND m.read_speed = t.read_speed AND m.binning = t.binning AND m.nframe > 10
AND (m.target LIKE '%bias%' OR m.run_type = 'bias' OR m.comment LIKE '%bias%')
AND ABS(JULIANDAY(m.night)-JULIANDAY(t.night)) <= {diff}"""

        print(f'Searching for biases in "{dbase}" with SQL code:\n\n{query}\n')

        res = pd.read_sql_query(query, conn)
        if len(res):
            print(res)
            results.append(res)

        # close connection
        conn.close()

    # Save the results.
    biglist = pd.concat(results,sort=False)
    biglist.to_csv(output)
