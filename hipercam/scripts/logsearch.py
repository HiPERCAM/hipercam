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
    "logsearch",
]

#######################################################################
#
# logsearch -- carries out a search for objects in data logs
#
#######################################################################

def logsearch(args=None):
    description = \
    """``logsearch target (dmax) regex (nocase) tmin noothers output``

    Searches for targets in the |hiper| and |ucam| logs (using the
    database files ending '.db'). It can carry out a coordinate lookup
    given a name and/or carry out a regular expression search. It uses
    the sqlite3 databases generated by |hlogger| which it will try to
    download from the log webpages hosted at Warwick. See below for some
    details on how it does this. 'logsearch' is the best way to look for
    runs on a given target.

    If a target name is entered, it will first be searched for in
    SIMBAD. If that fails, it will be searched for coordinates in the
    form "J123456.7-123456" or similar, so the latter is always the
    fallback for objects that don't appear in SIMBAD. It can also
    search by regular expression matching.

    Arguments::

       target : str
          Target name. On the command line, must be enclosed in quotes if it
          contains spaces. This will be used first to carry out a lookup in
          SIMBAD to find the RA and Dec. Failing this, it tries to identify
          coordinates from a final strength of the form JHHMMSS.S[+/-]DDMMSS
          Enter "none" to ignore.

       dmax : float
          Maximum distance from lookup position, arcminutes

       regex : str
          Regular expression to use to try to match target names in
          addition to the coordinate matching. "none" to ignore.

       nocase : bool [if regex is not "none"]
          True for case-insensitive matching, else case-sensitive used
          with regex

       tmin : float
          Minimum exposure duration (in seconds) to cut out short runs.

       noothers : bool
          True to ignore any runs called "Others" (non-GTO ULTRASPEC)

       output : str
          Name of CSV file to store the results. 'none' to
          ignore. Usually you will want to specify this since the
          results are too wide to print to screen. Assuming you do
          save the results, they are best viewed in an excel-type
          programme or topcat, or they can be read programatically
          into a pandas Dataframe using pd.read_csv('results.csv').
          Column names from all instruments are concatenated which for
          instance means that a column appropriate for hipercam, might
          be blank for ULTRACAM and vice versa. If reading into
          oocalc, make sure to switch off semi-colons as delimiters
          and use UTF-8.

    .. Note::

       The program will attempt to download the databases from the
       Warwick server. Since they are often updated, it will always
       check, but only actually download them if the server files are
       newer than the local versions. When it does find a newer
       version, then you may see a delay as the download happens,
       depending upon your network speed. The downloads are stored on
       your computer in $HOME/.hipercam/dbases [or
       $HIPERCAM_ENV/dbases if you have set HIPERCAM_ENV]. The
       databases on the Warwick server are on the password-protected
       log pages. Thus the first time you use logsearch, you will be
       asked for these, one each for the three cameras. The program
       will store these passwords in your system "keyring" in a folder
       called "Data logs" (look up the python module "keyring" to find
       out more) so once set, you won't need to re-enter them each time.
       If the passwords change (not often I hope), you will have to
       delete those stored in your keyring and re-enter them.

    """

    command, args = utils.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("target", Cline.LOCAL, Cline.PROMPT)
        cl.register("dmax", Cline.LOCAL, Cline.PROMPT)
        cl.register("regex", Cline.LOCAL, Cline.PROMPT)
        cl.register("nocase", Cline.LOCAL, Cline.PROMPT)
        cl.register("tmin", Cline.LOCAL, Cline.PROMPT)
        cl.register("noothers", Cline.LOCAL, Cline.PROMPT)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        target = cl.get_value(
            "target", "target name for simbad lookup ['none' to ignore]",
            "AR Sco", ignore="none"
        )

        if target is not None:
            dmax = cl.get_value(
                "dmax", "maximum distance from target [arcminutes]",
                12., 0.
            )

            regex = cl.get_value(
                "regex", "regular expression to match target name ['none' to ignore]",
                "none", ignore="none"
            )

        else:
            regex = cl.get_value(
                "regex", "regular expression to match target name",
                "ar\s*sco"
            )

        if regex is not None:
            nocase = cl.get_value(
                "nocase", "case insensitive match?", True
            )

        tmin = cl.get_value(
            "tmin", "minimum exposure duration for a run to be included [seconds]", -1.
        )
        noothers = cl.get_value(
            "noothers", "ignore runs called 'Others'", True
        )
        output = cl.get_value(
            "output", "name of spreadsheet of results ['none' to ignore]",
            cline.Fname('results', '.csv', cline.Fname.NEW), ignore="none"
        )

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
                start_time = ''
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
    if target is not None:
        name, ra, dec = target_lookup(target)
        if name == 'UNDEF':
            print(f'Coordinate lookup for "{target}" failed')
            exit(1)

        print(
            f'Coordinate lookup for "{target}" returned name = "{name}", RA [hrs] = {ra}, Dec [deg] = {dec}'
        )

        field = dmax/60.
        ra *= 15
        cdec = np.cos(np.radians(dec))
        ralo = ra - field/cdec
        rahi = ra + field/cdec
        declo = dec - field
        dechi = dec + field

    results = []
    for dbase, dtable in dbases:

        # connect to database
        conn = sqlite3.connect(f"file:{dbase}?mode=ro",uri=True)

        # build query string
        query = f'SELECT * FROM {dtable}\n'

        if target is not None:
            if dtable == 'ultracam':
                query += (
                    f"WHERE (ra_deg > {ralo} AND ra_deg < {rahi}\n"
                    f"AND dec_deg > {declo} AND dec_deg < {dechi}\n"
                    f"AND total > {tmin})\n"
                )
            else:
                query += (
                    f"WHERE (((ra_deg > {ralo} AND ra_deg < {rahi}\n"
                    f"AND dec_deg > {declo} AND dec_deg < {dechi}) OR\n"
                    f"(ra_tel_deg > {ralo} AND ra_tel_deg < {rahi}\n"
                    f"AND dec_tel_deg > {declo} AND dec_tel_deg < {dechi}))\n"
                    f"AND total > {tmin})\n"
                )

            if regex is not None:
                conn.create_function("REGEXP", 2, regexp)
                query += f'OR (REGEXP("{regex}",target) AND total > {tmin})\n'

        else:
            conn.create_function("REGEXP", 2, regexp)
            query += f"WHERE (REGEXP('{regex}',target) AND total > {tmin})\n"

        if noothers:
            query += f"AND obs_run != 'Others'\n"

        print(f'\nQuerying {dbase}\n')
        res = pd.read_sql_query(query, conn)
        if len(res):
            print(res)
            results.append(res)
        else:
            print('   no runs found')

        # close connection
        conn.close()

    if output is not None:
        total = pd.concat(results,sort=False)
        total.to_csv(output)

def regexp(expr, item):
    """Function to use in sqlite3 regex search"""
    reg = re.compile(expr,re.IGNORECASE)
    return reg.search(item) is not None
