"""Command line script to generate meta data"""

import sys
import os
import time
import re
import warnings
import traceback
import argparse

import numpy as np
import pandas as pd

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ = [
    "jtrawl",
]

#############################
#
# jtrawl -- trawls for junk
#
#############################

def jtrawl(args=None):
    description = \
    """jtrawl

    Attempts to identify junk data frames. This command is to be run
    in the "raw_data" directory containing night-by-night directories
    of data for |hipercam|, ULTRACAM or ULTRASPEC. It looks for stats
    data in the meta subdirectories of night directories containing
    runs and attempts simple-minded identification of bad data basically
    to try to reduce the amount of space.

    """

    cwd = os.getcwd()
    if os.path.basename(cwd) != "raw_data":
        print("** hmeta must be run in a directory called 'raw_data'")
        print("hmeta aborted", file=sys.stderr)
        return

    if cwd.find("ultracam") > -1:
        instrument = "ULTRACAM"
        itype = 'U'
        source = 'ul'
        cnams = ('1','2','3')
    elif cwd.find("ultraspec") > -1:
        instrument = "ULTRASPEC"
        itype = 'U'
        source = 'ul'
        cnams = ('1',)
    elif cwd.find("hipercam") > -1:
        instrument = "HiPERCAM"
        itype = 'H'
        source = 'hl'
        cnams = ('1','2','3','4','5')
    else:
        print("** jtrawl: cannot find either ultracam, ultraspec or hipercam in path")
        print("hmeta aborted", file=sys.stderr)
        return

    warnings.filterwarnings('ignore')

    linstrument = instrument.lower()

    # Now the actual work. Next are regular expressions to match run
    # directories, nights, and run files
    nre = re.compile("^\d\d\d\d-\d\d-\d\d$")
    ure = re.compile("^run\d\d\d\.xml$")
    hre = re.compile("^run\d\d\d\d\.fits$")

    # Get list of night directories
    nnames = [
        nname
        for nname in os.listdir(".")
        if nre.match(nname)
        and os.path.isdir(nname)
    ]
    nnames.sort()

    if len(nnames) == 0:
        print("no night directories found", file=sys.stderr)
        print("hmeta aborted", file=sys.stderr)
        return

    if instrument == 'ULTRASPEC':
        MEDNAMS = ('median',)
        JUNKLIMS = (65000,)
    elif instrument == 'ULTRACAM':
        MEDNAMS = ('median_L','median_R')
        JUNKLIMS = (50000, 30000, 30000)
    elif instrument == 'HIPERCAM':
        MEDNAMS = ('median_E','median_F','median_G','median_H')
        JUNKLIMS = (65000,65000,65000,65000,65000)

    MINSPREAD = 10

    for nname in nnames:

        print(f"Night {nname}")

        # load all the run names
        if itype == 'U':
            runs = [run[:-4] for run in os.listdir(nname) if ure.match(run) and
                    os.path.exists(os.path.join(nname,run[:-4]+'.dat'))]
        else:
            runs = [run[:-5] for run in os.listdir(nname) if hre.match(run)]
        runs.sort()

        if len(runs) == 0:
            print(f' No runs with data found in {nname}; skipping')
            continue

        # directory for any meta info such as the times, stats
        meta = os.path.join(nname, 'meta')

        for run in runs:
            dfile = os.path.join(nname, run)

            # check for existence of data files
            all_ok = True
            for cnam in cnams:
                oname = os.path.join(meta, f'{run}_{cnam}.csv')
                if not os.path.exists(oname):
                    all_ok = False
                    break

            if not all_ok:
                print(f'{dfile}, skipping as could not find all stats files')
                continue

            # now read and check
            for cnam, jlim in zip(cnams,JUNKLIMS):
                oname = os.path.join(meta, f'{run}_{cnam}.csv')

                # read pandas dataframe
                table = pd.read_csv(oname)

                for mnam in MEDNAMS:
                    med = table[mnam]
                    if len(med) < 5:
                        break
                    spread = table['p95'] - table['p5']
                    if len(med[med < jlim]) > 2 and len(spread[spread > MINSPREAD]) > 2:
                        break
                else:
                    print(f'{dfile} is probably JUNK')
