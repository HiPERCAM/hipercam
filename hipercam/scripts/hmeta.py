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
    "hmeta",
]

#############################
#
# hmeta -- generate meta data
#
#############################


def hmeta(args=None):
    description = \
    """hmeta

    This command is to be run in the "raw_data" directory containing
    night-by-night directories of data for |hipercam|, ULTRACAM or
    ULTRASPEC. It attempts to generate meta data on any run data it
    can find and write this to a file 'statistics.csv' in a
    sub-directory called meta. These data can be picked up by logging
    scripts. The sort of data it produces are means, medians, etc of
    the frames (or some of the frames -- up to a maximum of 100-200)
    of each run. The program only considers genuine frames, i.e. it
    copes with the nblue and nskip options of ULTRACAM and |hipercam|.
    It also know about the different amplifier outputs of the
    instruments and starts by subtracting the medians of each
    amplifier output to avoif being disturbed by variable mean bias
    offsets.

    hmeta takes a good while to run, so be nice when running it.

    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        dest="full",
        action="store_true",
        help="carry out full re-computation of stats for all valid nights",
    )
    args = parser.parse_args()

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
        print("** hmeta: cannot find either ultracam, ultraspec or hipercam in path")
        print("hmeta aborted", file=sys.stderr)
        return

    warnings.filterwarnings('ignore')

    linstrument = instrument.lower()

    # Now the actual work.  Next are regular expressions to match run
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

        # create directory for any meta info such as the times
        meta = os.path.join(nname, 'meta')
        os.makedirs(meta, exist_ok=True)

        # name of stats file
        stats = os.path.join(meta, 'statistics.csv')
        if not args.full and os.path.exists(stats):
            # if file already present and we are re-doing things
            # in full, don't attempt to re-compute
            continue

        # Accumulate results in an array
        barr = []
        for run in runs:
            dfile = os.path.join(nname, run)

            try:
                if itype == 'U':
                    rdat = hcam.ucam.Rdata(dfile)
                else:
                    rdat = hcam.hcam.Rdata(dfile, 1, False, False)
                print(f"  {dfile}")
            except hcam.ucam.PowerOnOffError:
                print(f'  {dfile} -- power on/off; skipping')
                continue
            except:
                # some other failure
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_tb(
                    exc_traceback, limit=1, file=sys.stderr
                )
                traceback.print_exc(file=sys.stderr)
                print(f'  {dfile} -- problem occurred; skipping')
                continue

            # For speed, analyse a maximum of 100-200
            # images of each CCD from any given run. Have to take
            # into account the skips / nblue parameters.
            ntotal = rdat.ntotal()
            if ntotal == 0:
                print(f'  {dfile} -- zero frames; skipping')
                continue

            ncframes = {}
            if instrument == 'ULTRASPEC':
                # just the one CCD here
                nstep = ntotal // min(ntotal, 100)
                ncframes['1'] = list(range(1,ntotal+1,nstep))
                nframes = ncframes['1']
                ns = {'1' : 0}

            elif instrument == 'ULTRACAM':
                # CCD 1, 2 read out each time, but 3 can be skipped
                nstep = ntotal // min(ntotal, 100)
                ncframes['1'] = list(range(1,ntotal+1,nstep))
                ncframes['2'] = ncframes['1']

                nb = rdat.nblue
                nbstep = nb*max(1, ntotal // min(ntotal, 100*nb))
                ncframes['3'] = list(range(nb,ntotal+1,nbstep))
                nframes = sorted(set(ncframes['1']+ncframes['3']))

                ns = {'1' : 0, '2' : 0, '3' : 0}

            else:
                raise NotImplementedError('HiPERCAM case not done yet')

            # define arrays for holding the stats
            medians, means, p1s, p16s, rmsps, p84s, p99s = {}, {}, {}, {}, {}, {}, {}
            for cnam, ncframe in ncframes.items():
                if instrument == 'ULTRASPEC':
                    medians[cnam] = np.empty_like(ncframe,dtype=np.float)
                elif instrument == 'ULTRACAM':
                    medians[cnam] = {
                        'L' : np.empty_like(ncframe,dtype=np.float),
                        'R' : np.empty_like(ncframe,dtype=np.float)
                    }
                else:
                    raise NotImplementedError('HiPERCAM case not done yet')

                means[cnam] = np.empty_like(ncframe,dtype=np.float)
                p1s[cnam] = np.empty_like(ncframe,dtype=np.float)
                p16s[cnam] = np.empty_like(ncframe,dtype=np.float)
                p84s[cnam] = np.empty_like(ncframe,dtype=np.float)
                rmsps[cnam] = np.empty_like(ncframe,dtype=np.float)
                p99s[cnam] = np.empty_like(ncframe,dtype=np.float)

            # now access the data and calculate stats. we have to remember that
            # ultracam and hipercam have different readout amps so we subtract
            # median values calculated for each amp separately, which is a little
            # painful.
            for n, nf in enumerate(nframes):
                mccd = rdat(nf)
                for cnam, ncframe in ncframes.items():
                    if nf in ncframe:
                        ccd = mccd[cnam]
                        nc = ns[cnam]
                        if instrument == 'ULTRASPEC':
                            medval = ccd.median()
                            ccd -= medval
                            medians[cnam][nc] = medval
                        elif instrument == 'ULTRACAM':
                            wl = hcam.Group(hcam.Window)
                            wr = hcam.Group(hcam.Window)
                            for nw, wnam in enumerate(ccd):
                                if nw % 2 == 0:
                                    wl[wnam] = ccd[wnam]
                                else:
                                    wr[wnam] = ccd[wnam]
                            ccdl = hcam.CCD(wl,ccd.nxtot,ccd.nytot)
                            medl = ccdl.median()
                            ccdl -= medl
                            ccdr = hcam.CCD(wr,ccd.nxtot,ccd.nytot)
                            medr = ccdr.median()
                            ccdr -= medr
                            medians[cnam]['L'][nc] = medl
                            medians[cnam]['R'][nc] = medr
                        else:
                            raise NotImplementedError('HiPERCAM case not done yet')

                        # At this stage the median value should have been subtracted
                        # from the CCD on a per output basis. Remaining stats calculated
                        # from these median subtracted images.
                        means[cnam][nc] = ccd.mean()
                        p1,p16,p84,p99 = ccd.percentile([1.0,15.865,84.135,99.0])
                        p1s[cnam][nc] = p1
                        p16s[cnam][nc] = p16
                        rmsps[cnam][nc] = (p84-p16)/2
                        p84s[cnam][nc] = p84
                        p99s[cnam][nc] = p99
                        ns[cnam] += 1

            # All extracted from the run; take and store medians
            # of the extracted stats
            brow = [run,]
            for cnam in cnams:
                if len(means[cnam]) == 0:
                    if instrument == 'ULTRASPEC':
                        brow += [0] + 21*[None]
                    elif instrument == 'ULTRACAM':
                        brow += [0] + 24*[None]
                    else:
                        raise NotImplementedError('HiPERCAM case not done yet')
                else:
                    if instrument == 'ULTRASPEC':
                        min_med = np.min(medians[cnam])
                        med_med = np.median(medians[cnam])
                        max_med = np.max(medians[cnam])
                        brow += [len(means[cnam]),min_med,med_med,max_med]
                    elif instrument == 'ULTRACAM':
                        min_medl = np.min(medians[cnam]['L'])
                        med_medl = np.median(medians[cnam]['L'])
                        max_medl = np.max(medians[cnam]['L'])
                        min_medr = np.min(medians[cnam]['R'])
                        med_medr = np.median(medians[cnam]['R'])
                        max_medr = np.max(medians[cnam]['R'])
                        brow += [len(means[cnam]),min_medl,med_medl,max_medl,min_medr,med_medr,max_medr]
                    else:
                        raise NotImplementedError('HiPERCAM case not done yet')

                    brow += [
                        np.min(means[cnam]),np.median(means[cnam]),np.max(means[cnam])
                    ]
                    brow += [
                        np.min(p1s[cnam]),np.median(p1s[cnam]),np.max(p1s[cnam])
                    ]
                    brow += [
                        np.min(p16s[cnam]),np.median(p16s[cnam]),np.max(p16s[cnam])
                    ]
                    brow += [
                        np.min(rmsps[cnam]),np.median(rmsps[cnam]),np.max(rmsps[cnam])
                    ]
                    brow += [
                        np.min(p84s[cnam]),np.median(p84s[cnam]),np.max(p84s[cnam])
                    ]
                    brow += [
                        np.min(p99s[cnam]),np.median(p99s[cnam]),np.max(p99s[cnam])
                    ]

            barr.append(brow)

        # Create pandas dataframe for easy output
        if instrument == 'ULTRASPEC':
            colnames = ULTRASPEC_META_COLNAMES
        elif instrument == 'ULTRACAM':
            colnames = ULTRACAM_META_COLNAMES
        else:
            raise NotImplementedError('HiPERCAM case not done yet')

        cnames, dtypes = [], {}
        for cname, dtype, defn in colnames:
            cnames.append(cname)
            dtypes[cname] = dtype
        table = pd.DataFrame(data=barr,columns=cnames)
        table = table.astype(dtypes)
        table.to_csv(stats,index=False)
        print(f'Written statistics to {stats}\n')

# Create and write out spreadsheet
ULTRASPEC_META_COLNAMES = (
    ('run_no', 'str', 'Run number'),
    ('nf_used', 'float32', 'number of frames used for stats. First of many stats computed by "hmeta"'),
    ('med_min', 'float32', 'minimum median'),
    ('med_med', 'float32', 'median median'),
    ('med_max', 'float32', 'maximum median'),
    ('mean_min', 'float32', 'minimum mean'),
    ('mean_med', 'float32', 'median mean'),
    ('mean_max', 'float32', 'maximum mean'),
    ('p1_min', 'float32', 'minimum 1%-ile'),
    ('p1_med', 'float32', 'median 1%-ile'),
    ('p1_max', 'float32', 'maximum 1%-ile'),
    ('p16_min', 'float32', 'minimum 15.865%-ile'),
    ('p16_med', 'float32', 'median 15.865%-ile'),
    ('p16_max', 'float32', 'maximum 15.865%-ile'),
    ('rmsp_min', 'float32', 'minimum RMS estimated from (84-16%-ile)/2'),
    ('rmsp_med', 'float32', 'median  RMS estimated from (84-16%-ile)/2'),
    ('rmsp_max', 'float32', 'maximum  RMS estimated from (84-16%-ile)/2'),
    ('p84_min', 'float32', 'minimum 84.135%-ile'),
    ('p84_med', 'float32', 'median 84.135%-ile'),
    ('p84_max', 'float32', 'maximum 84.135%-ile'),
    ('p99_min', 'float32', 'minimum 99%-ile'),
    ('p99_med', 'float32', 'median 99%-ile'),
    ('p99_max', 'float32', 'maximum 99%-ile'),
)

ULTRACAM_META_COLNAMES = (
    ('run_no', 'str', 'Run number'),

    ('nf_used_1', 'float32', 'number of frames used for stats, CCD 1. First of many stats computed by "hmeta"'),
    ('med_min_L1', 'float32', 'minimum median left, CCD 1'),
    ('med_med_L1', 'float32', 'median median left CCD 1'),
    ('med_max_L1', 'float32', 'maximum median left CCD 1'),
    ('med_min_R1', 'float32', 'minimum median right, CCD 1'),
    ('med_med_R1', 'float32', 'median median right CCD 1'),
    ('med_max_R1', 'float32', 'maximum median right CCD 1'),
    ('mean_min_1', 'float32', 'minimum mean CCD 1'),
    ('mean_med_1', 'float32', 'median mean CCD 1'),
    ('mean_max_1', 'float32', 'maximum mean CCD 1'),
    ('p1_min_1', 'float32', 'minimum 1%-ile CCD 1'),
    ('p1_med_1', 'float32', 'median 1%-ile CCD 1'),
    ('p1_max_1', 'float32', 'maximum 1%-ile CCD 1'),
    ('p16_min_1', 'float32', 'minimum 15.865%-ile CCD 1'),
    ('p16_med_1', 'float32', 'median 15.865%-ile CCD 1'),
    ('p16_max_1', 'float32', 'maximum 15.865%-ile CCD 1'),
    ('rmsp_min_1', 'float32', 'minimum RMS estimated from (84-16%-ile)/2 CCD 1'),
    ('rmsp_med_1', 'float32', 'median  RMS estimated from (84-16%-ile)/2 CCD 1'),
    ('rmsp_max_1', 'float32', 'maximum  RMS estimated from (84-16%-ile)/2 CCD 1'),
    ('p84_min_1', 'float32', 'minimum 84.135%-ile CCD 1'),
    ('p84_med_1', 'float32', 'median 84.135%-ile CCD 1'),
    ('p84_max_1', 'float32', 'maximum 84.135%-ile CCD 1'),
    ('p99_min_1', 'float32', 'minimum 99%-ile CCD 1'),
    ('p99_med_1', 'float32', 'median 99%-ile CCD 1'),
    ('p99_max_1', 'float32', 'maximum 99%-ile CCD 1'),

    ('nf_used_2', 'float32', 'number of frames used for stats, CCD 2'),
    ('med_min_L2', 'float32', 'minimum median left, CCD 2'),
    ('med_med_L2', 'float32', 'median median left CCD 2'),
    ('med_max_L2', 'float32', 'maximum median left CCD 2'),
    ('med_min_R2', 'float32', 'minimum median right, CCD 2'),
    ('med_med_R2', 'float32', 'median median right CCD 2'),
    ('med_max_R2', 'float32', 'maximum median right CCD 2'),
    ('mean_min_2', 'float32', 'minimum mean CCD 2'),
    ('mean_med_2', 'float32', 'median mean CCD 2'),
    ('mean_max_2', 'float32', 'maximum mean CCD 2'),
    ('p1_min_2', 'float32', 'minimum 1%-ile CCD 2'),
    ('p1_med_2', 'float32', 'median 1%-ile CCD 2'),
    ('p1_max_2', 'float32', 'maximum 1%-ile CCD 2'),
    ('p16_min_2', 'float32', 'minimum 15.865%-ile CCD 2'),
    ('p16_med_2', 'float32', 'median 15.865%-ile CCD 2'),
    ('p16_max_2', 'float32', 'maximum 15.865%-ile CCD 2'),
    ('rmsp_min_2', 'float32', 'minimum RMS estimated from (84-16%-ile)/2 CCD 2'),
    ('rmsp_med_2', 'float32', 'median  RMS estimated from (84-16%-ile)/2 CCD 2'),
    ('rmsp_max_2', 'float32', 'maximum  RMS estimated from (84-16%-ile)/2 CCD 2'),
    ('p84_min_2', 'float32', 'minimum 84.135%-ile CCD 2'),
    ('p84_med_2', 'float32', 'median 84.135%-ile CCD 2'),
    ('p84_max_2', 'float32', 'maximum 84.135%-ile CCD 2'),
    ('p99_min_2', 'float32', 'minimum 99%-ile CCD 2'),
    ('p99_med_2', 'float32', 'median 99%-ile CCD 2'),
    ('p99_max_2', 'float32', 'maximum 99%-ile CCD 2'),

    ('nf_used_3', 'float32', 'number of frames used for stats, CCD 3'),
    ('med_min_L3', 'float32', 'minimum median left, CCD 3'),
    ('med_med_L3', 'float32', 'median median left CCD 3'),
    ('med_max_L3', 'float32', 'maximum median left CCD 3'),
    ('med_min_R3', 'float32', 'minimum median right, CCD 3'),
    ('med_med_R3', 'float32', 'median median right CCD 3'),
    ('med_max_R3', 'float32', 'maximum median right CCD 3'),
    ('mean_min_3', 'float32', 'minimum mean CCD 3'),
    ('mean_med_3', 'float32', 'median mean CCD 3'),
    ('mean_max_3', 'float32', 'maximum mean CCD 3'),
    ('p1_min_3', 'float32', 'minimum 1%-ile CCD 3'),
    ('p1_med_3', 'float32', 'median 1%-ile CCD 3'),
    ('p1_max_3', 'float32', 'maximum 1%-ile CCD 3'),
    ('p16_min_3', 'float32', 'minimum 15.865%-ile CCD 3'),
    ('p16_med_3', 'float32', 'median 15.865%-ile CCD 3'),
    ('p16_max_3', 'float32', 'maximum 15.865%-ile CCD 3'),
    ('rmsp_min_3', 'float32', 'minimum RMS estimated from (84-16%-ile)/2 CCD 3'),
    ('rmsp_med_3', 'float32', 'median  RMS estimated from (84-16%-ile)/2 CCD 3'),
    ('rmsp_max_3', 'float32', 'maximum  RMS estimated from (84-16%-ile)/2 CCD 3'),
    ('p84_min_3', 'float32', 'minimum 84.135%-ile CCD 3'),
    ('p84_med_3', 'float32', 'median 84.135%-ile CCD 3'),
    ('p84_max_3', 'float32', 'maximum 84.135%-ile CCD 3'),
    ('p99_min_3', 'float32', 'minimum 99%-ile CCD 3'),
    ('p99_med_3', 'float32', 'median 99%-ile CCD 3'),
    ('p99_max_3', 'float32', 'maximum 99%-ile CCD 3'),
)
