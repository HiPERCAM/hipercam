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
    ULTRASPEC. It computes data on every runs it can find listing the
    mininum, maximum, mean and median and the following percentiles:
    1.0,5.0,15.865,84.135,95.0,99.0. It does this for all frames in a
    run up to a maximum of 100-200 frames. The data are listed along
    with the frame number separately for each CCD since this only looks
    at valid data, i.e. it takes into account nblue for ULTRACAM and nskips
    for |hiper|. This means there may be different numbers of frames listed
    for each CCD. It writes a file runXXX[X].csv for each run containing data
    to a sub-directory called meta.

    The program only considers genuine frames, i.e. it copes with the
    nblue and nskip options of ULTRACAM and |hipercam|. It also knows
    about the different amplifier outputs of the instruments and
    starts by subtracting the medians of each amplifier output to
    avoid being disturbed by variable mean bias offsets. These offsets
    are recorded as well.

    hmeta takes a good while to run, so be nice when running it.

    """

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-n",
        dest="night",
        help="name of specific night [for testing purposes]",
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
        print(
            "** hmeta: cannot find either ultracam, "
            "ultraspec or hipercam in path"
        )
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
        and os.path.isdir(nname) and
        (args.night is None or nname == args.night)
    ]
    nnames.sort()

    if len(nnames) == 0:
        print("no night directories found", file=sys.stderr)
        print("hmeta aborted", file=sys.stderr)
        return

    # This defines approx how many we will analyse; could be up to
    # double this in practice
    NLIM = 100

    for nname in nnames:

        print(f"Night {nname}")

        # load all the run names
        if itype == 'U':
            runs = [
                run[:-4] for run in os.listdir(nname) if ure.match(run) and
                os.path.exists(os.path.join(nname,run[:-4]+'.dat'))
            ]
        else:
            runs = [run[:-5] for run in os.listdir(nname) if hre.match(run)]
        runs.sort()

        if len(runs) == 0:
            print(f' No runs with data found in {nname}; skipping')
            continue

        # create directory for any meta info such as the times
        meta = os.path.join(nname, 'meta')
        os.makedirs(meta, exist_ok=True)

        # Accumulate results in an array
        for run in runs:
            dfile = os.path.join(nname, run)

            # check for existence of data files
            for cnam in cnams:
                oname = os.path.join(meta, f'{run}_{cnam}.csv')
                if not os.path.exists(oname):
                    break
            else:
                if len(cnams) > 1:
                    print(f'  {dfile}: stats files already exist and will not be re-created')
                else:
                    print(f'  {dfile}: stats file already exists and will not be re-created')
                continue

            try:
                if itype == 'U':
                    rdat = hcam.ucam.Rdata(dfile)
                else:
                    rdat = hcam.hcam.Rdata(dfile, 1, False, False)
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

            # For speed, analyse a maximum of NLIM to 2*NLIM-1
            # images of each CCD from any given run. Have to take
            # into account the skips / nblue parameters.
            ntotal = rdat.ntotal()
            if ntotal == 0:
                print(f'  {dfile} -- zero frames; skipping')
                continue

            ncframes = {}
            if linstrument == 'ultraspec':
                # just the one CCD here. nstep defines how often we sample.
                # note that ntotal has to be 2*NLIM before nstep > 1.
                nstep = ntotal // min(ntotal, NLIM)
                ncframes['1'] = list(range(1,ntotal+1,nstep))
                nframes = ncframes['1']

            elif linstrument == 'ultracam':
                # CCD 1, 2 read out each time, but 3 can be skipped
                nstep = ntotal // min(ntotal, NLIM)
                ncframes['1'] = list(range(1,ntotal+1,nstep))
                ncframes['2'] = ncframes['1']

                nb = rdat.nblue
                nbstep = nb*max(1, ntotal // min(ntotal, NLIM*nb))
                ncframes['3'] = list(range(nb,ntotal+1,nbstep))
                nframes = sorted(set(ncframes['1']+ncframes['3']))

            else:
                # All CCDs potentially skippable
                nskips = rdat.nskips
                nsteps = []
                for cnam, nskip in zip( cnams, nskips):
                    nstep = (nskip+1)*max(
                        1, ntotal // min(
                            ntotal, NLIM*(nskip+1)
                        )
                    )
                    ncframes[cnam] = list(range(nskip+1,ntotal+1,nstep))

                # sorted list of frames to access
                nframes = sorted(
                    set(
                        ncframes['1']+ncframes['2']+ncframes['3']+ncframes['4']+ncframes['5']
                    )
                )

            # index counters
            ns = dict(zip(cnams,len(cnams)*[0]))

            # define arrays for holding the stats
            medians, means, mins, p1s, p5s, p16s, p84s, p95s, p99s, maxs = \
                {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
            for cnam, ncframe in ncframes.items():
                if linstrument == 'ultraspec':
                    medians[cnam] = np.empty_like(ncframe,dtype=float)
                elif linstrument == 'ultracam':
                    medians[cnam] = {
                        'L' : np.empty_like(ncframe,dtype=float),
                        'R' : np.empty_like(ncframe,dtype=float)
                    }
                elif linstrument == 'hipercam':
                    medians[cnam] = {
                        'E' : np.empty_like(ncframe,dtype=float),
                        'F' : np.empty_like(ncframe,dtype=float),
                        'G' : np.empty_like(ncframe,dtype=float),
                        'H' : np.empty_like(ncframe,dtype=float)
                    }

                means[cnam] = np.empty_like(ncframe,dtype=float)
                mins[cnam] = np.empty_like(ncframe,dtype=float)
                p1s[cnam] = np.empty_like(ncframe,dtype=float)
                p5s[cnam] = np.empty_like(ncframe,dtype=float)
                p16s[cnam] = np.empty_like(ncframe,dtype=float)
                p84s[cnam] = np.empty_like(ncframe,dtype=float)
                p95s[cnam] = np.empty_like(ncframe,dtype=float)
                p99s[cnam] = np.empty_like(ncframe,dtype=float)
                maxs[cnam] = np.empty_like(ncframe,dtype=float)

            # now access the data and calculate stats. we have to
            # remember that ultracam and hipercam have different
            # readout amps so we subtract median values calculated for
            # each amp separately, which is a little painful.
            for n, nf in enumerate(nframes):
                try:
                    mccd = rdat(nf)
                    for cnam, ncframe in ncframes.items():
                        if nf in ncframe:
                            ccd = mccd[cnam]
                            nc = ns[cnam]
                            if linstrument == 'ultraspec':
                                medval = ccd.median()
                                ccd -= medval
                                medians[cnam][nc] = medval

                            elif linstrument == 'ultracam':
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

                            elif linstrument == 'hipercam':
                                we = hcam.Group(hcam.Window)
                                wf = hcam.Group(hcam.Window)
                                wg = hcam.Group(hcam.Window)
                                wh = hcam.Group(hcam.Window)
                                for nw, wnam in enumerate(ccd):
                                    if wnam.startswith('E'):
                                        we[wnam] = ccd[wnam]
                                    elif wnam.startswith('F'):
                                        wf[wnam] = ccd[wnam]
                                    elif wnam.startswith('G'):
                                        wg[wnam] = ccd[wnam]
                                    elif wnam.startswith('H'):
                                        wh[wnam] = ccd[wnam]

                                ccde = hcam.CCD(we,ccd.nxtot,ccd.nytot)
                                mede = ccde.median()
                                ccde -= mede

                                ccdf = hcam.CCD(wf,ccd.nxtot,ccd.nytot)
                                medf = ccdf.median()
                                ccdf -= medf

                                ccdg = hcam.CCD(wg,ccd.nxtot,ccd.nytot)
                                medg = ccdg.median()
                                ccdg -= medg

                                ccdh = hcam.CCD(wh,ccd.nxtot,ccd.nytot)
                                medh = ccdh.median()
                                ccdh -= medh

                                medians[cnam]['E'][nc] = mede
                                medians[cnam]['F'][nc] = medf
                                medians[cnam]['G'][nc] = medg
                                medians[cnam]['H'][nc] = medh

                            # At this stage the median value should
                            # have been subtracted from the CCD on a
                            # per output basis. Remaining stats
                            # calculated from these median subtracted
                            # images.
                            means[cnam][nc] = ccd.mean()
                            mins[cnam][nc] = ccd.min()
                            p1,p5,p16,p84,p95,p99 = ccd.percentile([1.0,5.0,15.865,84.135,95.0,99.0])
                            p1s[cnam][nc] = p1
                            p5s[cnam][nc] = p5
                            p16s[cnam][nc] = p16
                            p84s[cnam][nc] = p84
                            p95s[cnam][nc] = p95
                            p99s[cnam][nc] = p99
                            maxs[cnam][nc] = ccd.max()

                            ns[cnam] += 1
                except:
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_tb(
                        exc_traceback, limit=1, file=sys.stderr
                    )
                    traceback.print_exc(file=sys.stderr)
                    continue

            # All numbers extracted from the run. Now want to create
            # pandas dataframe for easy output

            # Create pandas dataframe for easy output. First names and
            # datatypes
            cnames, dtypes = get_cnames_dtypes(linstrument)

            for cnam in cnams:
                # one file per CCD because of potentially different
                # column lengths
                oname = os.path.join(meta, f'{run}_{cnam}.csv')

                # Form output array
                barr = [ncframes[cnam],]
                if linstrument == 'ultraspec':
                    barr.append(medians[cnam])
                elif linstrument == 'ultracam':
                    barr += [medians[cnam]['L'],medians[cnam]['R']]
                elif linstrument == 'hipercam':
                    barr += [
                        medians[cnam]['E'],medians[cnam]['F'],
                        medians[cnam]['G'],medians[cnam]['H']
                    ]

                barr += [
                    means[cnam], mins[cnam],
                    p1s[cnam], p5s[cnam], p16s[cnam],
                    p84s[cnam], p95s[cnam], p99s[cnam],
                    maxs[cnam]
                ]

                # Make pandas dataframe
                table = pd.DataFrame(data=np.column_stack(barr),columns=cnames)
                table = table.astype(dtypes)
                table.to_csv(oname,index=False)
                print(f'  {dfile}, written stats to {oname}')

def get_cnames_dtypes(instrument):
    if instrument == 'ultraspec':
        coldefs = ULTRASPEC_COLNAMES
    elif instrument == 'ultracam':
        coldefs = ULTRACAM_COLNAMES
    elif instrument == 'hipercam':
        coldefs = HIPERCAM_COLNAMES

    cnames, dtypes = [], {}
    for cname, dtype, defn in coldefs:
        cnames.append(cname)
        dtypes[cname] = dtype

    return (cnames, dtypes)

ULTRASPEC_COLNAMES = (
    ('frame_no', 'float32', 'Frame number'),
    ('median', 'float32', 'median offset'),
    ('mean', 'float32', 'mean'),
    ('min', 'float32', 'minimum'),
    ('p1', 'float32', '1%-ile'),
    ('p5', 'float32', '5%-ile'),
    ('p16', 'float32', '15.865%-ile'),
    ('p84', 'float32', '84.135%-ile'),
    ('p95', 'float32', '95%-ile'),
    ('p99', 'float32', '99%-ile'),
    ('max', 'float32', 'maximum')
)

ULTRACAM_COLNAMES = (
    ('frame_no', 'float32', 'Frame number'),
    ('median_L', 'float32', 'left window median offset'),
    ('median_R', 'float32', 'right window median offset'),
    ('mean', 'float32', 'mean'),
    ('min', 'float32', 'minimum'),
    ('p1', 'float32', '1%-ile'),
    ('p5', 'float32', '5%-ile'),
    ('p16', 'float32', '15.865%-ile'),
    ('p84', 'float32', '84.135%-ile'),
    ('p95', 'float32', '95%-ile'),
    ('p99', 'float32', '99%-ile'),
    ('max', 'float32', 'maximum')
)

HIPERCAM_COLNAMES = (
    ('frame_no', 'float32', 'Frame number'),
    ('median_E', 'float32', 'E-quadrant median offset'),
    ('median_F', 'float32', 'F-quadrant median offset'),
    ('median_G', 'float32', 'G-quadrant median offset'),
    ('median_H', 'float32', 'H-quadrant median offset'),
    ('mean', 'float32', 'mean'),
    ('min', 'float32', 'minimum'),
    ('p1', 'float32', '1%-ile'),
    ('p5', 'float32', '5%-ile'),
    ('p16', 'float32', '15.865%-ile'),
    ('p84', 'float32', '84.135%-ile'),
    ('p95', 'float32', '95%-ile'),
    ('p99', 'float32', '99%-ile'),
    ('max', 'float32', 'maximum')
)
