"""Command line script to generate meta data plots"""

import sys
import os
import time
import re
import traceback
import argparse

import numpy as np
import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = [
    "redplt",
]

##########################################################
#
# redplt -- generate standardised plots of reduce log data
#
##########################################################


def redplt(args=None):
    description = \
    """redplt

    This command is to be run in the "raw_data" directory containing
    night-by-night directories of data for |hipercam|, ULTRACAM or
    ULTRASPEC. It attempts to generate plots of any runs it finds with
    corresponding reduce logs files inside a sub-directory reduce and
    then stores these inside a sub-directory "meta". The purpose of
    these plots is so that they can be attached to the runs logs as a
    quick look on past runs.

    The code assumes that aperture 1 contains the target while
    aperture 2 has the best comparison. It produces plots in which the
    LHS panels show the target divided by the comparison while the RHS shows
    the comparison only as a crude check on clouds.

    Runs with fewer than 20 points or lasting less than 10 minutes will
    not be plotted.

    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        dest="full",
        action="store_true",
        help="carry out full re-computation of plots",
    )
    args = parser.parse_args()

    cwd = os.getcwd()
    if os.path.basename(cwd) != "raw_data":
        print("redplt must be run in a directory called 'raw_data'")
        print("redplt aborted", file=sys.stderr)
        return

    if cwd.find("ultracam") > -1:
        instrument = "ULTRACAM"
        itype = 'U'
        source = 'ul'
        cnams = ('1','2','3')
        cols = {'1' : "red", '2' : "green", '3' : "blue"}
    elif cwd.find("ultraspec") > -1:
        instrument = "ULTRASPEC"
        itype = 'U'
        source = 'ul'
        cnams = ('1',)
        cols = {'1' : "blue"}
    elif cwd.find("hipercam") > -1:
        instrument = "HiPERCAM"
        itype = 'H'
        source = 'hl'
        cnams = ('1','2','3','4','5')
        cols = {'1' : "blue", '2' : "green", '3' : "orange", '4' : "red", '5' : "mud"}
    else:
        print("cannot find either ultracam, ultraspec or hipercam in path")
        print("redplt aborted", file=sys.stderr)
        return

    linstrument = instrument.lower()

    # Now the actual work.  Next are regular expressions to match run
    # directories, nights, and run files
    nre = re.compile("^\d\d\d\d-\d\d-\d\d$")
    lre = re.compile("^run\d\d\d\\d?.log$")

    # Get list of night directories
    nnames = [
        nname
        for nname in os.listdir(".")
        if nre.match(nname)
        and os.path.isdir(nname)
        and os.path.exists(os.path.join(nname,'reduce'))
    ]
    nnames.sort()

    if len(nnames) == 0:
        print("no night directories found", file=sys.stderr)
        print("redplt aborted", file=sys.stderr)
        return

    for nname in nnames:

        print(f"Night {nname}")

        # reduce and meta directories
        rdir = os.path.join(nname,'reduce')
        mdir = os.path.join(nname,'meta')

        # load all the run names that can be found in reduce
        runs = [run[:-4] for run in os.listdir(rdir) if lre.match(run)]
        runs.sort()

        if len(runs) == 0:
            print(f' No run logs found in {rdir}; skipping')
            continue

        # ensure meta directory exists
        os.makedirs(mdir, exist_ok=True)

        # Minimum number of points / minutes to bother with
        NMIN, TMIN = 20, 10

        # Create plots, where possible.
        for run in runs:
            rlog = os.path.join(rdir, run + '.log')

            pname = os.path.join(mdir,run + '.png')
            if not args.full and os.path.exists(pname):
                print(f'  Plot {pname} exists and will not be re-computed')
                continue

            # OK attempt a plot
            try:

                hlog = hcam.hlog.Hlog.rascii(rlog)

                # LHS for targ/comp, RHS for comp
                cnams = sorted(list(hlog.keys()))
                fig, axs = plt.subplots(
                    len(cnams),2,sharex=True,
                    squeeze=False,
                    figsize=(16,len(cnams)*4)
                )

                made_a_plot = False
                apnames = hlog.apnames

                for nc, cnam in enumerate(cnams):
                    if cnam in apnames:
                        apnams = apnames[cnam]

                        if '1' in apnams:
                            # we have a target
                            targ = hlog.tseries(cnam,'1')

                            axl, axr = axs[nc,0], axs[nc,1]

                            if '2' in apnams:
                                # we have a comparison
                                comp = hlog.tseries(cnam,'2')

                                # run checks
                                ts = targ.t[~targ.get_mask(hcam.BAD_TIME) & ~comp.get_mask(hcam.BAD_TIME)]
                                if len(ts) < NMIN:
                                    print(f'{run}, CCD={cnam} has too few points ({len(ts)} < {NMIN})')
                                    continue
                                tmins = 1440*(ts.max()-ts.min())
                                if tmins < TMIN:
                                    print(f'{run}, CCD={cnam} is too short ({tmins} < {TMIN} mins)')
                                    continue

                                targ /= comp

                                ndat = len(targ)
                                if ndat > 3000:
                                    # avoid plotting too many points to keep
                                    # size down
                                    binsize = ndat // 1500
                                    targ.bin(binsize)
                                    comp.bin(binsize)

                                (_d,_d),(ylo,_d),(_d,yhi) = targ.percentile([2,98],  bitmask=hcam.BAD_TIME)
                                yrange = yhi-ylo
                                ylo -= 0.2*yrange
                                yhi += 0.2*yrange
                                targ.mplot(axl,utils.rgb(cols[cnam]),ecolor='0.5', bitmask=hcam.BAD_TIME)
                                axl.set_ylim(ylo,yhi)

                                (_d,_d),(ylo,_d),(_d,yhi) = comp.percentile([2,98],  bitmask=hcam.BAD_TIME)
                                yrange = yhi-ylo
                                ylo -= 0.2*yrange
                                yhi += 0.2*yrange
                                comp.mplot(axr,utils.rgb(cols[cnam]),ecolor='0.5', bitmask=hcam.BAD_TIME)
                                axr.set_ylim(ylo,yhi)

                                made_a_plot = True

                            else:
                                # target only

                                # run some checks
                                ts = targ.t[~targ.get_mask(hcam.BAD_TIME)]
                                if len(ts) < NMIN:
                                    print(f'{run}, CCD={cnam} has too few points ({len(ts)} < {NMIN})')
                                    continue
                                tmins = 1440*(ts.max()-ts.min())
                                if tmins < TMIN:
                                    print(f'{run}, CCD={cnam} is too short ({tmins} < {TMIN} mins)')
                                    continue

                                ndat = len(targ)
                                if ndat > 3000:
                                    # stop plotting too many points to keep
                                    # size down
                                    binsize = ndat // 1500
                                    targ.bin(binsize)

                                (_d,_d),(ylo,_d),(_d,yhi) = targ.percentile([2,98],  bitmask=hcam.BAD_TIME)
                                yrange = yhi-ylo
                                ylo -= 0.2*yrange
                                yhi += 0.2*yrange
                                targ.mplot(axl,utils.rgb(cols[cnam]),ecolor='0.5', bitmask=hcam.BAD_TIME)
                                axl.set_ylim(ylo,yhi)

                                made_a_plot = True

                            axl.set_ylabel('Targ / Comp')
                            axr.set_ylabel('Comp')

                if made_a_plot:
                    plt.tight_layout()
                    axs[0,0].set_title(f'{nname}, {run}')
                    axs[1,1].set_title(f'{nname}, {run}')
                    axs[-1,0].set_xlabel('Time [MJD]')
                    axs[-1,1].set_xlabel('Time [MJD]')
                    plt.savefig(pname)
                    print(f'Written {pname}')
                plt.close()

            except:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_tb(
                    exc_traceback, limit=1, file=sys.stderr
                )
                traceback.print_exc(file=sys.stderr)
                print("Problem reading log for run =", run)
