"""Command line script to generate meta data plots"""

import sys
import os
import time
import re
import traceback

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
    """redplt

    This command is to be run in the "raw_data" directory containing
    night-by-night directories of data for |hipercam|, ULTRACAM or
    ULTRASPEC. It attempts to generate plots of any runs it finds with
    corresponding reduce logs files inside a sub-directory reduce and
    then stores these inside a sub-directory "meta". It only does this
    for runs which appear in a file called meta/times because only runs
    with >= 20 frames and lasting >= 600 seconds will be plotted. The
    purpose of these plots is so that they can be attached to the runs
    logs as a quick look on past runs.

    The code assumes that aperture 1 contains the target while aperture 2
    has the best comparison.

    """

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
        and os.path.exists(os.path.join(nname,'meta','times'))
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

        # Read the times
        times = os.path.join(mdir, 'times')
        tdata = {}
        with open(times) as tin:
            for line in tin:
                arr = line.split()
                tdata[arr[0]] = arr[1:]
        print('  Read timing data from',times)

        # Create plots, where possible.
        for run in runs:
            rlog = os.path.join(rdir, run + '.log')

            # runs a few checks
            if run not in tdata:
                print(f'  Run {run} has a log but no entry in meta/times; skipping')
                continue

            ut_start, mjd_start, ut_end, mjd_end, cadence, expose, nok, ntotal = tdata[run]
            if mjd_end == 'UNDEF':
                print(f'  Run {run} has a log but undefined times in meta/times; skipping')
                continue

            ttotal = int(86400*(float(mjd_end)-float(mjd_start)))
            nok = int(nok)
            if ttotal < 300 or nok < 10:
                print(f'  Run {run} lasts {ttotal} secs, and has {nok} frames; one or both is too small -- skipping')
                continue

            # OK attempt a plot
            try:
                pname = os.path.join(mdir,run + '.png')

                # Two panels, target / comparison and comparison
                fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)

                hlog = hcam.hlog.Hlog.read(rlog)
                cnams = sorted(list(hlog.keys()))
                apnames = hlog.apnames

                # use the next to work out optimal plot ranges
                tymin, tymax, cymin, cymax = 4*[None]
                for cnam in cnams:
                    if cnam in apnames:
                        apnams = apnames[cnam]
                        if '1' in apnams:
                            targ = hlog.tseries(cnam,'1')
                            if '2' in apnams:
                                comp = hlog.tseries(cnam,'2')
                                targ /= comp
                                targ.normalise()
                                comp.normalise()

                                ndat = len(targ)
                                if ndat > 1000:
                                    # stop plotting too many points to keep
                                    # size down
                                    binsize = ndat // 1000
                                    targ.bin(binsize)
                                    comp.bin(binsize)

                                (_d,_d),(tylo,_d),(_d,tyhi) = targ.percentile([1,99],  bitmask=hcam.BAD_TIME)
                                (_d,_d),(cylo,_d),(_d,cyhi) = comp.percentile([1,99],  bitmask=hcam.BAD_TIME)
                                if tymax is not None:
                                    off = tymax - tylo
                                    targ += off
                                    tymax += tyhi-tylo
                                else:
                                    tymin, tymax = tylo, tyhi
                                if cymax is not None:
                                    off = cymax - cylo
                                    comp += off
                                    cymax += cyhi-cylo
                                else:
                                    cymin, cymax = cylo, cyhi

                                targ.mplot(ax1,utils.rgb(cols[cnam]),ecolor='0.5', bitmask=hcam.BAD_TIME)
                                comp.mplot(ax2,utils.rgb(cols[cnam]),ecolor='0.5',  bitmask=hcam.BAD_TIME)

                            else:
                                ndat = len(targ)
                                if ndat > 1000:
                                    # stop plotting too many points to keep
                                    # size down
                                    binsize = ndat // 1000
                                    targ.bin(binsize)

                                (_d,_d),(tylo,_d),(_d,tyhi) = targ.percentile([1,99],  bitmask=hcam.BAD_TIME)
                                if tymax is not None:
                                    off = tymax - tylo
                                    targ += off
                                    tymax += tyhi-tylo
                                else:
                                    tymin, tymax = tylo, tyhi
                                targ.mplot(ax1,utils.rgb(cols[cnam]),ecolor='0.5',  bitmask=hcam.BAD_TIME)

                yrange = tymax-tymin
                ax1.set_ylim(tymin-yrange/10, tymax+yrange/10)
                if cymin is not None:
                    yrange = cymax-cymin
                    ax2.set_ylim(cymin-yrange/10, cymax+yrange/10)
                ax1.set_ylabel('Normalised target data')
                ax1.set_title(f'{nname}, {run}')
                ax2.set_ylabel('Normalised comparison data')
                ax2.set_xlabel('Time [MJD]')
                plt.savefig(pname)
                plt.close()
                print(f'Written {pname}')

            except:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_tb(
                    exc_traceback, limit=1, file=sys.stderr
                )
                traceback.print_exc(file=sys.stderr)
                print("Problem reading log for run =", run)
