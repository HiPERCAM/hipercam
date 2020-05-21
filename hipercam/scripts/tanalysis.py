import sys
import os
import shutil
import struct

import numpy as np
import matplotlib.pylab as plt
from scipy import stats

import hipercam as hcam
from hipercam import cline, utils, spooler, ucam
from hipercam.cline import Cline

__all__ = [
    "tanalysis",
]

##################################
#
# tanalysis -- timing fixer script
#
##################################


def tanalysis(args=None):
    """``tanalysis [source] run [mintim dcmax] details (plot [check (output)])``

    .. Warning::

       This script should only be run if you know what you are doing.

    Analyses timing data in |hiper|, ULTRACAM or ULTRASPEC data.  You
    must first have made copy of the timing bytes (see ``tbytes'' and
    ``atbytes''). It looks for a timing bytes file in a sub-directory
    of the directory the script is run from called "tbytes". If the
    file does not exist, the script will exit.

    On occasion there are problems with ULTRACAM, ULTRASPEC and, very
    rarely, |hiper| timestamps. The main purpose of this script is to
    provide an analysis of the timing to spot these is as clear a way
    as possible. Some of the problems are::

      1) ULTRASPEC null timestamps: ULTRASPEC runs occasionally feature
         weird null timestamps. They don't replace the proper ones but
         push them down the line so they turn up later.

      2) Extra OK-looking timestamps: on three occasions, genuine but
         spurious timestamps have been added to the FIFO buffer. Again
         these are extra, and push the real ones down the line.

    Parameters:

        source : string [hidden]
           Data source, two options:

              | 'hl' : local HiPERCAM FITS file
              | 'ul' : ULTRA(CAM|SPEC) server

        run : string
           run number to access, e.g. 'run0034'. This will also be
           used to generate the name for the timing bytes file
           (extension '.tbts'). If a file of this name already exists,
           the script will attempt to read and compare the bytes of
           the two files and report any changes.  The timing bytes
           file will be written to the present working directory, not
           necessarily the location of the data file.

        mintim : int [hidden]
           Minimum number of times to attempt to do anything
           with. This must be at least 4 so that there are 3+ time
           differences to try to get a median time, but in practice it
           is probably advisable to use a larger number still.

        dcmax : float [hidden]
           Maximum differential in terms of cycle number from exact integers
           to use to indicate "failed" times. Pre 2010-02 needs to be looser
           than post 2010-02. Number of order 0.002 is fine post 2010-02,
           whereas 0.02 is more like it pre 2010-02.

        details : bool
           Print out detailed information on problem times

        plot : bool
           True to make diagnostic plots if problem runs found

        check : bool [if plot, hidden]
           True simply to check a run for timing problems, not try to
           fix them. Defaults to True, and is also hidden by
           default. If you set check=False, the code also attempts to
           fix the times and generate a new data file. At the moment
           it is set to try to correct for spurious extra timestamps
           as happened for the following runs 2011-10-31 run022,
           2013-12-31 run031 and 2015-06-22 run024.

        output : str [if not check]
           Name of the modified data file. Will not overwrite any existing file.

    .. Note::

       The final frame on many runs can be terminated early and thus ends with
       an out of sequence timestamp (it comes earlier than expected). The script
       regards such a run as being "OK" if there are not other issues.

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("source", Cline.LOCAL, Cline.HIDE)
        cl.register("run", Cline.GLOBAL, Cline.PROMPT)
        cl.register("mintim", Cline.LOCAL, Cline.HIDE)
        cl.register("dcmax", Cline.LOCAL, Cline.HIDE)
        cl.register("plot", Cline.LOCAL, Cline.PROMPT)
        cl.register("details", Cline.LOCAL, Cline.PROMPT)
        cl.register("check", Cline.LOCAL, Cline.HIDE)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value(
            "source", "data source [hl, ul]", "hl", lvals=("hl", "ul")
        )

        run = cl.get_value("run", "run name", "run005")
        if run.endswith(".fits"):
            run = run[:-5]

        mintim = cl.get_value("mintim", "minimum number of times needed", 6, 4)
        dcmax = cl.get_value(
            "dcmax", "maximum cycle number offset from an integer", 0.003, 0
        )
        details = cl.get_value(
            "details", "print detailed diagnostics of problems", True
        )
        plot = cl.get_value("plot", "make diagnostic plots", True)

        if plot:
            cl.set_default("check", True)
            check = cl.get_value("check", "check for, but do not fix, any problems", True)
            if not check:
                output = cl.get_value(
                    "output", "output file for data with modified times",
                    cline.Fname("modified",ftype=cline.Fname.NOCLOBBER)
                )
        else:
            check = True

    # create name of timing file
    tfile = os.path.join("tbytes", os.path.basename(run) + hcam.TBTS)
    if not os.path.isfile(tfile):
        raise hcam.HipercamError("Could not find timing file = {}".format(tfile))

    # create name of run file and the copy which will only get made
    # later if problems are picked up
    if source == "hl":
        rfile = run + ".fits"
    else:
        rfile = run + ".dat"

    # first load all old time stamps, and compute MJDs
    if source == "ul":

        # Load up the time stamps from the timing data file. Need also
        # the header of the original run to know how many bytes to
        # read and how to interpret them.
        rhead = hcam.ucam.Rhead(run)
        with open(tfile, "rb") as fin:
            atbytes, mjds, tflags, gflags = [], [], [], []
            nframe = 0
            while 1:
                tbytes = fin.read(rhead.ntbytes)
                if len(tbytes) != rhead.ntbytes:
                    break
                nframe += 1
                atbytes.append(tbytes)

                # interpret times
                mjd, tflag, gflag, uformat = u_tbytes_to_mjd(tbytes, rhead, nframe)
                mjds.append(mjd)
                tflags.append(tflag)
                gflags.append(gflag)

    else:
        raise NotImplementedError("source = {} not yet implemented".format(source))


    if len(mjds) < mintim:
        # Must have specified minimum of times to work with. This is
        # checked more stringently later, but this avoids silly crashes
        # when there are no times.
        print(
            run,
            "has too few frames to work with ({} vs minimum = {})".format(
                len(mjds), mintim
            ),
        )
        return

    # Independent of source, at this stage 'atbytes' is a list of all
    # timestamp bytes, while 'mjds' is a list of all equivalent MJDs,
    # and 'tflags' are bools which are True for OK times, False for
    # null timestamps.
    mjds = np.array(mjds)
    tflags = np.array(tflags)
    gflags = np.array(gflags)
    inds = np.arange(len(mjds))
    nulls = ~tflags
    nulls_present = nulls.any()
    btimes = ~gflags
    bad_times_present = btimes.any()

    # Remove null timestamps
    mjds_ok = mjds[tflags & gflags]
    inds_ok = inds[tflags & gflags]

    # Must have specified minimum of times to work with.
    if len(mjds_ok) < mintim:
        print(
            run,
            "has too few non-null times to work with ({} vs minimum = {})".format(
                len(mjds_ok), mintim
            ),
        )
        return

    # Median time difference of GPS timestamps with guard against identical
    # times that can dominate leading to mdiff = 0
    diffs = mjds_ok[1:] - mjds_ok[:-1]
    pos = diffs > 1.0e-9
    if len(diffs[pos]) == 0:
        print(run, "has zero positive time differences to work with")
        return
    mdiff = np.median(diffs[pos])

    # Maximum deviation from median separation to allow (days)
    # 1.e-9 ~ 1e-4 seconds, ~30x shorter than shortest ULTRACAM
    # cycle time. This could allow some bad stamps through but
    # they will be spotted later by looking at cycle numbers.
    MDIFF = 1.0e-9

    # Identify OK timestamps. Could miss some at this point, but
    # that's OK) kick off by marking all timestamps as False
    ok = mjds_ok == mjds_ok + 1
    for n in range(len(mjds_ok) - 1):
        if abs(mjds_ok[n + 1] - mjds_ok[n] - mdiff) < MDIFF:
            # if two timestamps differ by the right amount, mark both
            # as ok
            ok[n] = ok[n + 1] = True

    gmjds = mjds_ok[ok]
    ginds = inds_ok[ok]
    if len(gmjds) < 2:
        print("{}: found fewer than 2 good timestamps".format(run))
        return

    # Work out integer cycle numbers. First OK timestamp is given
    # cycle number = 0 automatically by the method used.  The cycle
    # numbers can go wrong if the median is not precise enough leading
    # to jumps in the cycle number on runs with large numbers of short
    # exposures, so we build up to the full fit in stages, doubling the
    # number fitted each time

    NMAX = 100
    cycles = (gmjds[:NMAX] - gmjds[0]) / mdiff
    moffset = (cycles - np.round(cycles)).mean()
    icycles = np.round(cycles - moffset).astype(int)

    while NMAX < len(gmjds):
        # fit linear trend of first NMAX times where hopefully NMAX is small
        # enough for there not to be a big error but large enough to allow extrapolation.
        slope, intercept, r, p, err = stats.linregress(icycles[:NMAX], gmjds[:NMAX])
        NMAX *= 1.5
        NMAX = int(NMAX)
        icycles = np.round((gmjds[:NMAX] - intercept) / slope).astype(int)

    # hope the cycle numbers should be good across the board for all the "gmjds"
    # at this point. Carry out final fits to try to refine the intercept and
    # slope because it is possible to be upset by jumps in time of less
    # than a cycle (i.e 2002-09-20, run014)

    cok = icycles == icycles
    cinds = np.arange(len(cok))
    nrej = 0
    while 1:
        # fit to ok ones
        slope, intercept, r, p, err = stats.linregress(icycles[cok], gmjds[cok])
        cycles = (gmjds - intercept) / slope

        # cdiffs, acdiffs -- cycle differences and their absolute
        # values
        cdiffs = cycles - icycles
        acdiffs = np.abs(cdiffs)

        dmax = np.argmax(acdiffs[cok])
        if acdiffs[cok][dmax] > dcmax:
            # at least one is out of spec. reject it
            cok[cinds[cok][dmax]] = False
            nrej += 1
        else:
            break
    print("Rejected {} times when fitting nominally good times".format(nrej))

    # now try to compute cycle numbers for *all* non-null timestamps
    acycles = (mjds_ok - intercept) / slope
    iacycles = np.round(acycles).astype(int)

    if iacycles[-1] == iacycles[-2]:
        # because the last frame can be terminated early, it is
        # not uncommon for it to end with same cycle number as the
        # penultimate one. Correct this here
        iacycles[-1] += 1

    cadiffs = acycles - iacycles
    monotonic = (iacycles[1:] - iacycles[:-1] > 0).all()

    # check for early termination on run causing last frame to appear
    # early
    terminated_early = cadiffs[-1] < -dcmax

    # check that all is OK
    if (
        not nulls_present
        and monotonic
        and (
            (terminated_early and (np.abs(cadiffs[:-1]) < dcmax).all())
            or (not terminated_early and (np.abs(cadiffs) < dcmax).all())
        )
    ):
        if terminated_early:
            mdev = np.abs(cadiffs[:-1]).max()
            print(
                "{} times are OK; {} frames; max. dev. all bar last = {:.2g} cyc, {:.2g} msec; last = {:.2g} cyc".format(
                    run, len(iacycles), mdev, 86400 * 1000 * slope * mdev, cadiffs[-1]
                )
            )
        else:
            mdev = np.abs(cadiffs).max()
            print(
                "{} times are OK; {} frames; max. dev. all = {:.2g} cyc, {:.2g} msec".format(
                    run, len(iacycles), mdev, 86400 * 1000 * slope * mdev
                )
            )
        # to prevent crash during plot
        fails = cadiffs != cadiffs
        run_ok = True

    else:
        # search for duplicate cycle numbers
        u, c = np.unique(iacycles, return_counts=True)
        dupes = u[c > 1]
        dupes_present = len(dupes) > 0

        ntot = len(mjds)
        ndupe = len(dupes)
        nnull = len(mjds[nulls])
        nbad = len(mjds[btimes])
        fails = np.abs(cadiffs) > dcmax
        if terminated_early:
            nfail = len(cadiffs[fails]) - 1
            mdev = np.abs(cadiffs[:-1]).max()
        else:
            nfail = len(cadiffs[fails])
            mdev = np.abs(cadiffs).max()

        if ndupe == 0 and nbad == 0 and nnull == 0 and inds_ok[fails][-1] < nfail:
            # save some runs with trivial level issues from being flagged
            if nfail < len(cadiffs):
                mdev = np.abs(cadiffs[nfail:]).max()
            comment = " [All fails are at start]"
            run_ok = True
        else:
            comment = ""
            run_ok = False

        # Summarise problems
        print(
            "{} timestamps corrupt. TOT,DUP,NULL,BAD,FAIL,FOK = {},{},{},{},{},{}; max dev = {:.2g} cyc, {:.2g} sec.{}".format(
                run,
                ntot,
                ndupe,
                nnull,
                nbad,
                nfail,
                inds_ok[ok][0] + 1,
                mdev,
                86400 * slope * mdev,
                comment,
            )
        )

        if details:

            # Some details
            fcadiffs = iacycles[1:] - iacycles[:-1]
            back = fcadiffs <= 0
            oinds = inds_ok[:-1][back]
            ninds = np.arange(len(fcadiffs))[back]
            cycs = iacycles[:-1][back]
            ncycs = iacycles[1:][back]
            tims = mjds_ok[:-1][back]
            if len(oinds):
                print("\nInverted order times:")
                for oind, nind, cyc, ncyc, mjd in zip(oinds, ninds, cycs, ncycs, tims):
                    print(
                        "  Old index = {}, new index = {}, cycle = {}, next cycle = {}, time = {}".format(
                            oind, nind, cyc, ncyc, mjd
                        )
                    )

            if nfail:
                print("\nFailed times:")
                inds = inds_ok[fails]
                cycs = iacycles[fails]
                cydiffs = cadiffs[fails]
                tims = mjds_ok[fails]
                for ind, cyc, cdiff, mjd in zip(inds, cycs, cydiffs, tims):
                    print(
                        "  Index = {}, cycle = {}, cycle diff. = {:.2g}, time = {}".format(
                            ind, cyc, cdiff, mjd
                        )
                    )

    if plot:
        # diagnostic plot: cycle differences vs cycle numbers
        def c2t(x):
            return 86400 * slope * x

        def t2c(x):
            return x / slope / 86400

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(iacycles[~fails], cadiffs[~fails], ".b")
        ax.plot(iacycles[fails], cadiffs[fails], ".r")
        ax.set_xlabel("Cycle number")
        ax.set_ylabel("Cycle difference")
        secxax = ax.secondary_xaxis("top", functions=(c2t, t2c))
        secxax.set_xlabel("Time [MJD - {}] (seconds)".format(intercept))
        secyax = ax.secondary_yaxis(
            "right", functions=(lambda x: 1000 * c2t(x), lambda x: t2c(x) / 1000)
        )
        secyax.set_xlabel("$\Delta t$ (msec)".format(intercept))
        plt.show()

    if not plot or run_ok or check:
        # go no further if run is OK, no plot has been made, or we are in check mode
        return

    # very specific code at this point for the post-2010 problem runs: remove duplicate
    # time stamps.

    if ginds[0] != 0 or nfail == 0  or (source == 'ul' and uformat != 2):
        # this case needs checking
        raise hcam.HipercamError(
            'Cannot yet handle cases where the first timestamp is no good, nfail == 0, or format==2 ultracam data'
        )

    # these frames have been ID-ed as junk
    skip_inds = inds_ok[fails]
    inds_keep = inds == inds
    inds_keep[skip_inds] = False

    # inds_keep : False for any Frame we want to forget
    new_mjds, new_atbytes = [], []
    for keep, mjd, tbytes in zip(inds_keep, mjds, atbytes):
        if keep:
            new_mjds.append(mjd)
            new_atbytes.append(tbytes)

    nskip = len(mjds)-len(new_mjds)

    # modify the last few times to make up for the skips
    for tbytes in atbytes[-nskip:]:
        mjd += slope
        new_mjds.append(mjd)
        unix_sec = ucam.DSEC*(mjd - ucam.UNIX)
        nsec = int(np.floor(unix_sec))
        nnsec = int(round(1.e7*(unix_sec-nsec)))
        tbytes = tbytes[:12] + struct.pack("<II", nsec, nnsec) + tbytes[20:]
        new_atbytes.append(tbytes)

    new_mjds_ok = np.array(new_mjds)
    new_mjds_ok = new_mjds_ok[tflags & gflags]
    acycles = (new_mjds_ok - intercept) / slope
    iacycles = np.round(acycles).astype(int)
    diffs = acycles-iacycles

    # search for duplicates
    u, c = np.unique(iacycles, return_counts=True)
    ndupes = len(u[c > 1])

    # search for gaps
    csteps = iacycles[1:]-iacycles[:-1]
    ngaps = len(iacycles[1:][csteps > 1])

    print('\nFound {} spurious time stamps; new times have {} duplications and {} gaps > 1'.format(nskip,ndupes,ngaps))
    print('Ideally, the last two numbers should be 0 and 0.')

    plt.plot(iacycles[:-nskip], diffs[:-nskip], '.g')
    plt.plot(iacycles[-nskip:], diffs[-nskip:], '.b')
    plt.xlabel('Cycle number')
    plt.ylabel('Cycle difference')
    plt.show()

    reply = input('Generate data with corrected time stamps? [no] ')
    if reply != 'yes':
        print('time correction aborted')

    # OK, finally, let's go for it!
    if source == "ul":
        with open(output,"wb") as fout:
            with open(rfile,"rb") as fin:
                nframe = 0
                while 1:
                    # Read a whole frame
                    buff = fin.read(rhead.framesize)
                    if len(buff) != rhead.framesize:
                        break

                    # write to output, writing modified timing
                    # bytes in place of old ones
                    fout.write(new_atbytes[nframe])
                    fout.write(buff[rhead.ntbytes:])
                    nframe += 1
                    print('Processed frame {} / {}'.format(nframe,len(atbytes)))

    else:
        raise NotImplementedError("source = {} not yet implemented".format(source))




def u_tbytes_to_mjd(tbytes, rtbytes, nframe):
    """Translates set of ULTRACAM timing bytes into an MJD.

    Marks ULTRASPEC null stamps as bad (32 bytes in length,
    the last 20 of which are 0)

    Returns (MJD, NotNull, GoodTime, format) where MJD is the MJD, NotNull is
    a flag to say whether the time stamp was not null (True=OK),
    GoodTime says whether the time is otherwise thought to be OK, and format
    is an integer code to distinguish different versions of the ultracam times.
    """
    try:
        ret = hcam.ucam.utimer(tbytes, rtbytes, nframe)
    except ValueError:
        return (0, True, False)

    if len(tbytes) == 32 and tbytes[12:] == 20 * b"\x00":
        return (ret[1]["gps"], False, False, ret[1]["format"])
    elif ret[1]["gps"] < 52395:
        return (ret[1]["gps"], True, False, ret[1]["format"])
    else:
        return (ret[1]["gps"], True, True, ret[1]["format"])


def h_tbytes_to_mjd(tbytes, nframe):
    """Translates set of HiPERCAM timing bytes into an MJD"""

    # number of seconds in a day
    DAYSEC = 86400.0

    (
        frameCount,
        timeStampCount,
        years,
        day_of_year,
        hours,
        mins,
        seconds,
        nanoseconds,
        nsats,
        synced,
    ) = htimer(tbytes)
    frameCount += 1

    if frameCount != nframe:
        if frameCount == nframe + 1:
            warnings.warn("frame count mis-match; a frame seems to have been dropped")
        else:
            warnings.warn(
                "frame count mis-match; {:d} frames seems to have been dropped".format(
                    frameCount - self.nframe
                )
            )

    try:
        imjd = gregorian_to_mjd(years, 1, 1) + day_of_year - 1
        fday = (hours + mins / 60 + (seconds + nanoseconds / 1e9) / 3600) / 24
    except ValueError:
        imjd = 51544
        fday = nframe / DAYSEC

    return imjd + fday
