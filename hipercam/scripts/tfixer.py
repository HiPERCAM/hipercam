import sys
import os
import shutil

import numpy as np
import matplotlib.pylab as plt
from scipy import stats

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ = ['tfixer',]

###############################
#
# tfixer -- timing fixer script
#
###############################

def tfixer(args=None):
    """``tfixer [source] run``

    .. Warning::

       This script should only be run if you know what you are doing.

    Fixes timestamps in |hiper|, ULTRACAM or ULTRASPEC data. This is
    carried out by first copying the data for safety and only
    optionally deleting the copy once the fixes have been made. It
    also requires that a copy of the timing bytes has been made
    previously (see ``tbytes'' and ``atbytes''). It looks for this
    file in a sub-directory of the directory the script is run from called
    "tbytes". If the file does not exist, the a script will exit.

    On occasion there are problems with ULTRACAM, ULTRASPEC and, very
    rarely, |hiper| timestamps. Most of these problems are
    fixable. The aim of this script is to try to accomplish such
    fixing is as bombproof a way as possible. These are the types of
    problems it tries to fix:

      1) ULTRASPEC null timestamps: ULTRASPEC runs occasionally feature
         weird null timestamps. They don't replace the proper ones but
         push them down the line so they turn up later.

      2) Extra OK-looking timestamps: we think that owing to glitches,
         sometimes genuine but spurious timestamps are added to the
         FIFO buffer. Again these are extra, and push the real ones down
         the line. They are however not as glaring as the null timestamps
         so a bit more care is needed to identify them.

    The program always defaults to safety: if there is anything that does
    not seem right, it will do nothing. If it does run, it will report either
    that the run times are "OK" or "corrupt".

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

        mintim : int
           Minimum number of times to attempt to do anything
           with. This must be at least 4 so that there are 3+ time
           differences to try to get a median time, but in practice it
           is probably advisable to use a larger number still.

        plot : bool
           True to make diagnostic plots if problem runs found

        check : bool
           True simply to check a run for timing problems, not try to
           fix them.

    .. Note::

       The final frame on many runs can be terminated early and thus ends with
       an out of sequence timestamp (it comes earlier than expected). The script
       regards such a run as being "OK".

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('source', Cline.LOCAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('mintim', Cline.LOCAL, Cline.PROMPT)
        cl.register('plot', Cline.LOCAL, Cline.PROMPT)
        cl.register('check', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value(
            'source', 'data source [hl, ul]',
            'hl', lvals=('hl','ul')
        )

        run = cl.get_value('run', 'run name', 'run005')
        if run.endswith('.fits'):
            run = run[:-5]

        mintim = cl.get_value('mintim', 'minimum number of times needed', 6, 4)
        plot = cl.get_value('plot', 'make diagnostic plots', True)
        check = cl.get_value('check', 'check for, but do not fix, any problems', True)

    # create name of timing file
    tfile = os.path.join('tbytes', os.path.basename(run) + hcam.TBTS)
    if not os.path.isfile(tfile):
        raise hcam.HipercamError('Could not find timing file = {} [OK]'.format(tfile))

    # create name of run file and the copy which will only get made
    # later if problems are picked up
    if source == 'hl':
        rfile = run + '.fits'
        cfile = run + '.fits.save'
    else:
        rfile = run + '.dat'
        cfile = run + '.dat.save'

    # first load all old time stamps, and compute MJDs
    if source == 'hl':
        raise NotImplementedError('HiPERCAM option not yet implemented')

    elif source == 'ul':

        # Load up the time stamps from the timing data file. Need also
        # the header of the original run to know how many bytes to
        # read and how to interpret them.
        rhead = hcam.ucam.Rhead(run)
        with open(tfile,'rb') as fin:
            atbytes, mjds, tflags = [], [], []
            nframe = 0
            while 1:
                tbytes = fin.read(rhead.ntbytes)
                if len(tbytes) != rhead.ntbytes:
                    break
                nframe += 1
                atbytes.append(tbytes)

                # interpret times
                mjd, tflag = u_tbytes_to_mjd(tbytes, rhead, nframe)
                mjds.append(mjd)
                tflags.append(tflag)

    # Independent of source, at this stage 'atbytes' is a list of all
    # timestamp bytes, while 'mjds' is a list of all equivalent MJDs,
    # and 'tflags' are bools which are True for OK times, False for
    # null timestamps.
    mjds = np.array(mjds)
    tflags = np.array(tflags)
    inds = np.arange(len(mjds))
    nulls = ~tflags
    nulls_present = nulls.any()

    # Remove null timestamps
    mjds_ok = mjds[tflags]
    inds_ok = inds[tflags]

    # Must have specified minimum of times to work with.
    if len(mjds_ok) < mintim:
        print(
            run,'has too few non-null times to work with ({} vs minimum = {}) [OK]'.format(
                len(mjds_ok), mintim)
        )
        return

    # Median time difference of GPS timestamps
    mdiff = np.median(mjds_ok[1:]-mjds_ok[:-1])

    # Maximum deviation from median separation to allow (days)
    MDIFF = 1.e-10

    # Identify OK timestamps ... (could miss some at this point, but
    # that's OK) kick off by marking all timestamps as False
    ok = mjds_ok == mjds_ok + 1
    for n in range(len(mjds_ok)-1):
        if abs(mjds_ok[n+1]-mjds_ok[n]-mdiff) < MDIFF:
            # if two timestamps differ by the right amount, mark both
            # as ok
            ok[n] = ok[n+1] = True

    gmjds = mjds_ok[ok]
    ginds = inds_ok[ok]

    # Work out integer cycle numbers. First OK timestamp is given
    # cycle number = 0 automatically by the method used.  The cycle
    # numbers can go wrong if the median is not precise enough leading
    # to jumps in the cycle number on runs with large numbers of short
    # exposures, so we build up to the full fit in stages, doubling the
    # number fitted each time

    NMAX = 100
    cycles = (gmjds[:NMAX]-gmjds[0])/mdiff
    moffset = (cycles - np.round(cycles)).mean()
    icycles = np.round(cycles-moffset).astype(int)

    while NMAX < len(gmjds):
        # fit linear trend of first NMAX times where hopefully NMAX is small
        # enough for there not to be a big error but large enough to allow extrapolation.
        print(NMAX)
        slope, intercept, r, p, err = stats.linregress(icycles[:NMAX],gmjds[:NMAX])
        NMAX *= 2
        icycles = np.round((gmjds[:NMAX]-intercept)/slope).astype(int)


    # hope the cycle numbers should be good across the board now final
    # fit
    slope, intercept, r, p, err = stats.linregress(icycles,gmjds)

    # now compute cycle numbers for *all* non-null timestamps
    cycles = (mjds_ok-intercept)/slope
    icycles = np.round(cycles).astype(int)
    cdiffs = cycles-icycles
    monotonic = (icycles[1:]-icycles[:-1] > 0).all()

    # Maximum deviation to allow [cycles]
    CDIFF = 2.e-3

    # if no large differences and cycle number monotically increase,
    # then we are good.
    if not nulls_present and monotonic and (np.abs(cdiffs) < CDIFF).all():
        print(
            'Cycle differences of the {} frames span {:.2e} to {:.2e} [OK]'.format(
                len(icycles),cdiffs.min(),cdiffs.max())
        )
        print(run,'times are OK')
        return

    # next, a common (but perfectly normal) problem is for the final
    # frame to be off in time to be curtailed early giving a final
    # false timestamp which comes early. If the rest are OK, then the
    # run is basically OK.
    bad = np.abs(cdiffs) > CDIFF
    terminated_early = bad[-1]

    if terminated_early and icycles[-1] == icycles[-2]:
        # final frame can end up with same cycle number as one before it
        # apply correction to make later steps easier
        icycles[-1] += 1
        cdiffs[-1] -= 1
        monotonic = (icycles[1:]-icycles[:-1] > 0).all()

    if not nulls_present and monotonic and terminated_early and len(cdiffs[bad]) == 1:
        print('Final frame has cycle difference = {:.2e} [OK]'.format(cdiffs[-1]))
        print(
            'Cycle differences of the {} other frames span {:.2e} to {:.2e} [OK]'.format(
                len(icycles)-1,cdiffs[:-1].min(),cdiffs[:-1].max())
        )
        print(run,'times are OK')
        return

    if plot:
        # diagnostic plot: cycle differences vs cycle numbers
        plt.plot(icycles,cdiffs,'.')
        plt.xlabel('Cycle number')
        plt.ylabel('Cycle difference')
        plt.show()

    # search for duplicate cycle numbers
    u, c = np.unique(icycles, return_counts=True)
    dupes = u[c > 1]
    dupes_present = len(dupes) > 0

    ntot = len(mjds)
    ndupe = len(dupes)
    nnull = len(mjds[nulls])
    bad = np.abs(cdiffs) > CDIFF
    if terminated_early:
        nbad = len(cdiffs[bad])-1
        mdev = np.abs(cdiffs[:-1]).max()
    else:
        nbad = len(cdiffs[bad])
        mdev = np.abs(cdiffs).max()

    # report problems
    print(
        '  {} timestamps are corrupt. TOT,DUP,NULL,BAD = {},{},{},{}; max dev = {:.4f} cycles'.format(
            run, ntot, ndupe, nnull, nbad, mdev)
    )

    if check:
        # go no further, pure check mode
        return

#    if ginds[0] != 0:
#        # this case needs checking
#        raise hcam.HipercamError(
#            'Cannot handle case where first timestamp is no good'
#        )

#    # make copy of data, if not already present
#    if not os.path.exists(cfile):
#        print('Copying',rfile,'to',cfile)
#        shutil.copyfile(rfile, cfile)
#    else:
#        print('Copy of',rfile,'called',cfile,'already exists')


def u_tbytes_to_mjd(tbytes, rtbytes, nframe):
    """Translates set of ULTRACAM timing bytes into an MJD.

    Marks ULTRASPEC null stamps as bad (32 bytes in length,
    the last 20 of which are 0)

    """
    ret = hcam.ucam.utimer(tbytes,rtbytes,nframe)
    if len(tbytes) == 32 and tbytes[12:] == 20*b'\x00':
        return (ret[1]['gps'], False)
    else:
        return (ret[1]['gps'], True)

def h_tbytes_to_mjd(tbytes, nframe):
    """Translates set of HiPERCAM timing bytes into an MJD"""

    # number of seconds in a day
    DAYSEC = 86400.

    frameCount, timeStampCount, years, day_of_year, hours, mins, \
        seconds, nanoseconds, nsats, synced = htimer(tbytes)
    frameCount += 1

    if frameCount != nframe:
        if frameCount == nframe + 1:
            warnings.warn(
                'frame count mis-match; a frame seems to have been dropped'
            )
        else:
            warnings.warn(
                'frame count mis-match; {:d} frames seems to have been dropped'.format(
                    frameCount-self.nframe)
            )

    try:
        imjd = gregorian_to_mjd(years, 1, 1) + day_of_year - 1
        fday = (hours+mins/60+(seconds+nanoseconds/1e9)/3600)/24
    except ValueError:
        imjd = 51544
        fday = nframe/DAYSEC

    return imjd+fday



