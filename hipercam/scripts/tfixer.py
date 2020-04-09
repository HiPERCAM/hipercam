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
    not seem right, it will do nothing.

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

        check : bool
           True simply to check a run for timing problems, not try to
           fix them.
    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('source', Cline.LOCAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('check', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value(
            'source', 'data source [hl, ul]',
            'hl', lvals=('hl','ul')
        )

        run = cl.get_value('run', 'run name', 'run005')
        if run.endswith('.fits'):
            run = run[:-5]

        check = cl.get_value('check', 'only check the run for problems', True)

    # create name of timing file
    tfile = os.path.join('tbytes', os.path.basename(run) + hcam.TBTS)
    if not os.path.isfile(tfile):
        raise hcam.HipercamError('Could not find timing file =',tfile)

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

    # Remove null timestamps
    mjds_ok = mjds[tflags]
    inds_ok = inds[tflags]

    # Median time difference of GPS timestamps
    mdiff = np.median(mjds_ok[1:]-mjds_ok[:-1])

    # Maximum deviation from median separation to allow (days)
    MDIFF = 1.e-9

    # Assuming mdiff is right, try to identify definitely OK
    # timestamps ... (could miss some at this point, but that's OK)
    # kick off by marking all timestamp as False
    ok = mjds_ok == mjds_ok + 1
    for n in range(len(mjds_ok)-1):
        if abs(mjds_ok[n+1]-mjds_ok[n]-mdiff) < MDIFF:
            # if two timestamps differ by the right amount, mark both
            # as ok
            ok[n] = ok[n+1] = True

    gmjds = mjds_ok[ok]
    ginds = inds_ok[ok]

    if len(mjds) == 1 or len(mjds) == len(gmjds):
        print(
            run,'has',len(mjds),'time stamps, no null or bad: >> RUN OK <<',
            file=sys.stderr
        )
    else:
        print(
            run,'has',len(mjds),'time stamps.',
            len(mjds[~tflags]),'null,',
            len(mjds_ok[~ok]),'bad. Null frames =',inds[~tflags]+1,
            ', bad frames & times =',inds_ok[~ok],'&',mjds_ok[~ok],
            '>> RUN NEEDS FIXING <<'
        )

    if check:
        return

    if ginds[0] != 0:
        # this case needs checking
        raise hcam.HipercamError('Cannot handle case where first timestamp is no good')

    # Now work out integer cycle numbers for all OK timestamps.
    # First OK timestamp is given cycle number = 0 automatically
    # by the method used.
    cycles = (gmjds-gmjds[0])/mdiff
    icycles = np.round(cycles).astype(int)

    # fit and remove linear trend
    slope, intercept, r, p, err = stats.linregress(icycles,gmjds)
    pmjds = slope*icycles + intercept
    diffs = gmjds-pmjds

    plt.hist(86400*diffs,51)
    plt.xlabel('Difference [secs]')
    plt.ylabel('Number')
    plt.show()

    plt.plot(icycles,86400*diffs,'.')
    plt.xlabel('Cycle number')
    plt.ylabel('Difference [secs]')
    plt.show()

    FRAC = 1.e-3
    if (abs(diffs) > FRAC*slope).any():
        iworst = abs(diffs).argmax()
        raise hcam.HipercamError(
            'Some timestamps deviate more than expected. Worst case = ',
            diffs[iworst],'days cf limit =',FRAC*slope,'which is',FRAC,
            'of mean cycle time'
        )

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
            warnings.warn('frame count mis-match; a frame seems to have been dropped')
        else:
            warnings.warn(
                'frame count mis-match; {:d} frames seems to have been dropped'.format(frameCount-self.nframe)
            )

    try:
        imjd = gregorian_to_mjd(years, 1, 1) + day_of_year - 1
        fday = (hours+mins/60+(seconds+nanoseconds/1e9)/3600)/24
    except ValueError:
        imjd = 51544
        fday = nframe/DAYSEC

    return imjd+fday



