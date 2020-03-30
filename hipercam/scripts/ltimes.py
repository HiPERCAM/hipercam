import sys
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time, TimeISO

from trm.pgplot import *
import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ = ['times',]

##############################################
#
# ltimes -- lists the times of multiple images
#
##############################################

def ltimes(args=None):
    """``ltimes [source] run first last [twait tmax tdigit edigit]``

    Lists timing information from a run. For each frame the associated GPS
    timestamp is listed (UTC, date then HMS), and then for each CCD the CCD
    number [1-5], the mid exposure time (UTC, HMS only), the exposure length
    (seconds) and a status flag are listed in parentheses (). The status flag
    is set by the NSKIP parameters for each CCD. Only the valid frames will
    have status = 'T'.

    Parameters:

        source : string [hidden]
           Data source, two options:

              | 'hs' : HiPERCAM server
              | 'hl' : local HiPERCAM FITS file
              | 'us' : local ULTRA(CAM|SPEC) file
              | 'ul' : ULTRA(CAM|SPEC) server

        run : string
           run number to access, e.g. 'run0034'

        first : int
           exposure number to start from. 1 = first frame; set = 0 to
           always try to get the most recent frame (if it changes)

        last : int
           last exposure number, 0 for all.

        twait : float [hidden]
           time to wait between attempts to find a new exposure, seconds.

        tmax : float [hidden]
           maximum time to wait between attempts to find a new exposure,
           seconds. Set to 0 to avoid waiting and extra output messages.

        tdigit : int [hidden]
           number of digits of precision for the seconds after the decimal
           point when reporting the times.

        edigit : int [hidden]
           number of digits of precision after the decimal point when
           reporting the exposure times

    .. Warning::

       This routine does not yet work with ULTRACAM data.

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('last', Cline.LOCAL, Cline.PROMPT)
        cl.register('twait', Cline.LOCAL, Cline.HIDE)
        cl.register('tmax', Cline.LOCAL, Cline.HIDE)
        cl.register('tdigit', Cline.LOCAL, Cline.HIDE)
        cl.register('edigit', Cline.LOCAL, Cline.HIDE)

        # get inputs
        source = cl.get_value('source', 'data source [hs, hl, us, ul]',
                              'hl', lvals=('hs','hl','us','ul'))

        resource = cl.get_value('run', 'run name', 'run005')
        if source == 'hs':
            first = cl.get_value('first', 'first frame to plot', 1)
        else:
            first = cl.get_value('first', 'first frame to plot', 1, 0)
        last = cl.get_value('last', 'last frame to list', 0, 0)
        if last < first and last != 0:
            sys.stderr.write('last must be >= first or 0')
            sys.exit(1)

        twait = cl.get_value(
            'twait', 'time to wait for a new frame [secs]', 1., 0.)
        tmax = cl.get_value(
            'tmax', 'maximum time to wait for a new frame [secs]', 10., 0.)

        tdigit = cl.get_value(
            'tdigit', 'digits after decimal point for times', 6, 1, 9
        )
        edigit = cl.get_value(
            'edigit', 'digits after decimal point for exposure times', 3, 1, 9
        )

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff

    # open the run as an Rtime with exact type depending on the source
    if source.startswith('h'):
        rtime = hcam.hcam.Rtime(resource, first, source.endswith('s'))
    else:
        rtime = hcam.ucam.Rtime(resource, first, source.endswith('s'))

    total_time = 0
    nframe = first
    for tdata in rtime:

        # Handle the waiting game ...
        give_up, try_again, total_time = spooler.hang_about(
            tdata, twait, tmax, total_time
            )

        if give_up:
            if tmax > 0:
                print('times stopped')
            break
        elif try_again:
            continue

        # indicate progress
        if len(tdata) == 3:
            # HiPERCAM
            tstamp, tinfo, tflag = tdata
            tstamp.precision = tdigit
            print(
                'Frame {:d}, GPS = {:s} [{:s}]'.format(
                    nframe, tstamp.iso, 'OK' if tflag else 'NOK'), end=''
            )

            message = ''
            for nccd, (mjd, exptime, flag) in enumerate(tinfo):
                ts = Time(mjd, format='mjd', precision=tdigit)
                message += ', ({:d}, {:s}, {:.{:d}f}, {:s})'.format(
                    nccd+1, ts.hms_custom, exptime, edigit, 'T' if flag else 'F'
                )
            print(message)

        elif len(tdata) == 2:
            # ULTRASPEC
            time, tinfo = tdata
            ts = Time(time.mjd, format='mjd', precision=tdigit)
            print(
                'Frame {:d}, HMS = {:s} [{:s}], exp = {:.5f}, MJD = {:.10f}'.format(
                    nframe, ts.hms_custom, 'OK' if time.good else 'NOK',
                    time.expose, time.mjd)
            )

        elif len(tdata) == 4:
            # ULTRACAM
            time, tinfo,btime,bgood = tdata
            ts = Time(time.mjd, format='mjd', precision=tdigit)
            bts = Time(time.mjd, format='mjd', precision=tdigit)
            print(
                'Frame {:d}, HMS = {:s} [{:s}], exp = {:.5f}, MJD = {:.9f}'.format(
                    nframe, ts.hms_custom, 'OK' if time.good else 'NOK',
                    time.expose, time.mjd)
            )

        # increment the frame
        nframe += 1

class TimeHMSCustom(TimeISO):
    """
    Just the HMS part of a Time as "<HH>:<MM>:<SS.sss...>".
    """
    name = 'hms_custom'  # Unique format name
    subfmts = (
        ('date_hms',
         '%H%M%S',
         '{hour:02d}:{min:02d}:{sec:02d}'),
        )
