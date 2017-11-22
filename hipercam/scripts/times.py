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

######################################
#
# times -- lists the times of multiple images
#
######################################

def times(args=None):
    """Lists timing information from a run. For each frame the associated GPS
    timestamp is listed (UTC, date then HMS), and then for each CCD the CCD
    number [1-5], the mid exposure time (UTC, HMS only), the exposure length
    (seconds) and a status flag are listed in parentheses (). The status flag
    is set by the NSKIP parameters for each CCD. Only the valid frames will have
    status = 'T'.

    Arguments::

        source  : string [hidden]
           Data source, two options::

               'hs' : HiPERCAM server
               'hl' : local HiPERCAM FITS file

        run     : string
           run number to access, e.g. 'run0034'

        first   : int
           exposure number to start from. 1 = first frame; set = 0 to
           always try to get the most recent frame (if it changes)

        last    : int
           last exposure number, 0 for all.

        twait   : float [hidden]
           time to wait between attempts to find a new exposure, seconds.

        tmax    : float [hidden]
           maximum time to wait between attempts to find a new exposure,
           seconds. Set to 0 to avoid waiting and extra output messages.



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

        # get inputs
        source = cl.get_value('source', 'data source [hs, hl]',
                              'hl', lvals=('hs','hl'))

        resource = cl.get_value('run', 'run name', 'run005')
        first = cl.get_value('first', 'first frame to list', 1, 0)
        last = cl.get_value('last', 'last frame to list', 0, 0)
        if last < first and last != 0:
            sys.stderr.write('last must be >= first or 0')
            sys.exit(1)

        twait = cl.get_value(
            'twait', 'time to wait for a new frame [secs]', 1., 0.)
        tmax = cl.get_value(
            'tmax', 'maximum time to wait for a new frame [secs]', 10., 0.)

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff

    # open the run as an Rtime
    rtime = hcam.hcam.Rtime(resource, first, source.endswith('s'))

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
        tstamp, tinfo = tdata
        print(
            'Frame {:d}, GPS = {:s}'.format(
                nframe, tstamp.iso), end=''
            )


        message = ''
        for nccd, (mjd, exptime, flag) in enumerate(tinfo):
            ts = Time(mjd, format='mjd', precision=6)
            message += ', ({:d}, {:s}, {:.4f}, {:s})'.format(nccd+1, ts.hms_custom, exptime, 'T' if flag else 'F')
        print(message)

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
