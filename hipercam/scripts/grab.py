"""Command line script to grab images"""

import sys
import os
import time

import numpy as np

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ =  ['grab',]

###########################################################
#
# grab -- downloads a series of images from a raw data file
#
###########################################################

def grab(args=None):
    """This downloads a sequence of images from a raw data file and writes them
    out to a series CCD / MCCD files.

    Arguments::

        source  : (string) [hidden]
           Data source, four options::

               'hs' : HiPERCAM server
               'hl' : local HiPERCAM FITS file
               'us' : ULTRACAM server
               'ul' : local ULTRACAM .xml/.dat files

           'hf' useful when looking at a set of frames generated
           by 'grab' or converted from some foreign data format.

      run    : (string)
         run name to access

      ndigit : (int)
         Files created will be written as 'run005_0013.fits' etc. `ndigit` is
         the number of digits used for the frame number (4 in this case)

      first  : (int)
         First frame to access

      last   : (int)
         Last frame to access, 0 for the lot

      twait  : (float) [hidden]
         time to wait between attempts to find a new exposure, seconds.

      tmax  : (float) [hidden]
         maximum time to wait between attempts to find a new exposure, seconds.

      bias   : (string)
         Name of bias frame to subtract, 'none' to ignore.

      dtype  : (string) [hidden, defaults to 'f']
         Data type on output. Options::

            'r' : "raw", do nothing; this is safe but could
                  be profligate as you could end up storing 64-bit
                  double precision values if you subtract a bias.

            'f' : Convert to 32-bit floats if already in 64-bit form. Leave
                  integers untouched. This implies some loss of precision,
                  but requires half the space.

            'i' : Convert to 16-bit unsigned integers. A warning will be
                  issued if loss of precision occurs; an exception will
                  be raised if the data are outside the range 0 to 65535.
    """

    if args is None:
        args = sys.argv[1:]

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'grab', args) as cl:

        # register parameters
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ndigit', Cline.LOCAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('last', Cline.LOCAL, Cline.PROMPT)
        cl.register('twait', Cline.LOCAL, Cline.HIDE)
        cl.register('tmax', Cline.LOCAL, Cline.HIDE)
        cl.register('bias', Cline.GLOBAL, Cline.PROMPT)
        cl.register('dtype', Cline.LOCAL, Cline.HIDE)

        # get inputs
        source = cl.get_value('source', 'data source [hs, hl, us, ul]',
                              'hl', lvals=('hs','hl','us','ul'))

        # set some flags
        server_or_local = source.endswith('s') or source.endswith('l')
        is_file_list = False
        server_on = source.endswith('s')

        # distinguish between the instruments HiPERCAM or
        # ULTRA(CAM/SPEC)
        instrument = 'HIPER' if source.startswith('h') else 'ULTRA'

        # OK, more inputs
        run = cl.get_value('run', 'run name', 'run005')

        ndigit = cl.get_value(
            'ndigit', 'number of digits in frame identifier', 3, 0)

        first = cl.get_value('first', 'first frame to grab', 1, 1)
        last = cl.get_value('last', 'last frame to grab', 0)
        if last < first and last != 0:
            sys.stderr.write('last must be >= first or 0')
            sys.exit(1)

        twait = cl.get_value(
            'twait', 'time to wait for a new frame [secs]', 1., 0.)
        tmax = cl.get_value(
            'tmax', 'maximum time to wait for a new frame [secs]', 10., 0.)

        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        if bias is not None:
            # read the bias frame
            bframe = hcam.MCCD.rfits(bias)

        cl.set_default('dtype', 'f')
        dtype = cl.get_value(
            'dtype', 'data type [r(aw), f(32-bit float), i(16-bit uint)]',
            'f', lvals=('r','f','i')
        )

    # Now the actual work.

    # strip off extensions
    if run.endswith(hcam.HRAW):
        run = run[:run.find(hcam.HRAW)]

    # initialisations
    total_time = 0 # time waiting for new frame
    nframe = first
    root = os.path.basename(run)

    # Finally, we can go
    with hcam.data_source(
            instrument, run, is_file_list, server_on, first) as spool:

        for mccd in spool:

            # Handle the waiting game ...
            give_up, try_again, total_time = hcam.hang_about(
                mccd, twait, tmax, total_time
            )

            if give_up:
                print('grab stopped')
                break
            elif try_again:
                continue

            # subtract bias
            if bias is not None:
                mccd -= bframe

            if dtype == 'i':
                mccd.uint16()
            elif dtype == 'f':
                mccd.float32()

            # write to disk
            fname = '{:s}_{:0{:d}}{:s}'.format(run,nframe,ndigit,hcam.HCAM)
            mccd.wfits(fname,True)

            print('Written frame {:d} to {:s}'.format(nframe,fname))
            nframe += 1
            if last and nframe > last:
                print('grab stopped')
                break

