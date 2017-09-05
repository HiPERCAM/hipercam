"""Command line script to grab images"""

import sys
import os

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

      source : (string)
         's' = server, 'l' = local

      inst   : (string) [hidden]
         Instrument involved. Two choices, 'u' for ULTRCAM or ULTRASPEC, 'h'
         for HiPERCAM. This is needed because of the different formats.

      run    : (string)
         run name to access

      ndigit : (int)
         Files created will be written as 'run005_0013.fits' etc. `ndigit` is
         the number of digits used for the frame number (4 in this case)

      first  : (int)
         First frame to access

      last   : (int)
         Last frame to access, 0 for the lot

      bias   : (string)
         Name of bias frame to subtract, 'none' to ignore.

      dtype  : (string)
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
        cl.register('source', Cline.LOCAL, Cline.HIDE)
        cl.register('inst', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ndigit', Cline.LOCAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('last', Cline.LOCAL, Cline.PROMPT)
        cl.register('bias', Cline.GLOBAL, Cline.PROMPT)
        cl.register('dtype', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value(
            'source', 'data source [s(erver), l(ocal)]',
            'l', lvals=('s','l')
        )

        inst = cl.get_value(
            'inst', 'instrument used [h(ipercam), u(ltracam/spec)]',
            'h', lvals=('h','u')
        )

        run = cl.get_value('run', 'run name', 'run005')

        ndigit = cl.get_value(
            'ndigit', 'number of digits in frame identifier', 3, 0)

        first = cl.get_value('first', 'first frame to grab', 1, 1)
        last = cl.get_value('last', 'last frame to grab', 0)
        if last < first and last != 0:
            sys.stderr.write('last must be >= first or 0')
            sys.exit(1)

        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        if bias is not None:
            # read the bias frame
            bframe = hcam.MCCD.rfits(bias)

        dtype = cl.get_value(
            'dtype', 'data type [r(aw), f(32-bit float), i(16-bit uint)]',
            'f', lvals=('r','f','i')
        )

    # Now the actual work. First set up the arguments for data_source
    server = source == 's'
    flist = None
    if inst == 'u':
        instrument = 'ULTRA'
    elif inst == 'h':
        instrument = 'HIPER'

    # strip off extensions
    if run.endswith(hcam.HRAW):
        run = run[:run.find(hcam.HRAW)]

    # Finally, we can go
    with hcam.data_source(instrument, run, flist, server, first) as spool:
        nframe = first
        root = os.path.basename(run)
        for frame in spool:

            # subtract bias
            if bias is not None:
                frame -= bframe

            if dtype == 'i':
                frame.uint16()
            elif dtype == 'f':
                frame.float32()

            # write to disk
            fname = '{:s}_{:0{:d}}{:s}'.format(run,nframe,ndigit,hcam.HCAM)
            frame.wfits(fname,True)

            print('Written frame {:d} to {:s}'.format(nframe,fname))
            nframe += 1
            if last and nframe > last: break

