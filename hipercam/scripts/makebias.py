"""Command line script to grab images"""

import sys
import os
import time
import tempfile

import numpy as np

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ =  ['grab',]

#########################################################
#
# makebias -- combines frames of a run into one using
# clipped mean averaging as appropriate for bias frames.
#
#########################################################

def makebias(args=None):
    """Combines the frames of a single run into a single frame using
    clipped-mean averaging appropriate for biases. This uses 'grab' to get the
    frames and 'combine' to combine them. If you already have the frames
    separately, then just use 'combine'.

    Arguments::

        source  : string [hidden]
           Data source, four options::

               'hs' : HiPERCAM server
               'hl' : local HiPERCAM FITS file
               'us' : ULTRACAM server
               'ul' : local ULTRACAM .xml/.dat files

      run    : string
         run name to access

      first  : int
         First frame to access

      last   : int
         Last frame to access, 0 for the lot

      sigma  : float
         The value of 'sigma' to pass to the clipped mean combination in
         'combine'

      plot   : bool
         make a plot of the mean level versus frame number for each CCD. This
         can provide a quick check that the frames are not too different. You will
         need explicitly to close the plot generated at the end of the script

      twait  : float [hidden]
         time to wait between attempts to find a new exposure, seconds.

      tmax   : float [hidden]
         maximum time to wait between attempts to find a new exposure, seconds.

      output : string
         name of final combined file

    .. notes::

       This routine writes the files returned by 'grab' to automatically generated
       files, typically in /tmp. These are removed at the end, but won't be if you
       ctrl-C. You might to check your /tmp disk for redundant files.
    """

    command, args = utils.script_args(args)

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('last', Cline.LOCAL, Cline.PROMPT)
        cl.register('sigma', Cline.LOCAL, Cline.PROMPT)
        cl.register('plot', Cline.LOCAL, Cline.PROMPT)
        cl.register('twait', Cline.LOCAL, Cline.HIDE)
        cl.register('tmax', Cline.LOCAL, Cline.HIDE)
        cl.register('output', Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value('source', 'data source [hs, hl, us, ul]',
                              'hl', lvals=('hs','hl','us','ul'))

        run = cl.get_value('run', 'run name', 'run005')

        first = cl.get_value('first', 'first frame to grab', 1, 0)
        last = cl.get_value('last', 'last frame to grab', 0)
        if last < first and last != 0:
            sys.stderr.write('last must be >= first or 0')
            sys.exit(1)

        sigma = cl.get_value(
            'sigma', 'number of RMS deviations to clip', 3., 1.
            )

        plot = cl.get_value(
            'plot', 'plot mean levels versus frame number?',
            False
            )

        twait = cl.get_value(
            'twait', 'time to wait for a new frame [secs]', 1., 0.)
        tmax = cl.get_value(
            'tmax', 'maximum time to wait for a new frame [secs]', 10., 0.)

        output = cl.get_value('output', 'output name', 'bias')

    # Now the actual work. 

    # We pass full argument lists to grab and combine because with None as the
    # command name, the default file mechanism is by-passed. 'PROMPT' is used
    # to expose the hidden parameters. Every argument needed is passed to
    # avoid any interactive prompting. Note that because of the Cline mechanisms,
    # in particular prompting of one parameter depending on another, it can be tricky
    # to get this right. The keyword 'LIST' can be helpful in case of difficulty. i.e.
    # rather than just 'PROMPT' you would use 'PROMPT','LIST'

    print("\nCalling 'grab' ...")
    args = [None,'PROMPT',source,run,'yes',str(first),str(last),str(twait),str(tmax),'none','f32']
    flist = hcam.scripts.grab(args)

    try:

        print("\nCalling 'combine' ...")
        args = [None,'PROMPT',flist,'none','c',str(sigma),'b','yes' if plot else 'no', 'yes', output]
        hcam.scripts.combine(args)

        # remove temporary files
        with open(flist) as fin:
            for fname in fin:
                fname = fname.strip()
                os.remove(fname)
        os.remove(flist)
        print('temporary files have been deleted')
        print('makebias finished')

    except KeyboardInterrupt:
        # this to ensure we delete the temporary files
        with open(flist) as fin:
            for fname in fin:
                fname = fname.strip()
                os.remove(fname)
        os.remove(flist)
        print('\ntemporary files have been deleted')
        print('makebias aborted')
