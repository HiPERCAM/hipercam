"""Command line script to create a dark frame"""

import sys
import os
import time
import tempfile
import warnings
import numpy as np

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ =  ['makedark',]

############################################################################
#
# makedark -- combines frames of a run into one using clipped mean averaging,
# with bias subtraction. The exposure of the dark is corrected by whatever
# exposure applies to the bias.
#
#############################################################################

def makedark(args=None):
    """``makedark [source] run first last bias sigma plot [twait tmax output]``

    Combines the frames of a single run into a single frame using clipped-mean
    averaging. This uses ``grab`` to get the frames and ``combine`` to combine
    them. It subtracts a bias and corrects the exposure time by the exposure
    time of the bias.

    Parameters:

       source  : string [hidden]
           Data source, four options:

             |   'hs' : HiPERCAM server
             |   'hl' : local HiPERCAM FITS file
             |   'us' : ULTRACAM server
             |   'ul' : local ULTRACAM .xml/.dat files

       run : string
           run name to access

       first : int
           First frame to access

       last : int
           Last frame to access, 0 for the lot

       bias : string
           Names of bias frame (made e.g. with |makebias|). This is
           so that the counts left in the dark frame are genuine dark
           counts which can be scaled by the ratio of exposure lengths
           during dark subtraction.

       sigma   : float
           The value of 'sigma' to pass to the clipped mean combination in
           'combine'

       plot    : bool
           make a plot of the mean level versus frame number for each
           CCD. This can provide a quick check that the frames are not too
           different. You will need explicitly to close the plot generated at
           the end of the script

       twait   : float [hidden]
           time to wait between attempts to find a new exposure, seconds.

       tmax    : float [hidden]
           maximum time to wait between attempts to find a new exposure,
           seconds.

       output  : string
           name of final combined file

      .. Note:: 

         This routine writes the files returned by 'grab' to automatically
         generated files, typically in /tmp, to avoid polluting the working
         directory. These are removed at the end, but won't be if you
         ctrl-C. You should check your /tmp disk for redundant files in this
         case.

    """

    command, args = utils.script_args(args)

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('last', Cline.LOCAL, Cline.PROMPT)
        cl.register('bias', Cline.LOCAL, Cline.PROMPT)
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

        bias = cl.get_value('bias', 'bias to subtract', 'bias')
        
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
    # command name, the default file mechanism is by-passed. 'prompt' is used
    # to expose the hidden parameters. Every argument needed is passed to
    # avoid any interactive prompting. Note that because of the Cline
    # mechanisms, in particular prompting of one parameter depending on
    # another, it can be tricky to get this right. The keyword 'list' can be
    # helpful in case of difficulty. i.e.  rather than just 'prompt' you would
    # use 'prompt', 'list'

    print("\nCalling 'grab' ...")
    args = [
        None, 'prompt', source, run, 'yes',
        str(first), str(last), str(twait), str(tmax),
        bias, 'f32'
    ]
    flist = hcam.scripts.grab(args)

    if first == 1:
        # test readout mode if the first == 1 as, with non clear modes, the
        # first file is different from all others. A warning is issued.
        with open(flist) as f:
            first_frame = f.readline().strip()

        mccd = hcam.MCCD.read(first_frame)
        instrument = mccd.head.get('INSTRUME','UNKNOWN')
        if instrument == 'ULTRACAM' or instrument == 'HIPERCAM' or \
           instrument == 'ULTRASPEC':
            if 'CLEAR' in mccd.head:
                if not mccd.head['CLEAR']:
                    warnings.warn(
                        """You should not include the first frame of a run when making a dark from
readout modes which do not have clear enabled since the first frame is
different from all others."""
                    )
            else:
                warnings.warn(
                        """Instrument = {:s} has readout modes with both clears enabled or not
between exposures. When no clear is enabled, the first frame is different
from all others and should normally not be included when making a dark.
This message is a temporary stop gap until the nature of the readout mode
has been determined with respect to clears."""
                    )
                
    try:

        print("\nCalling 'combine' ...")
        args = [
            None, 'prompt', flist, 'none', 'c', str(sigma),
            'i', 'yes', 'yes' if plot else 'no', 'yes', output
        ]
        hcam.scripts.combine(args)

        # remove temporary files
        with open(flist) as fin:
            for fname in fin:
                fname = fname.strip()
                os.remove(fname)
        os.remove(flist)
        print('temporary files have been deleted')

        # correct exposure time of dark frame by the exposure time of
        # the bias frame used
        dark = hcam.MCCD.read(add_extension(output,hcam.HCAM))
        bias = hcam.MCCD.read(add_extension(bias,hcam.HCAM))
        if 'EXPTIME' in dark.head and 'EXPTIME' in bias.head:
            dexpose = dark.head['EXPTIME']
            bexpose = bias.head['EXPTIME']
            dark.head['EXPTIME'] = dexpose-bexpose
            print('Corrected dark exposure time from {:.2f} to {:.2f}'.format(
                dexpose,bexpose)
            )
            dark.write(add_extension(output,hcam.HCAM), True)
        else:
            warnings.warn(
                'Could not find exposure time (EXPTIME) in the dark and/or'
                ' the bias hence could not correct it in the dark'
            )
        
        print('makedark finished')

    except KeyboardInterrupt:
        # this to ensure we delete the temporary files
        with open(flist) as fin:
            for fname in fin:
                fname = fname.strip()
                os.remove(fname)
        os.remove(flist)
        print('\ntemporary files have been deleted')
        print('makedark aborted')
