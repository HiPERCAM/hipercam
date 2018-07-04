"""Command line script to grab images"""

import sys
import os
import time
import tempfile
import getpass

import numpy as np

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ =  ['grab',]

###########################################################
#
# grab -- downloads a series of images from a raw data file
#
###########################################################

def grab(args=None):
    """``grab [source] run [temp] (ndigit) first last [twait tmax] bias
    [dtype]``

    This downloads a sequence of images from a raw data file and writes them
    out to a series CCD / MCCD files.

    Parameters:

       source  : string [hidden]
           Data source, four options:

              | 'hs' : HiPERCAM server
              | 'hl' : local HiPERCAM FITS file
              | 'us' : ULTRACAM server
              | 'ul' : local ULTRACAM .xml/.dat files

       run     : string
           run name to access

       temp    : bool [hidden, defaults to False]
           True to indicate that the frames should be written to temporary
           files with automatically-generated names in the expectation of
           deleting them later. This also writes out a file listing the names.
           The aim of this is to allow grab to be used as a feeder for other
           routines in larger scripts.  If temp == True, grab will return with
           the name of the list of hcm files. Look at 'makebias' for an
           example that uses this feature.

       ndigit  : int [if not temp]
           Files created will be written as 'run005_0013.fits' etc. `ndigit`
           is the number of digits used for the frame number (4 in this
           case). Any pre-existing files of the same name will be
           over-written.

       first   : int
           First frame to access

       last    : int
           Last frame to access, 0 for the lot

       twait   : float [hidden]
           time to wait between attempts to find a new exposure, seconds.

       tmax    : float [hidden]
           maximum time to wait between attempts to find a new exposure,
           seconds.

       bias    : string
           Name of bias frame to subtract, 'none' to ignore.

       dtype   : string [hidden, defaults to 'f32']
           Data type on output. Options:

            | 'f32' : output as 32-bit floats (default)

            | 'f64' : output as 64-bit floats.

            | 'u16' : output as 16-bit unsigned integers. A warning will be
                      issued if loss of precision occurs; an exception will
                      be raised if the data are outside the range 0 to 65535.

    """

    command, args = utils.script_args(args)

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('temp', Cline.GLOBAL, Cline.HIDE)
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

        # OK, more inputs
        resource = cl.get_value('run', 'run name', 'run005')

        cl.set_default('temp', False)
        temp = cl.get_value('temp', 'save to temporary automatically-generated file names?', True)

        if not temp:
            ndigit = cl.get_value(
                'ndigit', 'number of digits in frame identifier', 3, 0)

        first = cl.get_value('first', 'first frame to grab', 1, 0)
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

        cl.set_default('dtype', 'f32')
        dtype = cl.get_value(
            'dtype', 'data type [f32, f64, u16]',
            'f32', lvals=('f32','f64','u16')
        )

    # Now the actual work.

    # strip off extensions
    if resource.endswith(hcam.HRAW):
        resource = resource[:resource.find(hcam.HRAW)]

    # initialisations
    total_time = 0 # time waiting for new frame
    nframe = first
    root = os.path.basename(resource)
    bframe = None

    # Finally, we can go
    if temp:
        # create a directory on temp for the temporary file to avoid polluting
        # it too much
        fnames = []
        tdir = os.path.join(tempfile.gettempdir(), 'hipercam-{:s}'.format(getpass.getuser()))
        os.makedirs(tdir,exist_ok=True)

    with spooler.data_source(source, resource, first) as spool:

        try:

            for mccd in spool:

                # Handle the waiting game ...
                give_up, try_again, total_time = spooler.hang_about(
                    mccd, twait, tmax, total_time
                    )

                if give_up:
                    print('grab stopped')
                    break
                elif try_again:
                    continue

                if bias is not None:
                    # read bias after first frame so we can
                    # chop the format
                    if bframe is None:

                        # read the bias frame
                        bframe = hcam.MCCD.read(bias)

                        # reformat
                        bframe = bframe.crop(mccd)

                    mccd -= bframe

                if dtype == 'u16':
                    mccd.uint16()
                elif dtype == 'f32':
                    mccd.float32()
                elif dtype == 'f64':
                    mccd.float64()

                # write to disk
                if temp:
                    # generate name automatically
                    fd, fname = tempfile.mkstemp(suffix=hcam.HCAM, dir=tdir)
                    mccd.write(fname,True)
                    os.close(fd)
                    fnames.append(fname)
                else:
                    fname = '{:s}_{:0{:d}}{:s}'.format(root,nframe,ndigit,hcam.HCAM)
                    mccd.write(fname,True)

                print('Written frame {:d} to {:s}'.format(nframe,fname))


                # update the frame number
                nframe += 1
                if last and nframe > last:
                    break

        except KeyboardInterrupt:
            # trap ctrl-C so we can delete temporary files if temp
            if temp:
                for fname in fnames:
                    os.remove(fname)
                print('\ntemporary files deleted')
                print('grab aborted')
            else:
                print('\ngrab aborted')
            sys.exit(1)


    if temp:
        if len(fnames) == 0:
            raise hcam.HipercamError(
                'no files were grabbed; please check input parameters, especially "first"'
                )
        # write the file names to a list
        fd, fname = tempfile.mkstemp(suffix=hcam.LIST,dir=tdir)
        with open(fname,'w') as fout:
            for fnam in fnames:
                fout.write(fnam + '\n')
        os.close(fd)
        print('temporary file names written to {:s}'.format(fname))

        # return the name of the file list
        return fname
