import sys
import os

import numpy as np

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ = ['combine',]

############################
#
# combine -- combines images
#
############################

def combine(args=None):
    """Combines a series of images defined by a list using median
    combination. Only combines those CCDs for which is_data() is true (i.e. it
    skips blank frames caused by NSKIP / NBLUE options)

    Arguments::

        list   : (string)
           list of hcm files with images to combine. The formats of the
           images should all match

        output : (string)
           output file

    """

    if args is None:
        args = sys.argv[1:]

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'combine', args) as cl:

        # register parameters
        cl.register('list', Cline.GLOBAL, Cline.PROMPT)
        cl.register('output', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        flist = cl.get_value(
            'list', 'list of files to combine',
            cline.Fname('files', hcam.LIST)
        )

        outfile = cl.get_value(
            'output', 'output file',
            cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW)
        )

    # Read files into memory
    print('Loading files from {:s}'.format(flist))
    mccds = []
    with hcam.HcamListSpool(flist) as spool:
        for mccd in spool:
            mccds.append(mccd)
    print('Loaded {:d} files'.format(len(mccds)))

    # Combine
    ofile = mccds[0].copy()
    for cnam, occd in ofile.items():
        print('Combining CCD {:s} frames'.format(cnam))
        for wnam, owin in occd.items():
            # build lists of all windows, but only for those
            # CCD frames that return is_data() as True
            arrs = []
            for mccd in mccds:
                ccd = mccd[cnam]
                if ccd.is_data():
                    arrs.append(ccd[wnam].data)
            arr3d = np.stack(arrs)
            owin.data = np.median(arr3d,axis=0)

    # write out
    ofile.wfits(outfile)
    print('Written result to {:s}'.format(outfile))

