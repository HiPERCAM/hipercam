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

        bias    : (string)
           Name of bias frame to subtract, 'none' to ignore.

        adjust  : (string)
           adjustments to make: 'i' = ignore; 'n' = normalise the mean
           of all frames to match the first; 'b' = add offsets so that
           the mean of all frames is the same as the first.
           Option 'n' is useful for twilight flats; 'b' for combining
           biases.

        output : (string)
           output file

    """

    if args is None:
        args = sys.argv[1:]

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'combine', args) as cl:

        # register parameters
        cl.register('list', Cline.GLOBAL, Cline.PROMPT)
        cl.register('bias', Cline.LOCAL, Cline.PROMPT)
        cl.register('adjust', Cline.LOCAL, Cline.PROMPT)
        cl.register('output', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        flist = cl.get_value(
            'list', 'list of files to combine',
            cline.Fname('files', hcam.LIST)
        )

        # bias frame (if any)
        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        if bias is not None:
            # read the bias frame
            bframe = hcam.MCCD.rfits(bias)

        adjust = cl.get_value(
            'adjust', 'i(gnore), n(ormalise) b(ias offsets)',
            'i', lvals=('i','n','b')
        )

        outfile = cl.get_value(
            'output', 'output file',
            cline.Fname('hcam', hcam.HCAM, cline.Fname.NOCLOBBER)
        )

    # Read files into memory, insisting that they
    # all have the same set of CCDs
    print('Loading files from {:s}'.format(flist))
    mccds = []

    with hcam.HcamListSpool(flist) as spool:
        for n, mccd in enumerate(spool):
            if n == 0:
                cnams = set(mccd.keys())
                if bias is not None:
                    # crop the bias on the first frame only
                    bframe.crop(mccd)

            elif set(mccd.keys()) != cnams:
                raise hcam.HipercamError(
                    'Conflicting sets of CCDs'
                    ' found: {!s} vs {!s}'.format(
                        set(mccd.keys()), cnams)
                )

            if bias is not None:
                # subtract bias
                mccd -= bframe

            # keep the result
            mccds.append(mccd)

    print('Loaded {:d} files'.format(len(mccds)))

    # check that all CCDs have at least one valid frame
    for mccd in mccds:
        for cnam, ccd in mccd.items():
            if ccd.is_data() and cnam in cnams:
                cnams.remove(cnam)
        if len(cnams) == 0:
            break
    else:
        raise hcam.HipercamError(
            'Could not find a single valid frame for CCD(s) {:!s}'.format(
                cnams)
            )

    if adjust == 'b' or adjust == 'n':

        # adjust the means
        print('Computing and adjusting the mean levels')

        # first set reference means, a little tricky since we need to be sure
        # that we are dealing with OK data
        means = {}
        cnams = set(mccd.keys())
        for mccd in mccds:
            # add in a mean for each CCD on the first valid frame we
            # encounter
            for cnam, ccd in mccd.items():
                if ccd.is_data() and cnam not in means:
                    means[cnam] = ccd.mean()
                    cnams.remove(cnam)

            if len(cnams) == 0:
                # stop as soon as we can
                break

        # now carry out the adjustments
        for mccd in mccds:
            for cnam, ccd in mccd.items():
                if ccd.is_data():
                    if adjust == 'b':
                        ccd += means[cnam]-ccd.mean()
                    elif adjust == 'n':
                        ccd *= means[cnam]/ccd.mean()

    # Finally, combine
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

