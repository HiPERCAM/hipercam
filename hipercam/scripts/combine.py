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

        clobber : (bool) [hidden]
           clobber any pre-existing output files

        output  : (string)
           output file

    Notes: this routine reads all inputs into memory, so can be a bit of
    a hog. However, it does so one CCD at a time to alleviate this. It will
    fail if it cannot find a valid frame for any CCD
    """

    if args is None:
        args = sys.argv[1:]

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'combine', args) as cl:

        # register parameters
        cl.register('list', Cline.GLOBAL, Cline.PROMPT)
        cl.register('bias', Cline.LOCAL, Cline.PROMPT)
        cl.register('adjust', Cline.LOCAL, Cline.PROMPT)
        cl.register('clobber', Cline.LOCAL, Cline.HIDE)
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

        clobber = cl.get_value(
            'clobber', 'clobber any pre-existing files on output',
            False
        )

        outfile = cl.get_value(
            'output', 'output file',
            cline.Fname(
                'hcam', hcam.HCAM,
                cline.Fname.NEW if clobber else cline.Fname.NOCLOBBER
            )
        )

    # inputs done with

    # Read the first file of the list to act as a template
    # for the CCD names etc.
    with open(flist) as fin:
        for line in fin:
            if not line.startswith('#') and not line.isspace():
                template_name = line.strip()
                break
        else:
            raise hcam.HipercamError(
                'List = {:s} is empty'.format(flist)
            )

    template = hcam.MCCD.rfits(
        hcam.add_extension(template_name,hcam.HCAM)
    )

    if bias is not None:
        # crop the bias
        bframe.crop(template)

    # Now process each file CCD by CCD to reduce the memory
    # footprint
    for cnam in template:

        # Read files into memory, insisting that they
        # all have the same set of CCDs
        print(
            "\nLoading all CCDs labelled '{:s}' from {:s}".format(cnam, flist)
        )

        ccds = []
        with hcam.HcamListSpool(flist, cnam) as spool:

            if bias is not None:
                # extract relevant CCD from the bias
                bccd = bframe[cnam]

            mean = None
            for ccd in spool:

                if ccd.is_data():
                    if bias is not None:
                        # subtract bias
                        ccd -= bccd

                    # keep the result
                    ccds.append(ccd)

                    # store the first mean
                    if mean is None:
                        mean = ccd.mean()

            if len(ccds) == 0:
                raise hcam.HipercamError(
                    'Found no valid examples of CCD {:s}'
                    ' in list = {:s}'.format(cnam,flist)
                )

            else:
                print('Loaded {:d} CCDs'.format(len(ccds)))

        if adjust == 'b' or adjust == 'n':

            # adjust the means
            print('Computing and adjusting their mean levels')

            # now carry out the adjustments
            for ccd in ccds:
                if adjust == 'b':
                    ccd += mean - ccd.mean()
                elif adjust == 'n':
                    ccd *= mean / ccd.mean()

        # Finally, combine
        print('Combining them and storing the result')

        for wnam, wind in template[cnam].items():
            # build list of all data arrays
            arrs = [ccd[wnam].data for ccd in ccds]
            arr3d = np.stack(arrs)
            wind.data = np.median(arr3d,axis=0)

    # write out
    template.wfits(outfile, clobber)
    print('\nFinal result written to {:s}'.format(outfile))
