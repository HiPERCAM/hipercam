import sys
import os

import numpy as np
import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils, support, spooler
from hipercam.cline import Cline

__all__ = ['combine',]

############################
#
# combine -- combines images
#
############################

def combine(args=None):
    """``combine list bias dark flat method (sigma) adjust (usemean) [plot clobber]
    output``

    Combines a series of images defined by a list using median or clipped
    mean combination. Only combines those CCDs for which is_data() is true
    (i.e. it skips blank frames caused by NSKIP / NBLUE options)

    Parameters:

        list : string
           list of hcm files with images to combine. The formats of the
           images should all match

        bias : string
           Name of bias frame to subtract, 'none' to ignore.

        dark : string
           Name of dark frame to subtract, 'none' to ignore.

        flat : string
           Name of flat field frame to subtract, 'none' to ignore.

        method  : string
           'm' for median, 'c' for clipped mean. See below for pros and cons.

        sigma   : float [if method == 'c']
           With clipped mean combination, pixels that deviate by more than
           sigma RMS from the mean are kicked out. This is carried out in an
           iterative manner. sigma <= 0 implies no rejection, just a straight
           average. sigma=3 is typical.

        adjust  : string
           adjustments to make: 'i' = ignore; 'n' = normalise the mean or
           median of all frames to match the first; 'b' = add offsets so that
           the mean or median of all frames is the same as the first. Option
           'n' is useful for twilight flats and fringe frames; 'b' is good
           for combining biases.

        usemean : bool [if adjust == 'n' or 'b']
           True to use the mean rather than the median for normalisation or
           biass offsetting.

        plot    : bool [hidden, if adjust == 'n' or 'b'; defaults to False]
           make a plot of the mean versus frame number. This can provide a
           quick check that the frames are not too different.

        clobber : bool [hidden]
           clobber any pre-existing output files

        output  : string
           output file

       Clipped mean can work well for large numbers of frames but gets worse
       for small numbers as the RMS can be heavily influenced by a single bad
       value. The median can be better in such cases, but has the downside of
       digitisation noise. For instance, the average of 100 bias frames could
       have a noise level significantly below 1 count, depending upon the
       readout noise, and the +/- 0.5 count uncertainty of median combination
       may be worse than this.

    .. Note::

       This routine reads all inputs into memory, so can be a bit of a
       hog. However, it does so one CCD at a time to alleviate this. It will
       fail if it cannot find a valid frame for any CCD.

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('list', Cline.GLOBAL, Cline.PROMPT)
        cl.register('bias', Cline.LOCAL, Cline.PROMPT)
        cl.register('dark', Cline.LOCAL, Cline.PROMPT)
        cl.register('flat', Cline.LOCAL, Cline.PROMPT)
        cl.register('method', Cline.LOCAL, Cline.PROMPT)
        cl.register('sigma', Cline.LOCAL, Cline.PROMPT)
        cl.register('adjust', Cline.LOCAL, Cline.PROMPT)
        cl.register('usemean', Cline.LOCAL, Cline.PROMPT)
        cl.register('plot', Cline.LOCAL, Cline.HIDE)
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
            bias = hcam.MCCD.read(bias)

        # dark frame (if any)
        dark = cl.get_value(
            'dark', "dark frame ['none' to ignore]",
            cline.Fname('dark', hcam.HCAM), ignore='none'
        )
        if dark is not None:
            # read the dark frame
            dark = hcam.MCCD.read(dark)

        # flat frame (if any)
        flat = cl.get_value(
            'flat', "flat frame ['none' to ignore]",
            cline.Fname('flat', hcam.HCAM), ignore='none'
        )
        if flat is not None:
            # read the flat frame
            flat = hcam.MCCD.read(flat)

        method = cl.get_value(
            'method', 'c(lipped mean), m(edian)', 'c', lvals=('c','m')
        )

        if method == 'c':
            sigma = cl.get_value(
                'sigma', 'number of RMS deviations to clip', 3.
            )

        adjust = cl.get_value(
            'adjust', 'i(gnore), n(ormalise), b(ias offsets)',
            'i', lvals=('i','n','b')
        )

        if adjust == 'n' or adjust == 'b':

            usemean = cl.get_value(
                'usemean',
                'use the mean for normalisation / offsetting [else median]',
                True
            )

            plot = cl.get_value(
                'plot',
                'plot mean levels versus frame number?' if usemean else \
                'plot median levels versus frame number?', False
            )

        else:
            plot = False
            usemean = False

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

    template = hcam.MCCD.read(
        utils.add_extension(template_name,hcam.HCAM)
    )

    if bias is not None:
        # crop the bias
        bias = bias.crop(template)

    if dark is not None:
        # crop the dark
        dark = dark.crop(template)

    if flat is not None:
        # crop the flat
        flat = flat.crop(template)

    # Now process each file CCD by CCD to reduce the memory
    # footprint
    for cnam in template:

        # Read files into memory, insisting that they
        # all have the same set of CCDs
        print(
            "\nLoading all CCDs labelled '{:s}' from {:s}".format(cnam, flist)
        )

        ccds, means = [], []
        nrej, ntot = 0, 0
        with spooler.HcamListSpool(flist, cnam) as spool:

            if bias is not None:
                # extract relevant CCD from the bias
                bccd = bias[cnam]
                bexpose = bias.head.get('EXPTIME',0.)
            else:
                bexpose = 0.

            if dark is not None:
                # extract relevant CCD from the dark
                dccd = dark[cnam]
                dexpose = dark.head['EXPTIME']

            if flat is not None:
                # extract relevant CCD from the flat
                fccd = flat[cnam]

            mean = None
            for ccd in spool:

                if ccd.is_data():

                    if bias is not None:
                        # subtract bias
                        ccd -= bccd

                    if dark is not None:
                        # subtract dark
                        scale = (ccd.head['EXPTIME']-bexpose)/dexpose
                        ccd -= scale*dccd

                    if flat is not None:
                        # apply flat
                        ccd /= fccd

                    # keep the result
                    ccds.append(ccd)

                    if (adjust == 'b' or adjust == 'n') and mean is None:
                        # store the first mean [median]
                        if usemean:
                            mean = ccd.mean()
                        else:
                            mean = ccd.median()

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
                if usemean:
                    cmean = ccd.mean()
                else:
                    cmean = ccd.median()

                means.append(cmean)
                if adjust == 'b':
                    ccd += mean - cmean
                elif adjust == 'n':
                    ccd *= mean / cmean

            if plot:
                plt.plot(means)
                plt.plot(means,'.k')
                plt.text(len(means)+1,means[-1],cnam,va='center',ha='left')

        # Finally, combine
        if method == 'm':
            print('Combining them (median) and storing the result')
        elif method == 'c':
            print('Combining them (clipped mean) and storing the result')
        else:
            raise NotImplementedError('method = {:s} not implemented'.format(method))

        for wnam, wind in template[cnam].items():

            # build list of all data arrays
            arrs = [ccd[wnam].data for ccd in ccds]
            arr3d = np.stack(arrs)

            # at this point, arr3d is a 3D array, with the first dimension
            # (axis=0) running over the images. We want to average / median
            # over this axis.

            if method == 'm':
                # median
                wind.data = np.median(arr3d,axis=0)

            elif method == 'c':
                # Cython routine avgstd requires np.float32 input
                arr3d = arr3d.astype(np.float32)
                if sigma > 0.:
                    avg, std, num = support.avgstd(arr3d, sigma)
                    nrej += len(ccds)*num.size-num.sum()
                    ntot += len(ccds)*num.size
                else:
                    avg = np.mean(arr3d,axis=0)
                wind.data = avg

        # Add history
        if method == 'm':
            template[cnam].head.add_history(
                'Median combine of {:d} images'.format(len(ccds))
                )
        elif method == 'c':
            print('Rejected {:d} pixels = {:.3f}% of the total'.format(nrej,100*nrej/ntot))
            template[cnam].head.add_history(
                'Clipped mean combine of {:d} images, sigma = {:.1f}'.format(len(ccds), sigma)
                )

    # write out
    template.write(outfile, clobber)
    print('\nFinal result written to {:s}'.format(outfile))

    # finally
    if plot:
        plt.xlabel('Frame number')
        plt.ylabel('Mean counts')
        plt.show()
