import sys
import os

import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

#######################################
#
# hfilter -- filters a multi-CCD image.
#
#######################################


def hfilter(args=None):
    """``hfilter input ccd nx filter [fwhm] output``

    Filters a multi-CCD image. e.g. smooths it. Can be useful
    in some analysis steps, e.g. for picking out defects, division
    by a smoother version of an image can be useful.

    Parameters:

      input  : string
         name of MCCD input file

      ccd    : string
         CCD(s) to filter, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes).

      filter  : string [single character]
         type of filter. 'g' = gaussian. Uses an FFT-based approach which
         regards the boundaries as periodic, so you will get significant
         edge-effects if the values on opposite edges of a window are
         significantly different from each other.

      clobber : bool [hidden]
         clobber any pre-existing output files

      output  : string
         name of MCCD output file

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("input", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("filter", Cline.LOCAL, Cline.PROMPT)
        cl.register("fwhm", Cline.LOCAL, Cline.PROMPT)
        cl.register("clobber", Cline.LOCAL, Cline.HIDE)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        frame = cl.get_value("input", "frame to filter", cline.Fname("hcam", hcam.HCAM))
        mccd = hcam.MCCD.read(frame)

        max_ccd = len(mccd)
        if max_ccd > 1:
            ccd = cl.get_value("ccd", "CCD(s) to filter [0 for all]", "0")
            if ccd == "0":
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(mccd.keys())

        # define the display intensities
        filt = cl.get_value("filter", "filter to apply: g(aussian)", "g", lvals=["g",])

        if filt == "g":
            fwhm = cl.get_value(
                "fwhm", "FWHM for gaussian filter (binned pixels)", 4.0, 0.01
            )

        clobber = cl.get_value(
            "clobber", "clobber any pre-existing files on output", False
        )

        output = cl.get_value(
            "output",
            "output file",
            cline.Fname(
                "filtered",
                hcam.HCAM,
                cline.Fname.NEW if clobber else cline.Fname.NOCLOBBER,
            ),
        )

    # all inputs obtained, filter

    if filt == "g":
        # create filter kernal
        kern = Gaussian2DKernel(fwhm / np.sqrt(8 * np.log(2)))

    for cnam in ccds:

        ccd = mccd[cnam]
        for wnam, wind in ccd.items():
            if filt == "g":
                wind.data = convolve_fft(wind.data, kern, "wrap")

        print("Filtered CCD {:s}".format(cnam))

    # write out
    mccd.write(output, clobber)
    print("\nFiltered result written to {:s}".format(output))
