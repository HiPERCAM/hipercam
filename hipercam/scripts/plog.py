import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

##########################################
#
# plog -- simple plots of reduce log files
#
##########################################


def plog(args=None):
    """``plog log [device width height] ccd1 aper1 param1 ccd2 (aper2 param2
    scheme)``

    Provides quick-look plots of HiPERCAM |reduce| logs.

    Parameters:

      log : string
          name of |reduce| ASCII log file (text file with loads of columns)

      device : string [hidden, defaults to 'term']
         'term' for interactive plot, file name such as 'plot.pdf'
         for a hardcopy.

      width : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH
         width AND height must be non-zero to have any effect

      ccd1 : string
         first CCD to consider, e.g. '1'

      aper1 : string
         first aperture to consider

      param1 : string
         first parameter to consider. Choices are 'x' = X position,
         'y' = Y position, 'f' = FWHM, 'b' = Moffat beta, 's' = sky.

      ccd2 : string
         second CCD to consider; '!' to ignore. Can be (and typically
         would be) the same as ccd1.

      aper2 : string [if ccd2 != '!']
         second aperture to consider

      param2 : string [if ccd2 != '!']
         second parameter. See param1 for choices

      scheme : string [if ccd2 != '!']
         how to plot if both apertures are chosen. Choices:

            | 'd' = difference, i.e. plot 1-2
            | 'b' = both plotted on same panel
            | 'r' = ratio, i.e. 1 / 2, good for relative photom
            | 's' = scatter plot, 2 vs 1.

    .. Note::

       Points with negative errors are ignored. Be careful with linked
       apertures where all x, y, FWHM, beta automatically have negative errors
       since they are not fitted.

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("log", Cline.LOCAL, Cline.PROMPT)
        cl.register("device", Cline.LOCAL, Cline.HIDE)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)
        cl.register("ccd1", Cline.LOCAL, Cline.PROMPT)
        cl.register("aper1", Cline.LOCAL, Cline.PROMPT)
        cl.register("param1", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd2", Cline.LOCAL, Cline.PROMPT)
        cl.register("aper2", Cline.LOCAL, Cline.PROMPT)
        cl.register("param2", Cline.LOCAL, Cline.PROMPT)
        cl.register("scheme", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        log = cl.get_value(
            "log", "reduce log file to plot", cline.Fname("hcam", hcam.LOG)
        )

        device = cl.get_value("device", "plot device name", "term")
        width = cl.get_value("width", "plot width (inches)", 0.0)
        height = cl.get_value("height", "plot height (inches)", 0.0)

        ccd1 = cl.get_value("ccd1", "first CCD to plot", "1")
        aper1 = cl.get_value("aper1", "first aperture", "1")
        param1 = cl.get_value(
            "param1",
            "first parameter [x, y, f(whm), b(eta), c(ounts), s(ky)]",
            "c",
            lvals=("x", "y", "f", "b", "c", "s"),
        )
        lab1 = "[CCD={:s},ap={:s},p={:s}]".format(ccd1, aper1, param1)

        ccd2 = cl.get_value("ccd2", "second CCD to plot ['!' to ignore]", ccd1)
        if ccd2 != "!":
            aper2 = cl.get_value("aper2", "second aperture", "1")
            param2 = cl.get_value(
                "param2",
                "second parameter [x, y, f(whm), b(eta), c(ounts), s(ky)]",
                "c",
                lvals=("x", "y", "f", "b", "c", "s"),
            )
            lab2 = "[CCD={:s},ap={:s},p={:s}]".format(ccd2, aper2, param2)

        # map to names used in reduce log files
        MAP = {"x": "x", "y": "y", "f": "fwhm", "b": "beta", "c": "counts", "s": "sky"}
        pname1 = MAP[param1]
        if ccd2 != "!":
            pname2 = MAP[param2]
            scheme = cl.get_value(
                "scheme",
                "b(oth), d(ifference), r(atio), s(catter)",
                "b",
                lvals=("b", "d", "r", "s"),
            )

    # load the reduce log
    hlog = hcam.hlog.Hlog.read(log)

    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width, height))
    else:
        fig = plt.figure()

    dat1 = hlog.tseries(ccd1, aper1, pname1)

    if ccd2 != "!":
        dat2 = hlog.tseries(ccd2, aper2, pname2)
        if scheme == "b":
            # plots both together
            dat1.mplot(plt, "b", bitmask=hcam.BAD_TIME)
            dat2.mplot(plt, "r", bitmask=hcam.BAD_TIME)
            xlabel = "Time [MJD]"
            ylabel = "{:s} & {:s}".format(lab1, lab2)

        elif scheme == "r":
            # ratio
            ratio = dat1 / dat2
            ratio.mplot(plt, "b", bitmask=hcam.BAD_TIME)
            xlabel = "Time [MJD]"
            ylabel = "{:s} / {:s}".format(lab1, lab2)

        elif scheme == "d":
            # difference
            diff = dat1 - dat2
            diff.mplot(plt, "b", bitmask=hcam.BAD_TIME)
            xlabel = "Time [MJD]"
            ylabel = "{:s} - {:s}".format(lab1, lab2)

        elif scheme == "s":
            mask1 = dat1.get_mask(bitmask=hcam.BAD_TIME, invert=True)
            mask2 = dat1.get_mask(bitmask=hcam.BAD_TIME, invert=True)
            ok = ~mask1 & ~mask2
            plt.errorbar(
                dat1.y[ok], dat2.y[ok], dat2.ye[ok], dat1.ye[ok], ".", capsize=0
            )
            xlabel = lab1
            ylabel = lab2

    else:
        # just one
        dat1.mplot(plt, bitmask=hcam.BAD_TIME)
        xlabel = "Time [MJD]"
        ylabel = lab1

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(log)

    if device == "term":
        plt.show()
    else:
        plt.savefig(device)
