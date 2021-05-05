import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from trm.pgplot import *

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

###################################
#
# hplot -- plots a multi-CCD image.
#
###################################


def hplot(args=None):
    """``hplot input [device] ccd nx msub [cmap] hsbox iset (ilo ihi | plo
    phi) xlo xhi ylo yhi [width height]``

    Plots a multi-CCD image. Can use PGPLOT or matplotlib. The matplotlib
    version is slightly clunky in its choice of the viewing area but has some
    features that could be useful, in particular, the interactive plot
    (device='/mpl') allows one to pan and zoom and to compare the same part of
    multiple CCDs easily.

    Parameters:

      input : string
         name of MCCD file

      device : string [hidden]
         Plot device name. Uses characters after a final trailing '/' to
         identify the type in PGPLOT style. Thus:

          |  /xs : PGPLOT xserver interactive plot
          |  1/xs : PGPLOT xserver interactive plot called '1'
          |  plot.ps/cps : PGPLOT colour postscript called 'plot.ps'
          |  plot.ps/vps : PGPLOT B&W portrait oriented plot
          |  /mpl : matplotlib interactive plot
          |  plot.pdf/mpl : matplotlib PDF plot

      ccd : string
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      nx : int
         number of panels across to display, prompted if more than one CCD is
         to be plotted.

      msub : bool
         True/False to subtract median from each window before scaling

      cmap : str [if matplotlib]
         The colour map to use. "Greys" is the usual, but there are
         many others. Typing an incorrect one will give a list. "none"
         for matplotlib default.

      hsbox : int [if device = '/mpl'; hidden]
         half-width in binned pixels of stats box as offset from central pixel
         hsbox = 1 gives a 3x3 box; hsbox = 2 gives 5x5 etc.

      iset : string [single character]
         determines how the intensities are determined. There are three
         options: 'a' for automatic simply scales from the minimum to the
         maximum value found on a per CCD basis. 'd' for direct just takes two
         numbers from the user. 'p' for percentile dtermines levels based upon
         percentiles determined from the entire CCD on a per CCD bais.

      ilo : float [if iset=='d']
         lower intensity level

      ihi : float [if iset=='d']
         upper intensity level

      plo : float [if iset=='p']
         lower percentile level

      phi : float [if iset=='p']
         upper percentile level

      xlo : float
         left X-limit, for PGPLOT plots. Applies to matplotlib plots
         to restrict region used to compute percentile limits. This is
         useful in case where bias strips otherwise distort the plot
         limits (e.g. ultraspec full frame images)

      xhi : float
         right X-limit. See comments for xlo as well.

      ylo : float
         bottom Y-limit. See comments for xlo as well.

      yhi : float
         top Y-limit. See comments for xlo as well.

      width : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

    """

    global fig, mccd, caxes, hsbox

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("input", Cline.LOCAL, Cline.PROMPT)
        cl.register("device", Cline.LOCAL, Cline.HIDE)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("msub", Cline.GLOBAL, Cline.PROMPT)
        cl.register("cmap", Cline.GLOBAL, Cline.PROMPT)
        cl.register("hsbox", Cline.GLOBAL, Cline.HIDE)
        cl.register("iset", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ilo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ihi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("plo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("phi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("xlo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("xhi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ylo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("yhi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)

        # get inputs
        frame = cl.get_value("input", "frame to plot", cline.Fname("hcam", hcam.HCAM))
        mccd = hcam.MCCD.read(frame)

        device = cl.get_value("device", "plot device name", "/mpl")

        # set type of plot (PGPLOT or matplotlib) and the name of the file
        # if any in the case of matplotlib
        fslash = device.rfind("/")
        if fslash > -1:
            if device[fslash + 1 :] == "mpl":
                ptype = "MPL"
                hard = device[:fslash].strip()
            else:
                ptype = "PGP"

        else:
            raise ValueError(
                "Could not identify plot type from device = {:s}".format(device)
            )

        # define the panel grid
        try:
            nxdef = cl.get_default("nx")
        except:
            nxdef = 3

        max_ccd = len(mccd)
        if max_ccd > 1:
            ccd = cl.get_value("ccd", "CCD(s) to plot [0 for all]", "0")
            if ccd == "0":
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
            if len(ccds) > 1:
                nxdef = min(len(ccds), nxdef)
                cl.set_default("nx", nxdef)
                nx = cl.get_value("nx", "number of panels in X", 3, 1)
            else:
                nx = 1
        else:
            ccds = list(mccd.keys())
            nx = 1

        # define the display intensities
        msub = cl.get_value("msub", "subtract median from each window?", True)

        if ptype == "MPL":
            cmap = cl.get_value("cmap", "colour map to use ['none' for mpl default]", "Greys")
            cmap = None if cmap == "none" else cmap

        if ptype == "MPL" and hard == "":
            hsbox = cl.get_value(
                "hsbox", "half-width of stats box (binned pixels)", 2, 1
            )

        iset = cl.get_value(
            "iset",
            "set intensity a(utomatically), d(irectly) or with p(ercentiles)?",
            "a",
            lvals=["a", "A", "d", "D", "p", "P"],
        )
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == "d":
            ilo = cl.get_value("ilo", "lower intensity limit", 0.0)
            ihi = cl.get_value("ihi", "upper intensity limit", 1000.0)
        elif iset == "p":
            plo = cl.get_value(
                "plo", "lower intensity limit percentile", 5.0, 0.0, 100.0
            )
            phi = cl.get_value(
                "phi", "upper intensity limit percentile", 95.0, 0.0, 100.0
            )

        # region to plot
        for i, cnam in enumerate(ccds):
            ccd = mccd[cnam]
            nxtot, nytot, nxpad, nypad = ccd.nxtot, ccd.nytot, ccd.nxpad, ccd.nypad
            if i == 0:
                xmin, xmax = float(-nxpad), float(nxtot + nxpad + 1)
                ymin, ymax = float(-nypad), float(nytot + nypad + 1)
            else:
                xmin = min(xmin, float(-nxpad))
                xmax = max(xmax, float(nxtot + nxpad + 1))
                ymin = min(ymin, float(-nypad))
                ymax = max(ymax, float(nytot + nypad + 1))

        xlo = cl.get_value("xlo", "left-hand X value", xmin, xmin, xmax)
        xhi = cl.get_value("xhi", "right-hand X value", xmax, xmin, xmax)
        ylo = cl.get_value("ylo", "lower Y value", ymin, ymin, ymax)
        yhi = cl.get_value("yhi", "upper Y value", ymax, ymin, ymax)
        width = cl.get_value("width", "plot width (inches)", 0.0)
        height = cl.get_value("height", "plot height (inches)", 0.0)

    # all inputs obtained, plot
    if ptype == "MPL":
        if width > 0 and height > 0:
            fig = plt.figure(figsize=(width, height))
        else:
            fig = plt.figure()

        mpl.rcParams["xtick.labelsize"] = hcam.mpl.Params["axis.number.fs"]
        mpl.rcParams["ytick.labelsize"] = hcam.mpl.Params["axis.number.fs"]

        nccd = len(ccds)
        ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

        ax = None
        caxes = {}
        for n, cnam in enumerate(ccds):
            if ax is None:
                axes = ax = fig.add_subplot(ny, nx, n + 1)
                axes.set_aspect("equal", adjustable="box")
            else:
                axes = fig.add_subplot(ny, nx, n + 1, sharex=ax, sharey=ax)
                axes.set_aspect("equal")

            # store the CCD associated with these axes for the cursor callback
            caxes[axes] = cnam

            axes.set_xlim(xlo, xhi)
            axes.set_ylim(ylo, yhi)

            if msub:
                # subtract median from each window
                for wind in mccd[cnam].values():
                    wind -= wind.median()

            vmin, vmax = hcam.mpl.pCcd(
                axes,
                mccd[cnam],
                iset,
                plo,
                phi,
                ilo,
                ihi,
                "CCD {:s}".format(cnam),
                xlo,
                xhi,
                ylo,
                yhi, cmap=cmap
            )
            print("CCD =", cnam, "plot range =", vmin, "to", vmax)

        try:
            plt.tight_layout()
        except:
            pass

        if hard == "":

            # add in the callback
            fig.canvas.mpl_connect("button_press_event", buttonPressEvent)
            print(
                "\nClick points in windows for stats in a {:d}x{:d} box".format(
                    2 * hsbox + 1, 2 * hsbox + 1
                )
            )
            plt.subplots_adjust(wspace=0.1, hspace=0.1)
            plt.show()
        else:
            plt.savefig(hard, bbox_inches="tight", pad_inches=0)

    elif ptype == "PGP":
        # open the plot
        dev = hcam.pgp.Device(device)
        if width > 0 and height > 0:
            pgpap(width, height / width)

        nccd = len(ccds)
        ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

        # set up panels and axes
        pgsubp(nx, ny)

        for cnam in ccds:
            pgsci(hcam.pgp.Params["axis.ci"])
            pgsch(hcam.pgp.Params["axis.number.ch"])
            pgenv(xlo, xhi, ylo, yhi, 1, 0)

            vmin, vmax = hcam.pgp.pCcd(
                mccd[cnam], iset, plo, phi, ilo, ihi, "CCD {:s}".format(cnam)
            )
            print("CCD =", cnam, "plot range =", vmin, "to", vmax)


def buttonPressEvent(event):
    """
    callback 
    """
    global fig, mccd, caxes, hsbox

    pzoom = fig.canvas.manager.toolbar.mode == "pan/zoom"

    if not pzoom:
        # only when not in pan/zoom mode
        if event.inaxes is not None:
            cnam = caxes[event.inaxes]
            ccd = mccd[cnam]
            x, y = event.xdata, event.ydata
            utils.print_stats(ccd, cnam, x, y, hsbox)
