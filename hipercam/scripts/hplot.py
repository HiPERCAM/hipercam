import sys
import os
from functools import partial

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from trm.pgplot import *
from trm import cline
from trm.cline import Cline

import hipercam as hcam
from hipercam import utils
from hipercam.scripts.nrtplot import Fpar

###################################
#
# hplot -- plots a multi-CCD image.
#
###################################


def hplot(args=None):
    """``hplot input [device] ccd nx msub ([cmap]) hsbox iset (ilo ihi | plo
    phi) xlo xhi ylo yhi [width height]``

    Plots a multi-CCD image. Can use PGPLOT or matplotlib. The matplotlib
    version is slightly clunky in its choice of the viewing area but has some
    features that could be useful, in particular, the interactive plot
    (device='/mpl') allows one to pan and zoom and to compare the same part of
    multiple CCDs easily.

    If the interactive plot is used (device='/mpl'), then pressing the space
    bar will initiate a profile fit at the position of the cursor. The
    parameters of a successfull fit will be printed to the screen.
    A radial profile plot centered on the fit position can be plotted
    by pressing the 'r' key. This will open a new matplotlib window.

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

      cmap : str [if matplotlib; hidden]
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

      method : str [hidden]
          if the device is set to '/mpl' (matplotlib), then pressing the
          space bar will initiate a profile fit at the position of the cursor.
          This defines the profile fitting method, either a gaussian or a
          moffat profile. The latter is usually best.

      beta : float [method == 'm'; hidden]
          default Moffat exponent

      fwhm : float [hidden]
          default FWHM, unbinned pixels. Do try to get this about right as
          it affects whether the profile can be fitted at all.

      fwhm_min : float [hidden]
          minimum FWHM to allow, unbinned pixels.

      shbox : float [hidden]
          half width of box for searching for a star, unbinned
          pixels. The above-threshold target closest to the centre of
          the box in a region +/- shbox around an intial position
          will be selected. It may not be the brightest, depending
          upon your threshold settings, so use those to filter faint
          objects. If profit=True, 'shbox' should be large enough to
          allow for likely changes in position from frame to frame,
          but not too large to avoid jumping to brighter targets or
          possibly cosmic rays.

      smooth : float [hidden]
          FWHM for gaussian smoothing, binned pixels. The initial position
          for fitting is determined by finding the maximum flux in a smoothed
          version of the image in a box of width +/- shbox around the starter
          position. Typically should be comparable to the stellar width. Its
          main purpose is to combat cosmic rays which tend only to occupy a
          single pixel.

      fhbox : float [hidden]
          half width of box for profile fit, unbinned pixels. The fit box is
          centred on the position located by the initial search. It should
          normally be > ~2x the expected FWHM, and usually smaller than shbox

      hmin : float [hidden]
          height threshold to accept a fit. If the height is below this
          value, the position will not be updated. This is to help in cloudy
          conditions. The limit is applied to the image after it has been
          smoothed to make less vulnerable to seeing fluctuations. This
          can mean it can be quite small.

      read : float [hidden]
          readout noise, RMS ADU, for assigning uncertainties. Set to a -ve
          value to try to ascertain on the fly; this is advisable either if
          you don't have a bias or you apply 'msub', and it will default to
          this if not specified in this case. The value returned in this case
          includes sky noise, i.e it should be roughly sqrt(R**2+S/G) where
          R is the true read noise, S are the sky counts per pixel, and G the
          gain. If read is set -ve, two fits are carried out per target. The
          second of these should usually be pretty fast. The first is carried
          out with an assumed large read noise of 20 in order to soften the
          weights. The results reported apply to the second fit.

      gain : float [hidden]
          gain, ADU/count, for assigning uncertainties.

      thresh : float [hidden]
          sigma rejection threshold for fits

    """

    global fig, mccd, caxes, hsbox

    command, args = cline.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:
        # register parameters
        cl.register("input", Cline.LOCAL, Cline.PROMPT)
        cl.register("device", Cline.LOCAL, Cline.HIDE)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("msub", Cline.GLOBAL, Cline.PROMPT)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
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
        cl.register("method", Cline.LOCAL, Cline.HIDE)
        cl.register("beta", Cline.LOCAL, Cline.HIDE)
        cl.register("fwhm", Cline.LOCAL, Cline.HIDE)
        cl.register("fwhm_min", Cline.LOCAL, Cline.HIDE)
        cl.register("shbox", Cline.LOCAL, Cline.HIDE)
        cl.register("smooth", Cline.LOCAL, Cline.HIDE)
        cl.register("fhbox", Cline.LOCAL, Cline.HIDE)
        cl.register("hmin", Cline.LOCAL, Cline.HIDE)
        cl.register("read", Cline.LOCAL, Cline.HIDE)
        cl.register("gain", Cline.LOCAL, Cline.HIDE)
        cl.register("thresh", Cline.LOCAL, Cline.HIDE)

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
        nxdef = cl.get_default("nx", 3)

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
            cmap = cl.get_value(
                "cmap", "colour map to use ['none' for mpl default]", "Greys"
            )
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

        xlo = cl.get_value("xlo", "left-hand X value", xmin, xmin, xmax, enforce=False)
        xhi = cl.get_value("xhi", "right-hand X value", xmax, xmin, xmax, enforce=False)
        ylo = cl.get_value("ylo", "lower Y value", ymin, ymin, ymax, enforce=False)
        yhi = cl.get_value("yhi", "upper Y value", ymax, ymin, ymax, enforce=False)
        width = cl.get_value("width", "plot width (inches)", 0.0)
        height = cl.get_value("height", "plot height (inches)", 0.0)

        # many parameters for profile fits, although none are
        # requested by default
        method = cl.get_value(
            "method", "fit method g(aussian) or m(offat)", "m", lvals=["g", "m"]
        )
        if method == "m":
            beta = cl.get_value(
                "beta", "initial exponent for Moffat fits", 5.0, 0.5, 20.0
            )
        else:
            beta = 0.0

        fwhm_min = cl.get_value(
            "fwhm_min", "minimum FWHM to allow [unbinned pixels]", 1.5, 0.01
        )
        fwhm = cl.get_value(
            "fwhm",
            "initial FWHM [unbinned pixels] for profile fits",
            6.0,
            fwhm_min,
        )
        shbox = cl.get_value(
            "shbox",
            "half width of box for initial location" " of target [unbinned pixels]",
            11.0,
            2.0,
        )
        smooth = cl.get_value(
            "smooth",
            "FWHM for smoothing for initial object" " detection [binned pixels]",
            6.0,
        )
        fhbox = cl.get_value(
            "fhbox",
            "half width of box for profile fit" " [unbinned pixels]",
            21.0,
            3.0,
        )
        hmin = cl.get_value("hmin", "minimum peak height to accept the fit", 20.0)
        read = cl.get_value("read", "readout noise, RMS ADU", -1.0)
        gain = cl.get_value("gain", "gain, ADU/e-", 1.0)
        thresh = cl.get_value("thresh", "number of RMS to reject at", 4.0)

    # all inputs obtained, plot
    if ptype == "MPL":
        if width > 0 and height > 0:
            fig = plt.figure(figsize=(width, height))
        else:
            fig = plt.figure()

        # set plot label sizes
        mpl.rcParams["xtick.labelsize"] = hcam.mpl.Params["axis.number.fs"]
        mpl.rcParams["ytick.labelsize"] = hcam.mpl.Params["axis.number.fs"]
        # disable 'r' key as shortcut for reset zoom
        mpl.rcParams["keymap.home"].remove("r")

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

            vmin, vmax, _ = hcam.mpl.pCcd(
                axes,
                mccd[cnam],
                iset,
                plo,
                phi,
                ilo,
                ihi,
                "CCD {:s}".format(cnam),
                xlo=xlo,
                xhi=xhi,
                ylo=ylo,
                yhi=yhi,
                cmap=cmap,
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

            # add the keypress callback
            fitter = OnDemandFit(
                fig,
                mccd,
                caxes,
                shbox,
                fwhm,
                beta,
                method,
                smooth,
                fhbox,
                hmin,
                fwhm_min,
                read,
                gain,
                thresh,
            )
            fig.canvas.mpl_connect("key_press_event", fitter._keyPressEvent)

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


class OnDemandFit(object):
    """
    Callback
    """

    def __init__(
        self,
        fig,
        mccd,
        caxes,
        shbox=31,
        fwhm=6,
        beta=5,
        method="g",
        smooth=6,
        fhbox=21,
        hmin=50,
        fwhm_min=1.5,
        read=-1.0,
        gain=1.0,
        thresh=4.0,
    ):
        self.fig = fig
        self.mccd = mccd
        self.caxes = caxes
        self.shbox = shbox
        self.fwhm = fwhm
        self.beta = beta
        self.method = method
        self.smooth = smooth
        self.fhbox = fhbox
        self.hmin = hmin
        self.fwhm_min = fwhm_min
        self.read = read
        self.gain = gain
        self.thresh = thresh

    def _keyPressEvent(self, event):
        pzoom = self.fig.canvas.manager.toolbar.mode == "pan/zoom"
        # only when not in pan/zoom mode
        if not pzoom and event.key in [" ", "r"] and event.inaxes is not None:
            cnam = self.caxes[event.inaxes]
            ccd = self.mccd[cnam]
            x, y = event.xdata, event.ydata
            # check that the position is inside a window
            wnam = ccd.inside(x, y, 2)
            if wnam is not None:
                # store the position, Window label, target number,
                # box size fwhm, beta
                fpar = Fpar(cnam, wnam, x, y, self.shbox, self.fwhm, self.beta)
                results, message = fpar.fit(
                    ccd,
                    self.method,
                    self.smooth,
                    self.fhbox,
                    self.hmin,
                    self.fwhm_min,
                    self.read,
                    self.gain,
                    self.thresh,
                )
                if results is not None:
                    # fitted OK
                    print(
                        f"profile fit with initial x,y = "
                        f"{x:.2f}, {y:.2f} in CCD {cnam}, window {wnam}, [applies to previous frame]:"
                    )
                    print(f"   {message}\n")
                    if event.key == "r":
                        # plot the radial fit
                        (
                            _,
                            radius,
                            pixel_value,
                            _,
                            vmin,
                            vmax,
                            radius_fit,
                            fit,
                        ) = results

                        # needed to avoid error messages about event loop
                        plt.ion()

                        # open new figure (if not already open)
                        radial_fig = plt.figure("radial profile")
                        # clear in case already open
                        radial_fig.clf()

                        # make plot
                        ax = radial_fig.add_subplot(111)
                        ax.plot(radius, pixel_value, "k.")
                        ax.plot(radius_fit, fit, "k-")
                        ax.set_ylim(vmin, vmax)
                        ax.set_xlabel("Radius (pixels)")
                        ax.set_ylabel("Pixel value")

                else:
                    print(
                        f"\n   ** fit failed at position x,y = {x:.2f}, {y:.2f} in CCD {cnam}, window {wnam}"
                    )
                    print(f"   ** fit message = {message}\n")


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
