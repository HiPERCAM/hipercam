import matplotlib as mpl
import matplotlib.pyplot as plt
from trm import cline
from trm.cline import Cline

import hipercam as hcam
from hipercam import utils
from hipercam.scripts.nrtplot import Pfit

###################################
#
# hplot -- plots a multi-CCD image.
#
###################################


def hplot(args=None):
    """``hplot input [device] ccd nx msub ([cmap]) hsbox iset (ilo ihi | plo
    phi) xlo xhi ylo yhi [width height]``

    Plots a multi-CCD image using matplotlib. The default matplotlib interactive
    window allows one to pan and zoom and to compare the same part of
    multiple CCDs easily.

    Pressing the space bar will initiate a profile fit at the position of the cursor.
    The parameters of a successful fit will be printed to the screen.
    A radial profile plot centered on the fit position can be plotted
    by pressing the 'r' key. This will open a new matplotlib window.

    Parameters:

      input : string
         name of MCCD file

      device : string [hidden, defaults to 'term']
         'term' for interactive plot, file name such as 'plot.pdf'
         or 'object.png' for a hardcopy.

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

      cmap : str [hidden]
         The colour map to use. "Greys" is the usual, but there are
         many others. Typing an incorrect one will give a list. "none"
         for matplotlib default.

      hsbox : int [hidden]
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
         left X-limit to restrict region used to compute percentile
         limits. This allows you to only consider useful data if for
         example there are bias strips, and makes setting the percentile
         limits less fiddly.

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
          when displaying an interactive window, pressing the space bar will initiate
          a profile fit at the position of the cursor.
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

        cl.set_default("device", "term")
        device = cl.get_value("device", "plot device name", "term")
        if device.endswith(".png"):
            dpi = cl.get_value("dpi", "dots per inch", 240, 10)

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

        cmap = cl.get_value(
            "cmap", "colour map to use ['none' for mpl default]", "Greys"
        )
        cmap = None if cmap == "none" else cmap

        if device == "term":
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
            "half width of box for initial location of target [unbinned pixels]",
            11.0,
            2.0,
        )
        smooth = cl.get_value(
            "smooth",
            "FWHM for smoothing for initial object detection [binned pixels]",
            6.0,
        )
        fhbox = cl.get_value(
            "fhbox",
            "half width of box for profile fit [unbinned pixels]",
            21.0,
            3.0,
        )
        hmin = cl.get_value("hmin", "minimum peak height to accept the fit", 20.0)
        read = cl.get_value("read", "readout noise, RMS ADU", -1.0)
        gain = cl.get_value("gain", "gain, ADU/e-", 1.0)
        thresh = cl.get_value("thresh", "number of RMS to reject at", 4.0)

    # create the plot
    fig, caxes = hcam.mpl.plot_frame(
        mccd=mccd,
        cnams=ccds,
        ncols=nx,
        msub=msub,
        cmap=cmap,
        iset=iset,
        ilim=(plo, phi) if iset == "p" else (ilo, ihi),
        region=(xlo, xhi, ylo, yhi),
        title=frame,
        figsize=(width, height),
    )

    # print the plot ranges for each CCD
    for cnam in caxes:
        imgs = caxes[cnam].get_images()
        vmin, vmax = imgs[0].get_clim()  # the first window, but should all be the same
        print(f"CCD = {cnam}, plot range = {vmin} to {vmax}")

    if device != "term":
        # save plot to file
        if device.endswith(".png"):
            plt.savefig(device, dpi=dpi, bbox_inches="tight", pad_inches=0)
        else:
            plt.savefig(device, bbox_inches="tight", pad_inches=0)
        return

    # print out instructions for the user
    box_size = 2 * hsbox + 1
    print(
        f"\nClick points in windows for stats in a {box_size:d}x{box_size:d} box"
        f"\nPress space bar to fit a profile at the cursor position"
        f"\nPress 'r' to plot a radial profile of the fit in a new window"
    )

    # disable 'r' key as shortcut for reset zoom
    mpl.rcParams["keymap.home"].remove("r")

    # create the picker class and connect it to the figure
    picker = PickProfit(
        mccd,
        fig,
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
        hsbox,
    )
    fig.canvas.mpl_connect("key_press_event", picker._keyPressEvent)
    fig.canvas.mpl_connect("button_press_event", picker._buttonPressEvent)

    # finally set up the plot for interactive mode and display
    hcam.mpl.setup_plot(fig, disable_shortcuts=False)
    plt.show()


class PickProfit:
    """Class to pick targets for profile fitting."""

    def __init__(
        self,
        mccd,
        fig,
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
        hsbox=2,
    ):
        # image content data
        self.mccd = mccd

        # figure and axes
        self.fig = fig
        self.caxes = caxes
        self.axnams = {axes: cnam for cnam, axes in caxes.items()}  # reverse for events

        # fitting params
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

        # other params
        self.hsbox = hsbox

        # list to hold future plot artists
        self.artists = None

    def _keyPressEvent(self, event):
        """Fit a profile at the position of the cursor when the key is pressed."""
        # matplotlib 3.11 changes callback order so for "q" there is no
        # canvas manager any more
        if event.key == "q":
            return

        # only when curser is within the axes
        if event.inaxes is None:
            return

        # valid keys are space bar for fitting and 'r' for radial profile plot
        if not event.key in [" ", "r"]:
            return

        # check that the position is inside a window
        cnam = self.axnams[event.inaxes]
        ccd = self.mccd[cnam]
        x, y = event.xdata, event.ydata
        wnam = ccd.inside(x, y, 2)
        if wnam is None:
            return

        # profile fit at the position of the cursor
        pfit = Pfit(cnam, wnam, x, y, self.shbox, self.fwhm, self.beta)
        results, message = pfit.fit(
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
        if self.artists is not None:
            for artist in self.artists:
                artist.remove()
        self.artists = pfit.plot(event.inaxes)
        self.fig.canvas.draw_idle()

        print(
            f"\nProfile fit at x,y = {x:.2f},{y:.2f} in CCD {cnam}, window {wnam}"
            " [applies to previous frame]\n",
            message,
        )

        if results is None:
            # fit failed
            self.artists[0].set_color("red")
            self.artists[1].set_color("red")
            self.fig.canvas.draw_idle()
            return

        # fitted OK
        if event.key == "r":
            # plot the radial fit
            (
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

    def _buttonPressEvent(self, event):
        """Print statistics when the window is clicked."""
        # only when curser is within the axes and not in pan/zoom mode
        if self.fig.canvas.manager.toolbar.mode == "pan/zoom":
            return
        if event.inaxes is None:
            return

        cnam = self.axnams[event.inaxes]
        utils.print_stats(
            self.mccd[cnam], cnam, event.xdata, event.ydata, self.hsbox
        )
