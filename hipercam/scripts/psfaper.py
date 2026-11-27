import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import SigmaClip
from astropy.table import Table

try:
    from photutils.background import MADStdBackgroundRMS, MMMBackground
    from photutils.detection import DAOStarFinder
    from photutils.psf import (
        GaussianPRF,
        IterativePSFPhotometry,
        PSFPhotometry,
        SourceGrouper,
    )
except (ModuleNotFoundError, ImportError):
    raise ValueError("'psf photometry requries photutils to be installed")
from trm import cline

import hipercam as hcam
from hipercam.cline import Cline
from hipercam.psf_reduction import MoffatPSF, create_psf_model

# DEBUG writes  out stars to FITS table for debugging purposes. It is not needed for normal operation.
DEBUG = False


__all__ = [
    "psfaper",
]


#############################################################
#
# psfaper -- defines apertures for PSF fitting given an image
#
#############################################################
def psfaper(args=None):
    """``psfaper mccd aper ccd [linput imwidth imheight imxpos imypos memory] nx method
    thresh niters gfac msub [cmap] iset (ilo ihi | plo phi) [method (beta betafix
    betamax) fwhm fwfix (fwmin) shbox smooth splot fhbox read gain
    rejthresh]``

    Definition of apertures for PSF photometry. This occurs in three
    steps - first the user draws a box around the region of interest,
    then the user selects a reference star which is used to measure the
    PSF. Finally, all the stars in the region of interest are located
    following an iterative application of FINDing stars, FITting the image
    using a PSF model, SUBTRACTing the model.

    An aperture file is automatically created, which can be edited
    with |setaper| if necessary.

    Parameters:

      mccd   : string
         name of an MCCD file, as produced by e.g. 'grab'

      aper   : string
         the name of an aperture file. If it exists it will be read so that
         apertures can be added to it. If it does not exist, it will be
         created on exiting the routine. The aperture files are is a fairly
         readable / editiable text format

      ccd    : string
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      linput  : string [hidden]
         sets the way in which the apertures are labelled. 'n' = numerical
         input, with the program just incrementing by 1 for each successive
         aperture; 's' = single character (without requiring the user to hit
         hit <CR>); 'm' = multi-character, ending with <CR>.
         Allowed characters are 0-9, a-z, A-Z, no spaces or punctuation, but a
         single '0' on its own is not permitted.

      imwidth : float [hidden]
         image display plot width in inches
         (0 for default, and note that this will also set the height to its default)

      imheight : float [hidden]
         image display plot height in inches
         (0 for default, and note that this will also set the height to its default)

      imxpos : int [hidden]
         image display plot X position in pixels

      imypos : int [hidden]
         image display plot Y position in pixels

      memory : bool [hidden]
         stores the plot display dimensions (width, height) and position
         (xpos, ypos) at the end as defaults for next time.
         This helps if you are repeatedly re-starting the script.

      nx     : int
         number of panels across to display, prompted if more than one CCD is
         to be plotted.

      method : string
         this defines the profile function. Either a gaussian or a moffat
         profile, 'g' or 'm'.  The latter should usually be best.

      thresh : float
         this sets the threshold for detection of stars in the image,
         in multiples of the background RMS

      niters : int
         when detecting stars, this sets the number of iterations of the
         FIND-FIT-SUBTRACT loop.

      gfac : float
         in PSF fitting, stars are split into groups, and each group is fit
         seperately. This is to avoid fitting models with large numbers of
         free parameters. This number, multiplied by the FWHM, gives the
         maximum seperation for stars to be fitted within the same group.

      msub   : bool
         True/False to subtract median from each window before scaling

      cmap : str [hidden]
         The colour map to use. "Greys" is the usual; "Greys_r" reverses it.
         There are many others; typing an incorrect one will give a list. "none"
         for matplotlib default.

      iset   : string [single character]
         determines how the intensities are determined. There are three
         options: 'a' for automatic simply scales from the minimum to the
         maximum value found on a per CCD basis. 'd' for direct just takes two
         numbers from the user. 'p' for percentile dtermines levels based upon
         percentiles determined from the entire CCD on a per CCD bais.

      ilo    : float [if iset=='d']
         lower intensity level

      ihi    : float [if iset=='d']
         upper intensity level

      plo    : float [if iset=='p']
         lower percentile level

      phi    : float [if iset=='p']
         upper percentile level

      method : string [hidden]
         this defines the profile fitting method, if profile fitting is used
         to refine the aperture position. Either a gaussian or a moffat
         profile, 'g' or 'm'.  The latter should usually be best.

      beta : float [if method == 'm'; hidden]
         default Moffat exponent

      betafix : bool [if method == 'm'; hidden]
         fix beta or not

      beta_max : bool [if method == 'm' and not betafix; hidden]
         maximum value of beta. Moffat profiles are degenerate
         with gaussians at large beta so the idea is to avoid wandering
         to huge beta.

      fwhm : float [hidden]
         default FWHM, unbinned pixels.

      fwfix: bool [hidden]
         don't fit the FWHM. Can be more robust; the position is still fitted.

      fwmin : float [if not fwfix; hidden]
         minimum FWHM to allow, unbinned pixels.

      shbox  : float [hidden]
         half width of box for searching for a star, unbinned pixels. The
         brightest target in a region +/- shbox around an intial position will
         be found. 'shbox' should be large enough to allow for likely changes
         in position from frame to frame, but try to keep it as small as you
         can to avoid jumping to different targets and to reduce the chances
         of interference by cosmic rays.

      smooth : float [hidden]
         FWHM for gaussian smoothing, binned pixels. The initial position for
         fitting is determined by finding the maximum flux in a smoothed
         version of the image in a box of width +/- shbox around the starter
         position. Typically should be comparable to the stellar width. Its
         main purpose is to combat cosmic rays which tend only to occupy a
         single pixel.

      fhbox  : float [hidden]
         half width of box for profile fit, unbinned pixels. The fit box is
         centred on the position located by the initial search. It should
         normally be > ~2x the expected FWHM.

      read   : float [hidden]
         readout noise, RMS ADU, for assigning uncertainties

      gain   : float [hidden]
         gain, ADU/count, for assigning uncertainties

      rejthresh  : float [hidden]
         thresh rejection threshold

    Various standard keyboard shortcuts (e.g. 's' to save) are disabled as
    they just confuse things and are of limited use in setaper in any case.

    """

    command, args = cline.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:
        # register parameters
        cl.register("mccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("aper", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("linput", Cline.LOCAL, Cline.HIDE)
        cl.register("imwidth", Cline.LOCAL, Cline.HIDE)
        cl.register("imheight", Cline.LOCAL, Cline.HIDE)
        cl.register("imxpos", Cline.LOCAL, Cline.HIDE)
        cl.register("imypos", Cline.LOCAL, Cline.HIDE)
        cl.register("memory", Cline.LOCAL, Cline.HIDE)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("method", Cline.LOCAL, Cline.PROMPT)
        cl.register("thresh", Cline.LOCAL, Cline.PROMPT)
        cl.register("niters", Cline.LOCAL, Cline.PROMPT)
        cl.register("gfac", Cline.LOCAL, Cline.PROMPT)
        cl.register("msub", Cline.GLOBAL, Cline.PROMPT)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
        cl.register("iset", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ilo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ihi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("plo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("phi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("beta", Cline.LOCAL, Cline.HIDE)
        cl.register("betafix", Cline.LOCAL, Cline.HIDE)
        cl.register("betamax", Cline.LOCAL, Cline.HIDE)
        cl.register("fwhm", Cline.LOCAL, Cline.HIDE)
        cl.register("fwfix", Cline.LOCAL, Cline.HIDE)
        cl.register("fwmin", Cline.LOCAL, Cline.HIDE)
        cl.register("shbox", Cline.LOCAL, Cline.HIDE)
        cl.register("smooth", Cline.LOCAL, Cline.HIDE)
        cl.register("fhbox", Cline.LOCAL, Cline.HIDE)
        cl.register("read", Cline.LOCAL, Cline.HIDE)
        cl.register("gain", Cline.LOCAL, Cline.HIDE)
        cl.register("rejthresh", Cline.LOCAL, Cline.HIDE)

        # get inputs
        mccd = cl.get_value("mccd", "frame to plot", cline.Fname("hcam", hcam.HCAM))
        mccd = hcam.MCCD.read(mccd)

        aper = cl.get_value(
            "aper", "name of aperture file", cline.Fname("hcam", hcam.APER, exist=False)
        )

        if os.path.exists(aper):
            # read in old apertures
            mccdaper = hcam.MccdAper.read(aper)
            print("Loaded existing file = {:s}".format(aper))
        else:
            # create empty container
            mccdaper = hcam.MccdAper()
            print("No file called {:s} exists; will create from scratch".format(aper))

        # define the panel grid
        nxdef = cl.get_default("nx", 3)

        max_ccd = len(mccd)
        if max_ccd > 1:
            ccd = cl.get_value("ccd", "CCD(s) to plot [0 for all]", "0")
            if ccd == "0":
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(mccd.keys())

        imwidth = cl.get_value(
            "imwidth", "image plot width [inches, 0 for default]", 0.0, 0.0
        )

        imheight = cl.get_value(
            "imheight", "image plot height [inches, 0 for default]", 0.0, 0.0
        )

        imxpos = cl.get_value(
            "imxpos", "image plot X position [screen pixels]", 0,
        )

        imypos = cl.get_value(
            "imypos", "image plot Y position [screen pixels]", 0,
        )

        memory = cl.get_value(
            "memory", "remember final plot geometry for next time?", True
        )

        # number of panels in X
        if len(ccds) > 1:
            nxdef = min(len(ccds), nxdef)
            cl.set_default("nx", nxdef)
            nx = cl.get_value("nx", "number of panels in X", 3, 1)
        else:
            nx = 1

        # PSF fitting stuff
        method = cl.get_value(
            "method", "fit method g(aussian) or m(offat)", "m", lvals=["g", "m"]
        )
        if method == "m":
            beta = cl.get_value("beta", "initial exponent for Moffat fits", 5.0, 0.5)
            beta_fix = cl.get_value("betafix", "fix beta at start value?", False)
            if beta_fix:
                beta_max = beta
            else:
                beta_max = cl.get_value(
                    "betamax", "maximum beta to allow", max(beta, 20)
                )
            beta = cl.get_value("beta", "initial exponent for Moffat fits", 5.0, 0.5)
        else:
            beta = 0.0
            beta_fix = True
            beta_max = 0.0

        fwhm = cl.get_value(
            "fwhm", "initial FWHM [unbinned pixels] for profile fits", 6.0, 1.0
        )
        fwhm_fix = cl.get_value("fwfix", "fix the FWHM at start value?", False)
        if fwhm_fix:
            fwhm_min = fwhm
        else:
            fwhm_min = cl.get_value(
                "fwmin",
                "minimum FWHM to allow [unbinned pixels]",
                min(1.5, fwhm),
                0.01,
                fwhm,
            )

        gfac = cl.get_value(
            "gfac", "multiple of FWHM used to group stars for fitting", 2.0, 0.1
        )

        thresh = cl.get_value(
            "thresh", "threshold for object detection (multiple of sky RMS)", 3.0, 0.1
        )

        niters = cl.get_value(
            "niters", "number of iterations of FIND-FIT-SUBTRACT", 2, 1
        )

        # define the display intensities
        msub = cl.get_value("msub", "subtract median from each window?", True)

        cmap = cl.get_value("cmap", "colour map to use ['none' for mpl default]", "Greys")
        cmap = None if cmap == "none" else cmap

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

        nxmax, nymax = 0, 0
        for cnam in ccds:
            nxmax = max(nxmax, mccd[cnam].nxtot)
            nymax = max(nymax, mccd[cnam].nytot)

        # might be worth trying to improve this at some point
        xlo, xhi, ylo, yhi = 0, nxmax + 1, 0, nymax + 1

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
            "fhbox", "half width of box for profile fit [unbinned pixels]", 21.0, 3.0
        )
        read = cl.get_value("read", "readout noise, RMS ADU", 3.5)
        gain = cl.get_value("gain", "gain, ADU/e-", 1.1)
        rejthresh = cl.get_value(
            "rejthresh", "RMS rejection threshold for sky fitting", 4.0
        )

    # create the plot
    fig, axs = hcam.mpl.plot_frame(
        mccd=mccd,
        cnams=ccds,
        ncols=nx,
        msub=msub,
        cmap=cmap,
        iset=iset,
        ilim=(plo, phi) if iset == "p" else (ilo, ihi),
        region=(xlo, xhi, ylo, yhi),
        title='psfaper',
        figsize=(imwidth, imheight),
        position=(imxpos, imypos)
    )

    # add any existing apertures to the plot, storing a reference to the plot objects
    pobjs = {}
    for cnam, ax in axs.items():
        if cnam in mccdaper:
            pobjs[cnam] = hcam.mpl.pCcdAper(ax, mccdaper[cnam])
        else:
            # create in an empty group to be added to later
            mccdaper[cnam] = hcam.CcdAper()
            pobjs[cnam] = {}

    # print out instructions for the user
    print(
        """
Now use the mouse and the pan/zoom tools to zoom in to the region of interest (ROI) for
PSF photometry. All existing apertures inside this region will be retained, but apertures
outside this region will be deleted.

Once you are happy with your selection, use the key commands listed below to select one or
more bright, isolated, stars to use as references to determine the PSF. These stars will also be
set as reference objects in the final aperture file.

Key commands:

  a(dd)      : add an aperture
  d(elete)   : delete an aperture
  u(nlink)   : unlink axes for all CCDs so that you can tweak ROI for CCDs independently
  q(uit)     : quit interactive step and find all stars in ROI.

Hitting 'd' will delete the aperture nearest to the cursor, as long as it is
close enough.
"""
    )

    # create the picker class and connect it to the figure
    picker = PickRef(
        mccd,
        mccdaper,
        fig,
        axs,
        pobjs,
        aper,
        method,
        beta,
        beta_max,
        beta_fix,
        fwhm,
        fwhm_min,
        fwhm_fix,
        gfac,
        niters,
        thresh,
        shbox,
        smooth,
        fhbox,
        read,
        gain,
        rejthresh,
    )
    fig.canvas.mpl_connect("key_press_event", picker._keyPressEvent)

    # save window geometry on close
    if memory:
        def save_geometry(event):
            imwidth, imheight, imxpos, imypos = hcam.mpl.store_window(
                fig, command, args,
            )
            print(
                "Image display memory updated [memory=True]:",
                f"imwidth={imwidth:.1f}, imheight={imheight:.1f}, "
                f"imxpos={imxpos}, imypos={imypos}"
            )
        fig.canvas.mpl_connect("close_event", save_geometry)

    # finally set up the plot for interactive mode and display
    hcam.mpl.setup_plot(fig)
    plt.show()


class PickRef:
    """Class to zoom into ROI and select reference stars."""

    def __init__(
        self,
        mccd,
        mccdaper,
        fig,
        axs,
        pobjs,
        apernam,
        method,
        beta,
        beta_max,
        beta_fix,
        fwhm,
        fwhm_min,
        fwhm_fix,
        gfac,
        niters,
        thresh,
        shbox,
        smooth,
        fhbox,
        read,
        gain,
        rejthresh,
    ):
        # image content data
        self.mccd = mccd
        self.mccdaper = mccdaper

        # figure and axes
        self.fig = fig
        self.axs = axs
        self.axnams = {ax: cnam for cnam, ax in axs.items()}  # reverse for events
        self.pobjs = pobjs

        # output file
        self.apernam = apernam

        # fitting params
        self.method = method
        self.beta = beta
        self.beta_max = beta_max
        self.beta_fix = beta_fix
        self.fwhm = fwhm
        self.fwhm_min = fwhm_min
        self.fwhm_fix = fwhm_fix
        self.gfac = gfac
        self.niters = niters
        self.thresh = thresh
        self.fwhm_min = fwhm_min
        self.shbox = shbox
        self.smooth = smooth
        self.fhbox = fhbox
        self.read = read
        self.gain = gain
        self.rejthresh = rejthresh

        # dictionary to store fit data for PSFs of reference stars
        self.psf_data = dict()

    def _keyPressEvent(self, event):
        # shortcuts
        key, x, y, ax = event.key, event.xdata, event.ydata, event.inaxes
        if ax is not None:
            # store information in attributes accessible to all methods, for
            # later access: name, the key hit, the axes instance, x, y
            self._cnam = self.axnams[ax]
            self._key = key
            self._ax = ax
            self._x = x
            self._y = y

            if key == "a":
                self._add()
            elif key == "d":
                self._delete()
            elif key == "u":
                self._unlinkaxes()
            elif key == "q":
                self._find_stars_and_quit()

    def _unlinkaxes(self):
        for ax in self.axs.values():
            shared_x = ax.get_shared_x_axes()
            shared_y = ax.get_shared_y_axes()
            for child in shared_x.get_siblings(ax):
                shared_x.remove(child)
            for child in shared_y.get_siblings(ax):
                shared_y.remove(child)

    def _add(self):
        """
        Add a star to be used as a reference aperture for PSF fitting
        """
        self._label = self.mccdaper[self._cnam].next_label()

        # fitting
        ccd = self.mccd[self._cnam]

        # check that the selected position is inside a window
        wnam = ccd.inside(self._x, self._y, 2)
        if wnam is None:
            print(
                "  *** selected position ({:.1f},{:.1f}) not in a window;"
                " should not occur".format(self._x, self._y),
                file=sys.stderr,
            )
            return

        # get Window around the selected position
        wind = ccd[wnam].window(
            self._x - self.shbox,
            self._x + self.shbox,
            self._y - self.shbox,
            self._y + self.shbox,
        )

        try:
            # carry out initial search
            x, y, peak = wind.search(self.smooth, 0, 0, 0, False, True, 0)

            # now for a more refined fit. First extract fit Window
            fwind = ccd[wnam].window(
                x - self.fhbox, x + self.fhbox, y - self.fhbox, y + self.fhbox
            )
            sky = np.percentile(fwind.data, 25)

            # uncertainties
            sigma = np.sqrt(self.read**2 + np.maximum(0, fwind.data) / self.gain)

            # refine the Aperture position by fitting the profile
            (
                (sky, height, x, y, fwhm, beta),
                (esky, eheight, ex, ey, efwhm, ebeta),
                (wfit, X, Y, chisq, nok, nrej, npar, nfev, message),
            ) = hcam.fitting.combFit(
                fwind,
                sigma,
                self.method,
                sky,
                peak - sky,
                x,
                y,
                self.fwhm,
                self.fwhm_min,
                self.fwhm_fix,
                self.beta,
                self.beta_max,
                self.beta_fix,
                self.rejthresh,
            )

            print("Aperture {:s}: {:s}".format(self._label, message))
            self._x = x
            self._y = y

        except hcam.HipercamError as err:
            print(err, file=sys.stderr)
            # fit failed.
            return

        # create and add reference aperture
        aper = hcam.Aperture(self._x, self._y, 5, 10, 12, True)
        self.mccdaper[self._cnam][self._label] = aper

        # add fit params to object
        wf = 1.0 / efwhm**2 if efwhm > 0 else 0
        wb = 1.0 / ebeta**2 if ebeta > 0 else 0
        if self._cnam not in self.psf_data:
            self.psf_data[self._cnam] = dict(
                fsum=wf * fwhm, wfsum=wf, bsum=wb * beta, wbsum=wb
            )
        else:
            self.psf_data[self._cnam]["fsum"] += wf * fwhm
            self.psf_data[self._cnam]["wfsum"] += wf
            self.psf_data[self._cnam]["bsum"] += wb * beta
            self.psf_data[self._cnam]["wbsum"] += wb

        # add aperture to the plot, store plot objects
        self.pobjs[self._cnam][self._label] = hcam.mpl.pAper(
            self._ax, aper, self._label
        )

        # make sure it appears
        plt.draw()

        print(
            "added aperture {:s} to CCD {:s} at x,y = {:.2f},{:.2f}".format(
                self._label, self._cnam, self._x, self._y
            )
        )

    def _delete(self):
        """This deletes the nearest aperture to the currently selected
        position, if it is near enough. It will also break any links
        from other Apertures to the Aperture that is deleted.

        'Near enough' is defined as within max(rtarg, min(100, max(20,2*rsky2)))
        of the aperture centre.

        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is not None and dmin < 60:
            # near enough for deletion
            for obj in self.pobjs[self._cnam][apnam]:
                obj.remove()

            # delete Aperture from containers
            del self.pobjs[self._cnam][apnam]
            del self.mccdaper[self._cnam][apnam]

            # update plot
            plt.draw()
            print('  deleted aperture "{:s}"'.format(apnam))

        else:
            print("  found no aperture near enough the cursor position for deletion")

    def _find_aper(self):
        """Finds the nearest aperture to the currently selected position,

        It returns (aper, apnam, dmin) where aper is the Aperture, apnam its
        label, and dmin is the minimum distance. These are all returned as
        None if no suitable Aperture is found.

        """
        dmin = None
        apmin = None
        anmin = None
        for apnam, aper in self.mccdaper[self._cnam].items():
            dist = np.sqrt((aper.x - self._x) ** 2 + (aper.y - self._y) ** 2)
            if dmin is None or dist < dmin:
                dmin = dist
                apmin = aper
                anmin = apnam

        return (apmin, anmin, dmin)

    def _find_stars_and_quit(self):
        """
        Runs PSF photometry on the region of interest, finds the stars, adds them as apertures.

        For each CCD the PSF shape parameters are estimated from the reference stars,
        at which point the PSF shape parameters are held fixed and multiple iterations
        of the FIND-FIT-SUBTRACT loop are run to find more stars.

        The final apertures are then saved to the aperture file and the program exits.
        """
        for cnam, ccdaper in self.mccdaper.items():
            xlo, xhi = self.axs[cnam].get_xlim()
            ylo, yhi = self.axs[cnam].get_ylim()

            if cnam not in self.psf_data:
                warnings.warn(
                    "no (new) reference stars for CCD{} - skipping".format(cnam)
                )
                continue

            fwhm = -1
            beta = -1
            if self.psf_data[cnam]["wfsum"] > 0:
                fwhm = self.psf_data[cnam]["fsum"] / self.psf_data[cnam]["wfsum"]
            if self.psf_data[cnam]["wbsum"] > 0:
                beta = self.psf_data[cnam]["bsum"] / self.psf_data[cnam]["wbsum"]

            if fwhm < 0:
                warnings.warn("no fwhm measured for CCD{} - skipping".format(cnam))
                continue

            if self.method == "m" and beta < 0:
                warnings.warn("no beta measured for CCD{} - skipping".format(cnam))
                continue

            print("fitting CCD{}:".format(cnam))
            xlocs, ylocs = daophot(
                cnam,
                self.mccd[cnam],
                ccdaper,
                xlo,
                xhi,
                ylo,
                yhi,
                self.niters,
                self.method,
                fwhm,
                beta,
                self.gfac,
                self.thresh,
                self.rejthresh,
                self.read,
                self.gain,
            )

            # let's find which of our positions best matches our reference apertures
            # save ref apertures
            ref_apertures = [aper for aper in ccdaper.values()]

            # add apertures for each position
            self.mccdaper[cnam] = hcam.CcdAper()
            buffer = 0
            for x, y in zip(xlocs, ylocs):
                if x > xlo and x < xhi and y > ylo and y < yhi:
                    aper = hcam.Aperture(x, y, 5, 10, 15, False)
                    self.mccdaper[cnam][str(buffer + 1)] = aper
                    buffer += 1

            # find closest aperture to each ref aperture and make it
            # a reference aperture for PSF reduction
            for aper in ref_apertures:
                self._x = aper.x
                self._y = aper.y
                self._cnam = cnam
                aper, apnam, dmin = self._find_aper()
                if dmin is not None and dmin < 5:
                    # close enough, let's make this an aperture
                    self.mccdaper[cnam][apnam].ref = True

        # done. write aperture file
        plt.close()
        self.mccdaper.write(self.apernam)
        print("\nApertures saved to {:s}.\nBye".format(self.apernam))


def daophot(
    cnam,
    ccd,
    ccdaper,
    xlo,
    xhi,
    ylo,
    yhi,
    niters,
    method,
    fwhm,
    beta,
    gfac,
    thresh,
    rejthresh,
    read,
    gain,
):
    """
    Perform PSF photometry on region of CCD.

    The PSF shape parameters are determined from the reference stars, at which point
    the PSF shape parameters are held fixed and multiple iterations of the FIND-FIT-SUBTRACT
    loop are run to find more stars.

    Parameters
    ----------
    cnam : str
        name of CCD
    ccd : hcam.Ccd
        CCD to perform photometry on
    ccdaper : hcam.CcdAper
        reference apertures to use for determining PSF shape parameters
    xlo, xhi, ylo, yhi : float
        limits of region to perform photometry on, in device coordinates
    niters : int
        number of iterations of FIND-FIT-SUBTRACT to perform
    method : str
        PSF model to use, either 'm' for Moffat, 'g' for Gaussian
    fwhm : float
        initial FWHM to use for PSF model, in unbinned pixels
    beta : float
        initial beta to use for Moffat PSF model, ignored if method is 'g'
    gfac : float
        multiple of FWHM used to group stars for fitting.
        Stars within gfac*FWHM of each other are fitted simultaneously in the FIND-FIT-SUBTRACT loop.
    thresh : float
        threshold for object detection in FIND-FIT-SUBTRACT, in multiples of the background RMS
    rejthresh : float
        RMS rejection threshold for sky fitting, in multiples of the background RMS.
    read : float
        readout noise, RMS ADU, for assigning uncertainties
    gain : float
        gain, ADU/count, for assigning uncertainties
    """
    # first check we are within a single window
    wnam1 = ccd.inside(xlo, ylo, 2)
    wnam2 = ccd.inside(xhi, yhi, 2)
    if wnam1 != wnam2:
        raise hcam.HipercamError(
            "PSF photometry cannot currently be run across seperate windows"
        )
    wnam = wnam1

    # background stats from whole windpw
    # estimate background RMS
    wind = ccd[wnam]

    rms_func = MADStdBackgroundRMS(sigma_clip=SigmaClip(sigma=rejthresh))
    bkg_rms = rms_func(wind.data)
    bkg_func = MMMBackground(sigma_clip=SigmaClip(sigma=rejthresh))
    bkg = bkg_func(wind.data)
    print("  Background estimate = {}, BKG RMS = {}".format(bkg, bkg_rms))

    # crop window to ROI
    wind = ccd[wnam].window(xlo, xhi, ylo, yhi)
    # uncertainties
    sigma = np.sqrt(read**2 + np.maximum(0, wind.data) / gain)

    # correct FWHM for binning
    fwhm /= wind.xbin
    if method == "m":
        psf_model = MoffatPSF(beta=beta, x_fwhm=fwhm, y_fwhm=fwhm)
        psf_model_name = "moffat"
        # no need to set fixed and bound, they have sensible defaults.
    else:
        psf_model = GaussianPRF(x_fwhm=fwhm, y_fwhm=fwhm)
        psf_model.x_fwhm.fixed = False
        psf_model.y_fwhm.fixed = False
        psf_model.theta.fixed = False
        psf_model.theta.bounds = (-90, 90)
        psf_model_name = "gaussian"

    # define region to extract around positions for fits
    fitshape = int(5 * fwhm)
    # ensure odd
    if fitshape % 2 == 0:
        fitshape += 1

    # Step 1: fit the reference stars to determine the PSF parameters.
    # get pixel positions of reference apertures inside window
    xpos, ypos = list(
        zip(
            *[
                (wind.x_pixel(aper.x), wind.y_pixel(aper.y))
                for aper in ccdaper.values()
                if aper.ref and wind.distance(aper.x, aper.y) > 0.0
            ]
        )
    )
    if len(xpos) == 0:
        raise hcam.HipercamError(
            f"No reference apertures found in CCD {cnam} for PSF fitting"
        )
    reference_positions = Table(names=["x_0", "y_0"], data=(xpos, ypos))

    photometry_task = PSFPhotometry(
        psf_model=psf_model,
        aperture_radius=1.7 * fwhm,
        fit_shape=fitshape,
    )
    # do the PSF photometry
    photom_results = photometry_task(
        wind.data - bkg, error=sigma, init_params=reference_positions
    )

    # now create the PSF model from fitting the reference stars
    psf_model = create_psf_model(photom_results, psf_model_name, fixed_positions=False)
    output_params = ["x_fwhm", "y_fwhm", "theta"]
    if method == "m":
        output_params.append("beta")
    print("Fitted PSF parameters:", end=" ")
    for param in output_params:
        value = getattr(psf_model, param).value
        print(f"{param}: {value:.3f}", end=", ")
    print()

    # Step 2: run the FIND-FIT-SUBTRACT loop to find all stars in the region.
    photometry_task = IterativePSFPhotometry(
        psf_model=psf_model,
        fit_shape=fitshape,
        finder=DAOStarFinder(threshold=thresh * bkg_rms, fwhm=fwhm),
        grouper=SourceGrouper(gfac * fwhm),
        maxiters=niters,
        mode="all",
        aperture_radius=1.7 * fwhm,
        local_bkg_estimator=None,
    )

    results = photometry_task(wind.data - bkg, error=sigma)

    # filter out junk fits
    results = results[results["flags"] == 0]

    tiny = 1e-30
    bad_errs = (
        (results["flux_err"] < tiny)
        | (results["x_err"] < tiny)
        | (results["y_err"] < tiny)
    )
    results = results[~bad_errs]

    if DEBUG:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results.write("table_{}.fits".format(cnam), overwrite=True)
    print(f"  found {len(results)} stars")

    xlocs, ylocs = results["x_fit"], results["y_fit"]

    # convert to device coordinates
    xlocs = wind.x(xlocs)
    ylocs = wind.y(ylocs)
    return xlocs, ylocs


if __name__ == "__main__":
    psfaper()
