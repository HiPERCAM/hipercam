import sys
import os
import warnings
from copy import deepcopy

import numpy as np
import matplotlib as mpl

# re-configure the cursors: backend specific.
# aim to get rid of irritating 'hand' icon in
# favour of something pointier.

backend = mpl.get_backend()

if backend == "Qt4Agg" or "Qt5Agg":
    from matplotlib.backends.backend_qt5 import cursord as curs
elif backend == "GTK3agg":
    from matplotlib.backends.backend_gtk3 import cursord as curs
else:
    curs = None

if curs is not None:
    from matplotlib.backend_bases import cursors

    try:
        curs[cursors.HAND] = curs[cursors.POINTER]
        curs[cursors.WAIT] = curs[cursors.POINTER]
        curs[cursors.SELECT_REGION] = curs[cursors.POINTER]
        curs[cursors.MOVE] = curs[cursors.POINTER]
    except AttributeError:
        pass

import matplotlib.pyplot as plt

from photutils.psf import DAOPhotPSFPhotometry, IntegratedGaussianPRF
from photutils.background import MADStdBackgroundRMS, MMMBackground
from astropy.stats import gaussian_fwhm_to_sigma, SigmaClip

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline
from hipercam.scripts.psf_reduce import MoffatPSF

__all__ = [
    "psfaper",
]


#############################################################
#
# psfaper -- defines apertures for PSF fitting given an image
#
#############################################################
def psfaper(args=None):
    """``psfaper mccd aper ccd [linput width height] nx method thresh
    niters gfac msub iset (ilo ihi | plo phi) [method (beta betafix
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

      width  : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

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

    Some aspects of the usage of matplotlib in psfaper are tricky. It is
    possible that particular 'backends' will cause problems. I have tested
    this with Qt4Agg, Qt5agg and GTK3Agg. One aspect is the cursor icon in pan
    mode is a rather indistinct hand where one can't tell what is being
    pointed at. I have therefore suppressed this, but only for the tested
    backends. Others would need require further investigation.

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("mccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("aper", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("linput", Cline.LOCAL, Cline.HIDE)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("method", Cline.LOCAL, Cline.PROMPT)
        cl.register("thresh", Cline.LOCAL, Cline.PROMPT)
        cl.register("niters", Cline.LOCAL, Cline.PROMPT)
        cl.register("gfac", Cline.LOCAL, Cline.PROMPT)
        cl.register("msub", Cline.GLOBAL, Cline.PROMPT)
        cl.register("iset", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ilo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ihi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("plo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("phi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("method", Cline.LOCAL, Cline.HIDE)
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
            print(
                "No file called {:s} exists; " "will create from scratch".format(aper)
            )

        # define the panel grid
        try:
            nxdef = cl.get_default("nx")
        except KeyError:
            nxdef = 3

        max_ccd = len(mccd)
        if max_ccd > 1:
            ccd = cl.get_value("ccd", "CCD(s) to plot [0 for all]", "0")
            if ccd == "0":
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(mccd.keys())

        width = cl.get_value("width", "plot width (inches)", 0.0)
        height = cl.get_value("height", "plot height (inches)", 0.0)

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
                    "betamax", "maximum beta to allow", max(beta,20) 
                )
            beta = cl.get_value("beta", "initial exponent for Moffat fits", 5.0, 0.5)
        else:
            beta = 0.
            beta_fix = True
            beta_max = 0.

        fwhm = cl.get_value(
            "fwhm", "initial FWHM [unbinned pixels] for profile fits", 6.0, 1.0
        )
        fwhm_fix = cl.get_value("fwfix", "fix the FWHM at start value?", False)
        if fwhm_fix:
            fwhm_min = fwhm
        else:
            fwhm_min = cl.get_value(
                "fwmin", "minimum FWHM to allow [unbinned pixels]", min(1.5,fwhm), 0.01, fwhm
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

        iset = cl.get_value(
            "iset",
            "set intensity a(utomatically)," " d(irectly) or with p(ercentiles)?",
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
            "half width of box for initial" " location of target [unbinned pixels]",
            11.0,
            2.0,
        )
        smooth = cl.get_value(
            "smooth",
            "FWHM for smoothing for initial object" " detection [binned pixels]",
            6.0,
        )

        fhbox = cl.get_value(
            "fhbox", "half width of box for profile fit" " [unbinned pixels]", 21.0, 3.0
        )
        read = cl.get_value("read", "readout noise, RMS ADU", 3.0)
        gain = cl.get_value("gain", "gain, ADU/e-", 1.0)
        rejthresh = cl.get_value(
            "rejthresh", "RMS rejection threshold for sky fitting", 4.0
        )

    # Inputs obtained.

    # re-configure keyboard shortcuts to avoid otherwise confusing behaviour
    # quit_all does not seem to be universal, hence the try/except
    try:
        mpl.rcParams["keymap.back"] = ""
        mpl.rcParams["keymap.forward"] = ""
        mpl.rcParams["keymap.fullscreen"] = ""
        mpl.rcParams["keymap.grid"] = ""
        mpl.rcParams["keymap.home"] = ""
        mpl.rcParams["keymap.pan"] = ""
        mpl.rcParams["keymap.quit"] = ""
        mpl.rcParams["keymap.save"] = ""
        mpl.rcParams["keymap.pan"] = ""
        mpl.rcParams["keymap.save"] = ""
        mpl.rcParams["keymap.xscale"] = ""
        mpl.rcParams["keymap.yscale"] = ""
        mpl.rcParams["keymap.zoom"] = ""
    except KeyError:
        pass

    # start plot
    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width, height))
    else:
        fig = plt.figure()

    # get the navigation toolbar.
    toolbar = fig.canvas.manager.toolbar

    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    # we need to store some stuff
    ax = None
    cnams = {}
    anams = {}

    # this is a container for all the objects used to plot apertures to allow
    # deletion. This is Group of Group objects supporting tuple storage. The
    # idea is that pobjs[cnam][anam] returns the objects used to plot aperture
    # anam of CCD cnam. It is initially empty,
    pobjs = hcam.Group(hcam.Group)

    for n, cnam in enumerate(ccds):
        if ax is None:
            axes = ax = fig.add_subplot(ny, nx, n + 1)
            axes.set_aspect("equal", adjustable="box")
            axes.set_xlim(xlo, xhi)
            axes.set_ylim(ylo, yhi)
        else:
            axes = fig.add_subplot(ny, nx, n + 1, sharex=ax, sharey=ax)
            axes.set_aspect("equal")

        if msub:
            # subtract median from each window
            pccd = deepcopy(mccd)
            for wind in pccd[cnam].values():
                wind -= wind.median()
        else:
            pccd = mccd

        hcam.mpl.pCcd(
            axes, pccd[cnam], iset, plo, phi, ilo, ihi, "CCD {:s}".format(cnam)
        )

        # keep track of the CCDs associated with each axes
        cnams[axes] = cnam

        # and axes associated with each CCD
        anams[cnam] = axes

        if cnam in mccdaper:
            # plot any pre-existing apertures, keeping track of
            # the plot objects
            pobjs[cnam] = hcam.mpl.pCcdAper(axes, mccdaper[cnam])

        else:
            # add in an empty CcdApers for any CCD not already present
            mccdaper[cnam] = hcam.CcdAper()

            # and an empty container for any new plot objects
            pobjs[cnam] = hcam.Group(tuple)

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

    try:
        plt.tight_layout()
    except:
        pass

    # create a class for picking reference star and zooming in
    picker = PickRef(
        mccd,
        cnams,
        anams,
        toolbar,
        fig,
        mccdaper,
        method,
        beta,
        beta_max,
        fwhm,
        fwhm_min,
        fwhm_fix,
        gfac,
        niters,
        thresh,
        fwhm_min,
        shbox,
        smooth,
        fhbox,
        read,
        gain,
        rejthresh,
        pobjs,
        aper,
    )

    # squeeze space a bit
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    # finally show stuff ....
    plt.show()


# the next class is where all the action occurs. A rather complicated matter
# of handling events. note that standard terminal input with 'input' becomes
# impossible, explaining some of the weirdness. Effectively the class is used
# here to define a scope for variables that would otherwise be treated as globals
class PickRef:
    """
    Class to zoom into ROI and select reference stars.
    """

    def __init__(
        self,
        mccd,
        cnams,
        anams,
        toolbar,
        fig,
        mccdaper,
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
        fwhm_min,
        shbox,
        smooth,
        fhbox,
        read,
        gain,
        rejthresh,
        pobjs,
        apernam,
    ):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect("key_press_event", self._keyPressEvent)
        self.mccd = mccd
        self.cnams = cnams
        self.anams = anams
        self.toolbar = toolbar
        self.mccdaper = mccdaper
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
        self.pobjs = pobjs
        self.apernam = apernam
        # dictionary to store fit data for PSFs of reference stars
        self.psf_data = dict()

    def _keyPressEvent(self, event):
        # shortcuts
        key, x, y, axes = event.key, event.xdata, event.ydata, event.inaxes
        if axes is not None:
            # store information in attributes accessible to all methods, for
            # later access: name, the key hit, the axes instance, x, y
            self._cnam = self.cnams[axes]
            self._key = key
            self._axes = axes
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
        for axes in self.anams.values():
            shared_x = axes.get_shared_x_axes()
            shared_y = axes.get_shared_y_axes()
            for child in shared_x.get_siblings(axes):
                shared_x.remove(child)
            for child in shared_y.get_siblings(axes):
                shared_y.remove(child)

    def _add(self):
        """
        Add a star to be used as a reference aperture for PSF fitting
        """
        # numerical sequence input. Try to calculate
        # the largest number, label the new aperture
        # with one more
        high = 0
        for aper in self.mccdaper[self._cnam]:
            try:
                high = max(high, int(aper))
            except ValueError:
                pass
        self._buffer = str(high + 1)

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

            # refine the Aperture position by fitting the profile
            (
                (sky, height, x, y, fwhm, beta),
                (esky, eheight, ex, ey, efwhm, ebeta),
                (wfit, X, Y, sigma, chisq, nok, nrej, npar, nfev, message),
            ) = hcam.fitting.combFit(
                fwind,
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
                self.read,
                self.gain,
                self.rejthresh,
            )

            print("Aperture {:s}: {:s}".format(self._buffer, message))
            self._x = x
            self._y = y

        except hcam.HipercamError as err:
            print(err, file=sys.stderr)
            # fit failed.
            return

        # create and add reference aperture
        aper = hcam.Aperture(self._x, self._y, 5, 10, 12, True)
        self.mccdaper[self._cnam][self._buffer] = aper

        # add fit params to object
        wf = 1.0 / efwhm ** 2 if efwhm > 0 else 0
        wb = 1.0 / ebeta ** 2 if ebeta > 0 else 0
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
        self.pobjs[self._cnam][self._buffer] = hcam.mpl.pAper(
            self._axes, aper, self._buffer
        )

        # make sure it appears
        plt.draw()

        print(
            "added aperture {:s} to CCD {:s} at x,y = {:.2f},{:.2f}".format(
                self._buffer, self._cnam, self._x, self._y
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
            print("  found no aperture near enough " "the cursor position for deletion")

    def _find_aper(self):
        """Finds the nearest aperture to the currently selected position,

        It returns (aper, apnam, dmin) where aper is the Aperture, apnam its
        label, and dmin is the minimum distance. These are all returned as
        None if no suitable Aperture is found.

        """
        dmin = None
        apmin = None
        anmin = None
        for anam, aper in self.mccdaper[self._cnam].items():
            dist = np.sqrt((aper.x - self._x) ** 2 + (aper.y - self._y) ** 2)
            if dmin is None or dist < dmin:
                dmin = dist
                apmin = aper
                anmin = anam

        return (apmin, anmin, dmin)

    def _find_stars_and_quit(self):
        """
        Runs the PSF photometry on each ROI, and writes the aperture file
        """
        for cnam, ccdaper in self.mccdaper.items():
            xlo, xhi = self.anams[cnam].get_xlim()
            ylo, yhi = self.anams[cnam].get_ylim()

            if cnam not in self.psf_data:
                warnings.warn("no reference stars for CCD{} - skipping".format(cnam))
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
    cnam, ccd, xlo, xhi, ylo, yhi, niters, method, fwhm, beta, gfac, thresh, rejthresh
):
    """
    Perform iterative PSF photometry and star finding on region of CCD
    """
    print(xlo, ylo, xhi, yhi)
    # first check we are within a single window
    wnam1 = ccd.inside(xlo, ylo, 2)
    wnam2 = ccd.inside(xhi, yhi, 2)
    if wnam1 != wnam2:
        raise hcam.HipercamError(
            "PSF photometry cannot currently be run across seperate windows"
        )
    wnam = wnam1
    print(wnam)
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

    # correct FWHM for binning
    fwhm /= wind.xbin
    if method == "m":
        psf_model = MoffatPSF(fwhm, beta)
        print("  FWHM = {:.1f}, BETA={:.1f}".format(fwhm, beta))
    else:
        psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)
        print("  FWHM = {:.1f}".format(fwhm))

    # region to extract around positions for fits
    fitshape = int(5 * fwhm)
    # ensure odd
    if fitshape % 2 == 0:
        fitshape += 1

    photometry_task = DAOPhotPSFPhotometry(
        gfac * fwhm,
        thresh * bkg_rms,
        fwhm,
        psf_model,
        fitshape,
        niters=niters,
        sigma=rejthresh,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = photometry_task(wind.data - bkg)

    # filter out junk fits
    tiny = 1e-30
    bad_errs = (
        (results["flux_unc"] < tiny)
        | (results["x_0_unc"] < tiny)
        | (results["y_0_unc"] < tiny)
    )
    results = results[~bad_errs]

    results.write("table_{}.fits".format(cnam))
    print("  found {} stars".format(len(results)))

    xlocs, ylocs = results["x_fit"], results["y_fit"]

    # convert to device coordinates
    xlocs = wind.x(xlocs)
    ylocs = wind.y(ylocs)
    return xlocs, ylocs


if __name__ == "__main__":
    psfaper()
