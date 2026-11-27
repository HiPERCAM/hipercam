#!/usr/bin/env python
from __future__ import division, print_function

import os
import time
from contextlib import suppress

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interpolate
import scipy.spatial as sp
from astropy.stats import sigma_clip
from astropy.time import Time
from matplotlib.ticker import NullFormatter
from sklearn.cluster import KMeans
from trm import cline

import hipercam as hcam
from hipercam import spooler
from hipercam.ccd import crop_calib_frame_to
from hipercam.cline import Cline
from hipercam.extraction import findStars
from hipercam.instruments import available_instruments

__all__ = [
    "aligntool",
]


def findBestRigidTransform(x, y, xref, yref):
    """
    Given a set of points and *matching* reference points, find the optimal Rigid transform.

    A Rigid transform is a rotation and translation, i.e:

    B = R.A + t

    See https://scicomp.stackexchange.com/questions/6878/fitting-one-set-of-points-to-another-by-a-rigid-motion
    or http://nghiaho.com/?page_id=671

    for method.
    """

    A = np.column_stack((xref, yref))
    B = np.column_stack((x, y))

    # find centre of points
    centroidA = np.mean(A, axis=0)
    centroidB = np.mean(B, axis=0)

    # shift centres to origin
    Ashifted = A - centroidA
    Bshifted = B - centroidB

    # find dot product
    H = np.dot(Bshifted.T, Ashifted)

    # use SVD to find rotation
    u, s, v = np.linalg.svd(H)
    R = np.dot(v, u.T)
    if np.linalg.det(R) < 0:
        R = np.dot(v, np.array([1, -1]) * u.T)

    T = np.median(B - np.dot(A, R), axis=0)
    theta1 = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
    theta2 = np.degrees(-np.arctan2(R[0, 1], R[1, 1]))
    theta = 0.5 * (theta1 + theta2)
    return theta, T


def matchSpots(x, y, xref, yref):
    """
    Match pairs of spots, looking for the nearest neighbour.

    Parameters
    ----------
    x, y: `numpy.ndarray`
        spots to match
    xref, yref: `numpy.ndarray`
        reference positions

    Returns
    -------
    theta, offset: `numpy.ndarray`
        rotation and offset that matches points to reference points
    """
    tree = sp.cKDTree(np.column_stack((xref, yref)))
    seps, indices = tree.query(np.column_stack((x, y)))
    if np.all(seps < 0.01):
        # probably the same data
        mask = np.zeros_like(seps).astype("bool")
    else:
        sep_sigma = (seps - seps.mean()) / seps.std()
        mask = sep_sigma < 1.5
    return findBestRigidTransform(
        x[mask], y[mask], xref[indices][mask], yref[indices][mask]
    )


def measureSpots(ccd, thresh, use_hfd=True):
    """
    Detect spots in all windows of a ccd.
    """
    objects = np.concatenate(
        [findStars(wind, thresh, 3) for (wname, wind) in ccd.items()]
    )
    xpos = objects["x"]
    ypos = objects["y"]
    fwhm = objects["fwhm"]
    hfd = objects["hfd"]

    if use_hfd:
        mask = np.logical_and(np.isfinite(hfd).data, objects["peak"] < 20000)
        mask = mask & (fwhm <= 20)
        return xpos[mask], ypos[mask], hfd[mask], objects["cflux"][mask]
    else:
        mask = fwhm <= 20
        mask = np.logical_and(mask, objects["peak"] < 20000)
        return xpos[mask], ypos[mask], fwhm[mask], objects["cflux"][mask]


def separateBigSmallSpots(xpos, ypos, fwhm, flux):
    """
    Use K-Means clustering to seperate large and small starfield tube spots based on
    FWHM and flux.

    Parameters
    ----------
    xpos, ypos : np.ndarray
        spot locations
    fwhm : np.ndarray
        fwhm of spots
    flux : np.ndarray
        flux contained within spots

    Returns
    -------
    lx, ly, lfw, lfl, bx, by, bfw, bfl : np.ndarray
        positions, fwhm and flux of small and large spots, respectively
    """
    X = np.array((fwhm, flux)).T
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(X)
    size_ids = kmeans.predict(X)
    spot_mask0 = size_ids == 0
    spot_mask1 = size_ids == 1
    if fwhm[spot_mask0].mean() > fwhm[spot_mask1].mean():
        big_spot_mask = spot_mask0
        little_spot_mask = spot_mask1
    else:
        big_spot_mask = spot_mask1
        little_spot_mask = spot_mask0
    return (
        xpos[little_spot_mask],
        ypos[little_spot_mask],
        fwhm[little_spot_mask],
        flux[little_spot_mask],
        xpos[big_spot_mask],
        ypos[big_spot_mask],
        fwhm[big_spot_mask],
        flux[big_spot_mask],
    )


def aligntool(args=None):
    """
    Measures alignment of HiperCAM CCDs w.r.t a reference image of the blue CCD.

    Arguments::

        ref     : (string)
           name of an MCCD file to use as reference , as produced
           by e.g. 'grab'

        rccd    : (string)
           CCD to use as reference

        thresh  : (float)
            threshold for object detection, in multiples of background RMS

        source  : (string) [hidden]
           Data source, five options for built-in sources:

               'hs' : HiPERCAM server
               'hl' : local HiPERCAM FITS file
               'us' : ULTRACAM server
               'ul' : local ULTRACAM .xml/.dat files
               'hf' : list of HiPERCAM hcm FITS-format files

           'hf' is used to look at sets of frames generated by 'grab' or
           converted from foreign data formats. The standard
           start-off default for ``source`` can be set using the
           environment variable HIPERCAM_DEFAULT_SOURCE. e.g. in bash
           :code:`export HIPERCAM_DEFAULT_SOURCE="us"` would ensure it
           always started with the ULTRACAM server by default. If
           unspecified, it defaults to 'hl'.

           Any instruments installed via plugins can be used as a source by
           specifying either ``<instrument>:local`` or ``<instrument>:server``.

        cmap : str [hidden]
           The matplotlib colour map to use. "Greys" gives the usual greyscale.
           "none" will give whatever the current default is. Many other choices:
           "viridis", "jet", "hot", "Oranges", etc. Enter an invalid one and
           the program will fail but return a huge list of possibles in the
           process. Note that some will not work well with the colours used
           to plot annotations.

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

        fwwidth : float [hidden]
           FWHM display plot width in inches
           (0 for default, and note that this will also set the height to its default)

        fwheight : float [hidden]
           FWHM display plot height in inches
           (0 for default, and note that this will also set the height to its default)

        fwxpos : int [hidden]
            FWHM display plot X position in pixels

        fwypos : int [hidden]
            FWHM display plot Y position in pixels

        memory : bool [hidden]
            stores the plot display dimensions (width, height) and position
            (xpos, ypos) at the end as defaults for next time.
            This helps if you are repeatedly re-starting the script.

        run     : (string) [if source == 's' or 'l']
           run number to access, e.g. 'run034'

        flist   : (string) [if source == 'f']
           name of file list

        first   : (int) [if source='s' or 'l']
           exposure number to start from. 1 = first frame; set = 0 to
           always try to get the most recent frame (if it has changed)

        twait   : (float) [if source == 's'; hidden]
           time to wait between attempts to find a new exposure, seconds.

        tmax    : (float) [if source == 's'; hidden]
           maximum time to wait between attempts to find a new exposure, seconds.

        pause : float [hidden]
           seconds to pause between frames (defaults to 0). This is in addition
           to any time taken per frame to read and display it, so is just to slow
           things not to set an exact rate.

        tupdate : float [hidden]
           target number of seconds between plot updates. This is to
           make the plots more responsive and get round issues with
           plots freezing until the next frame comes through which can
           be especially annoying with virtual desktops. Experiments
           with ULTRACAM@NTT indicated 0.2 to 0.6 secs or so works OK.
           Main impact seems to be on drpc CPU.

        ccd     : (string)
           CCD to measure aligment of.

        bias    : (string)
           Name of bias frame to subtract, 'none' to ignore.

        msub    : (bool)
           subtract the median from each window before scaling for the
           image display or not. This happens after any bias subtraction.

        iset    : (string) [single character]
           determines how the intensities are determined. There are three
           options: 'a' for automatic simply scales from the minimum to the
           maximum value found on a per CCD basis. 'd' for direct just takes
           two numbers from the user. 'p' for percentile dtermines levels
           based upon percentiles determined from the entire CCD on a per CCD
           basis.

        ilo     : (float) [if iset='d']
           lower intensity level

        ihi     : (float) [if iset='d']
           upper intensity level

        plo     : (float) [if iset='p']
           lower percentile level

        phi     : (float) [if iset='p']
           upper percentile level

        small_spots : (bool) [hidden]
            use the small spots for alignment instead of the big spots
    """

    command, args = cline.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:
        # register parameters
        cl.register("ref", Cline.LOCAL, Cline.PROMPT)
        cl.register("thresh", Cline.LOCAL, Cline.PROMPT)
        cl.register("small_spots", Cline.LOCAL, Cline.HIDE)
        cl.register("source", Cline.LOCAL, Cline.HIDE)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
        cl.register("imwidth", Cline.LOCAL, Cline.HIDE)
        cl.register("imheight", Cline.LOCAL, Cline.HIDE)
        cl.register("imxpos", Cline.LOCAL, Cline.HIDE)
        cl.register("imypos", Cline.LOCAL, Cline.HIDE)
        cl.register("fwwidth", Cline.LOCAL, Cline.HIDE)
        cl.register("fwheight", Cline.LOCAL, Cline.HIDE)
        cl.register("fwxpos", Cline.LOCAL, Cline.HIDE)
        cl.register("fwypos", Cline.LOCAL, Cline.HIDE)
        cl.register("memory", Cline.LOCAL, Cline.HIDE)
        cl.register("run", Cline.GLOBAL, Cline.PROMPT)
        cl.register("first", Cline.LOCAL, Cline.PROMPT)
        cl.register("twait", Cline.LOCAL, Cline.HIDE)
        cl.register("tmax", Cline.LOCAL, Cline.HIDE)
        cl.register("flist", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("rccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("pause", Cline.LOCAL, Cline.HIDE)
        cl.register("tupdate", Cline.LOCAL, Cline.HIDE)
        cl.register("bias", Cline.GLOBAL, Cline.PROMPT)
        cl.register("msub", Cline.GLOBAL, Cline.PROMPT)
        cl.register("iset", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ilo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ihi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("plo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("phi", Cline.LOCAL, Cline.PROMPT)
        cl.register("xlo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("xhi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ylo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("yhi", Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        reference = cl.get_value(
            "ref", "reference frame to align to", cline.Fname("hcam", hcam.HCAM)
        )

        ref_mccd = hcam.MCCD.read(reference)

        thresh = cl.get_value("thresh", "detection threshold (sigma)", 5.0)
        small_spots = cl.get_value(
            "small_spots", "use small spots for alignment", False
        )

        # get inputs
        default_source = os.environ.get("HIPERCAM_DEFAULT_SOURCE", "hl")
        source = cl.get_value(
            "source",
            "data source [hs, hl, us, ul, hf, or <instrument>:local/<instrument>:server]",
            default_source,
            lvals=available_instruments(),
        )

        # set some flags
        server_or_local = source.endswith("s") or source.endswith("l")
        inst = source[0]

        # Some settings for the image plots
        cmap = cl.get_value(
            "cmap", "colour map to use ['none' for mpl default]", "Greys"
        )
        cmap = None if cmap == "none" else cmap

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

        fwwidth = cl.get_value(
            "fwwidth", "FWHM plot width [inches, 0 for default]", 0.0, 0.0
        )

        fwheight = cl.get_value(
            "fwheight", "FWHM plot height [inches, 0 for default]", 0.0, 0.0
        )

        fwxpos = cl.get_value(
            "fwxpos", "FWHM plot X position [screen pixels]", 0,
        )

        fwypos = cl.get_value(
            "fwypos", "FWHM plot Y position [screen pixels]", 0,
        )

        memory = cl.get_value(
            "memory", "remember final plot geometry for next time?", True
        )

        if server_or_local:
            resource = cl.get_value("run", "run name", "run005")
            first = cl.get_value("first", "first frame to plot", 1, 1)
            twait = cl.get_value(
                "twait", "time to wait for a new frame [secs]", 1.0, 0.0
            )
            tmax = cl.get_value(
                "tmax", "maximum time to wait for a new frame [secs]", 10.0, 0.0
            )
        else:
            # set inst = 'h' as only lists of HiPERCAM files are supported
            inst = "h"
            resource = cl.get_value(
                "flist", "file list", cline.Fname("files.lis", hcam.LIST)
            )
            first = 1

        flist = source.endswith("f")
        server = source.endswith("s")
        if inst == "u":
            instrument = "ULTRA"
        elif inst == "h":
            instrument = "HIPER"

        # get the labels and maximum dimensions and padding factors
        ccdinf = spooler.get_ccd_pars(source, resource)
        if len(ccdinf) > 1:
            data_cnam = cl.get_value("ccd", "CCD to plot alignment of", "1")
            ref_cnam = cl.get_value("rccd", "CCD to use as reference", "5")

        cl.set_default("pause", 0.0)
        pause = cl.get_value(
            "pause", "time delay to add between frame plots [secs]", 0.0, 0.0
        )

        tupdate = cl.get_value(
            "tupdate", "target time between plot updates [secs]", 0.2, 0.01
        )

        # bias frame (if any)
        bias = cl.get_value(
            "bias",
            "bias frame ['none' to ignore]",
            cline.Fname("bias", hcam.HCAM),
            ignore="none",
        )
        if bias is not None:
            # read the bias frame
            bias = hcam.MCCD.read(bias)

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

        # region to plot
        nxtot, nytot, nxpad, nypad = ccdinf[data_cnam]
        xmin, xmax = float(-nxpad), float(nxtot + nxpad + 1)
        ymin, ymax = float(-nypad), float(nytot + nypad + 1)

        xlo = cl.get_value("xlo", "left-hand X value", xmin, xmin, xmax)
        xhi = cl.get_value("xhi", "right-hand X value", xmax, xmin, xmax)
        ylo = cl.get_value("ylo", "lower Y value", ymin, ymin, ymax)
        yhi = cl.get_value("yhi", "upper Y value", ymax, ymin, ymax)

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff

    # Calibrate and measure spots in reference image
    # Since this doesn't change, we only need to do it once outside the loop
    ref_ccd = ref_mccd[ref_cnam]
    if msub:
        for wind in ref_ccd.values():
            wind -= wind.median()
    try:
        xpos, ypos, fwhm, rflux = measureSpots(ref_ccd, thresh)
        lx, ly, _lf, _lp, bx, by, _bf, _bp = separateBigSmallSpots(
            xpos, ypos, fwhm, rflux
        )
        if not small_spots:
            # use big spots
            xref, yref = bx, by
        else:
            xref, yref = lx, ly
    except Exception as err:
        raise ValueError(f"Failed to find spots in reference CCD {ref_cnam}") from err

    # Initialise the plot managers
    imanager = ImageManager(
        ref_ccd=ref_ccd,
        iset=iset,
        ilim=(plo, phi) if iset == "p" else (ilo, ihi),
        region=(xlo, xhi, ylo, yhi),
        cmap=cmap,
        figsize=(imwidth, imheight),
        position=(imxpos, imypos),
        memory_args=(command, args) if memory else None,
    )
    fwhmmanager = FwhmManager(
        fix_fwhm_scale=False,
        figsize=(fwwidth, fwheight),
        position=(fwxpos, fwypos),
        memory_args=(command, args) if memory else None,
    )

    # Start the main loop
    total_time = 0.0
    with  (
        suppress(KeyboardInterrupt),  # catch Ctrl-C to exit the loop cleanly
        spooler.data_source(source, resource, first) as spool
    ):
        for n, mccd in enumerate(spool):

            if server_or_local:
                # Handle the waiting game ...
                # Includes some awkward stuff involving updating on a cycle
                # faster than twait to make the plots more responsive,
                # if twait is long.
                updaters = [imanager.refresh, fwhmmanager.refresh]
                give_up, try_again, total_time = spooler.hang_about(
                    mccd, twait, tmax, total_time, updaters, tupdate
                )

                if give_up:
                    print("alignment tool stopped")
                    break
                elif try_again:
                    continue

            # first time setup
            if n == 0 and bias is not None:
                # crop the bias on the first frame only
                bias = crop_calib_frame_to(mccd, bias, "bias")

            # calibration
            ccd = mccd[data_cnam]
            if bias is not None:
                ccd -= bias[data_cnam]
            if msub:
                for wind in ccd.values():
                    wind -= wind.median()
                for wind in ref_ccd.values():
                    wind -= wind.median()

            # Update the image plot
            imanager.update(ccd)

            # Now print progress indicators
            nframe = mccd.head.get("NFRAME", n)
            tstamp = Time(mccd.head["TIMSTAMP"], format="isot", precision=3)
            tok = "ok" if mccd.head.get("GOODTIME", True) else "nok"
            print(f"{nframe}, utc= {tstamp.iso} ({tok}), ", end="")

            # Get statistics to print out
            message = f"exp= {mccd.head['EXPTIME']:.4f}"
            for cnam, data in {'reference': ref_ccd, 'data': ccd}.items():
                message += f", CCD {cnam}: "

                if not data.is_data():
                    # skipped due to nskip > 0, no data in this frame
                    message += "none"
                    continue

                # add the image limits
                vmin, vmax = imanager.get_ilim(cnam)
                message += f"{vmin:.1f} to {vmax:.1f}"

            # Print messages
            print(message)

            # time for measurement of spots!
            try:
                xpos, ypos, fwhm, flux = measureSpots(ccd, thresh)
                lx, ly, lf, lp, bx, by, bf, bp = separateBigSmallSpots(
                    xpos, ypos, fwhm, flux
                )
                if not small_spots:
                    # use big spots
                    x, y, f, flux = bx, by, bf, bp
                else:
                    x, y, f, flux = lx, ly, lf, lp
            except Exception as err:
                print(f"Failed to find spots in frame: {err}")
                continue

            # We found spots in both images, so let's plot the FWHM and offsets
            fwhm_data = {
                "x": x,
                "y": y,
                "xref": xref,
                "yref": yref,
                "fwhm": f,
                "flux": flux,
            }
            fwhmmanager.update(fwhm_data)

            if pause > 0.0:
                # pause between frames
                if pause > 2 * tupdate:
                    # run the updaters if pause is a little long.
                    # designed to make the plot more responsive
                    nupdate = int(pause / tupdate)
                    for _ in range(nupdate):
                        time.sleep(pause / nupdate)
                        if imanager is not None:
                            imanager.refresh()
                        if fwhmmanager is not None:
                            fwhmmanager.refresh()
                else:
                    # just pause because the updates will come
                    time.sleep(pause)

    # At then end the plots will close automatically,
    # so block this by asking the user to confirm
    input("\nHit carriage return to close the plot(s): ")


class ImageManager(hcam.mpl.AnimationManager):
    """Class to control the image animation.

    Initialising it creates and opens the empty figure.

    After that, call ``update`` once per frame to update the artists,
    and ``refresh`` as often as needed to keep the window responsive.
    """

    def __init__(
        self,
        ref_ccd,
        iset,
        ilim,
        region,
        cmap,
        figsize,
        position,
        memory_args,
    ):
        """Initialises and manages the plot.

        Arguments:

        ref_ccd : :class:`CCD` | None
            The calibrated reference CCD.

        ncols : int
            The number of columns in the grid of subplots

        iset : str
            How to set the intensity scale for the images:
                'a' for automatic min/max range
                'p' for percentiles (set using `ilim`)
                'd' for direct values (set using `ilim`)

        ilim : tuple[float, float] | None
            Lower and upper limits for the intensity scaling.

        region : tuple[int, int, int, int] | None
            The region to plot, values are (xmin, xmax, ymin, ymax)
            If None, the entire CCD is displayed.
            If iset='p', only this region is used to compute the percentiles for
            intensity scaling.

        cmap : str  | None
            Matplotlib colour map. None gives default

        figsize : tuple[float, float] | None
            The width and height of the figure, in inches.
            If None, the default size is used.

        position : tuple[int, int] | None
            The position of the figure window on the screen, in pixels.
            If None, the default position is used.

        memory_args : tuple[str, list] | None
            If not None, the first element is the command name and the second
            is the list of arguments. The window geometry will be saved to the
            command's Cline file on close, and restored on next run.

        """
        # create empty figure and axes
        self.fig, self.axs = hcam.mpl.create_figure(
            cnams = ['reference', 'data'],
            ncols=1,
            region=region,
            title="CCD display",
            figsize=figsize,
            position=position,
        )

        # initialise the animation manager with the figure
        super().__init__(self.fig)

        # store the reference CCD to display
        self.ref_ccd = ref_ccd

        # inputs for the image plots
        self.iset = iset
        self.ilim = ilim
        self.region = region
        self.cmap = cmap

        # save window geometry on close
        if memory_args:
            def save_geometry(event):
                imwidth, imheight, imxpos, imypos = hcam.mpl.store_window(
                    self.fig, memory_args[0], memory_args[1],
                )
                print(
                    "Image display memory updated [memory=True]:",
                    f"imwidth={imwidth:.1f}, imheight={imheight:.1f}, "
                    f"imxpos={imxpos}, imypos={imypos}"
                )
            self.fig.canvas.mpl_connect("close_event", save_geometry)

        # open the plot
        self.show()

    def draw_artists(self):
        """Draw all of the animated artists."""
        # Draw the CCD artists
        for cnam in self.artists:
            # each CCD has a dictionary of artists (returned by pCcd)
            for artists in self.artists[cnam].values():
                for artist in artists:
                    self.fig.draw_artist(artist)

    def update(self, ccd):
        """Create or update the artists with new data.

        Arguments:

        ccd : :class:`CCD` | None
            The frame to display.

        """
        # Keep track of if new artists are created
        new_artists = False

        # First, plot the reference CCD (only done once)
        if 'reference' not in self.artists:
            new_artists = True

            # get scale limits
            vmin, vmax = hcam.mpl.get_scale_limits(
                self.ref_ccd,
                mode=self.iset,
                limits=self.ilim,
                region=self.region,
            )

            # plot the reference CCD
            self.artists['reference'] = hcam.mpl.pCcd(
                ax=self.axs['reference'],
                ccd=self.ref_ccd,
                vmin=vmin,
                vmax=vmax,
                cmap=self.cmap,
                tlabel='reference',
                animated=True,
                artists=None,
            )

        # Update CCD artists for the new data frame
        if ccd.is_data():
            # get scale limits
            vmin, vmax = hcam.mpl.get_scale_limits(
                ccd=ccd,
                mode=self.iset,
                limits=self.ilim,
                region=self.region,
            )

            # plot the CCD
            # pCcd will update the existing artists if they exist,
            # otherwise it will create new ones
            if "data" not in self.artists:
                new_artists = True
            self.artists['data'] = hcam.mpl.pCcd(
                ax=self.axs['data'],
                ccd=ccd,
                vmin=vmin,
                vmax=vmax,
                cmap=self.cmap,
                tlabel='data',
                animated=True,
                artists=self.artists['data'] if 'data' in self.artists else None,
            )

        if new_artists:
            # In this case we are setting up for the first time
            # We need to set the layout and re-draw to avoid irritating distorted images
            self.fig.tight_layout()
            self.fig.canvas.draw()

        # Refresh the display
        self.refresh(redraw=True)

    def get_ilim(self, cnam):
        """Returns the current intensity limits for the specified CCD."""
        if self.artists[cnam] is None:
            return None
        imgs = self.axs[cnam].get_images()
        return imgs[0].get_clim()  # for the first window, but should be the same


class FwhmManager(hcam.mpl.AnimationManager):
    """Class to start and control the FWHM history animation."""

    def __init__(
        self,
        fix_fwhm_scale=False,
        bins=20,
        figsize=None,
        position=None,
        memory_args=None,
    ):
        """Initialises the FWHM history plot.

        Arguments:

        fix_fwhm_scale : bool
            If True, the FWHM axes are fixed to a scale of 1.9 to 2.8 arcsec.
            If False, the axes are automatically scaled to the data.

        bins : int
            The number of bins to use for the FWHM histograms.

        figsize : tuple[float, float] | None
            The width and height of the figure, in inches.
            If None, the default size is used.

        position : tuple[int, int] | None
            The position of the figure window on the screen, in pixels.
            If None, the default position is used.

        memory_args : tuple[str, list] | None
            If not None, the first element is the command name and the second
            is the list of arguments. The window geometry will be saved to the
            command's Cline file on close, and restored on next run.

        """
        # create empty figure
        self.fig = plt.figure()
        self.fig.canvas.manager.set_window_title("FWHM plot")
        if figsize is not None and figsize[0] > 0 and figsize[1] > 0:
            self.fig.set_size_inches(figsize)
        if position is not None:
            hcam.mpl.set_window_position(self.fig, position[0], position[1])

        # initialise the animation manager with the figure
        super().__init__(self.fig)

        # definitions for the axes
        self.left = 0.05
        self.bottom = 0.1
        self.width = 0.38
        self.height = 0.5
        self.right = self.left + self.width + 0.02
        self.left_h = self.right + self.width + 0.02
        self.bottom_h = self.bottom + self.height + 0.02

        rect_im = [self.right, self.bottom, self.width, self.height]
        rect_xcut = [self.right, self.bottom_h, self.width, 0.2]
        rect_ycut = [self.left_h, self.bottom, 0.1, self.height]
        rect_pos = [self.left, self.bottom, self.width, self.height]
        rect_corr = [self.left, self.bottom_h + 0.1, self.width, 0.1]

        self.axs = {}
        self.axs['im'] = self.fig.add_axes(rect_im, animated=True)
        self.axs['x'] = self.fig.add_axes(rect_xcut, animated=True)
        self.axs['y'] = self.fig.add_axes(rect_ycut, animated=True)
        self.axs['pos'] = self.fig.add_axes(
            rect_pos, sharex=self.axs['im'], sharey=self.axs['im'], animated=True
        )
        self.axs['corr'] = self.fig.add_axes(rect_corr, animated=True)

        # params
        self.fix_fwhm_scale = fix_fwhm_scale
        self.bins = bins

        # plot formatting
        self.axs['y'].yaxis.set_major_formatter(NullFormatter())
        self.axs['x'].xaxis.set_major_formatter(NullFormatter())
        self.axs['x'].yaxis.set_label_position("right")
        self.axs['x'].yaxis.tick_right()
        self.axs['pos'].set_xlabel("x")
        self.axs['pos'].set_ylabel("y")
        self.axs['im'].set_xlabel("x")
        self.axs['x'].set_ylabel("FWHM")
        self.axs['y'].set_xlabel("FWHM")
        if self.fix_fwhm_scale:
            self.axs['y'].set_xlim(1.9, 2.8)
            self.axs['x'].set_ylim(1.9, 2.8)

        # save window geometry on close
        if memory_args:
            def save_geometry(event):
                fwwidth, fwheight, fwxpos, fwypos = hcam.mpl.store_window(
                    self.fig, memory_args[0], memory_args[1], prefix="fw",
                )
                print(
                    "FWHM display memory updated [memory=True]:",
                    f"fwwidth={fwwidth:.1f}, fwheight={fwheight:.1f}, "
                    f"fwxpos={fwxpos}, fwypos={fwypos}"
                )
            self.fig.canvas.mpl_connect("close_event", save_geometry)

        # open the plot
        self.show()

    def draw_artists(self):
        """Draw all of the animated artists."""
        # Draw the axes
        for ax in self.axs.values():
            self.fig.draw_artist(ax)

            # We also need to manually rescale the axes when animating
            if not self.fix_fwhm_scale:
                ax.relim()
                ax.autoscale_view()

        # Draw the artists
        for artist in self.artists.values():
            if artist is not None:
                self.fig.draw_artist(artist)

    def update(self, data):
        """Create or update the artists with new data."""
        xpos = data["x"]
        ypos = data["y"]
        xref = data["xref"]
        yref = data["yref"]
        fwhm = data["fwhm"]
        peak = data["flux"]

        # remove outliers
        mask = ~sigma_clip(peak).mask
        if len(mask) > 0 and len(mask) == len(xpos):
            xpos = xpos[mask]
            ypos = ypos[mask]
            fwhm = fwhm[mask]
            peak = peak[mask]

        # plot the FWHM as a function of position and flux
        if self.artists.get('x') is None:
            self.artists['x'] = self.axs['x'].plot(xpos, fwhm, "r.", animated=True)[0]
            self.artists['y'] = self.axs['y'].plot(fwhm, ypos, "r.", animated=True)[0]
        else:
            self.artists['x'].set_data(xpos, fwhm)
            self.artists['y'].set_data(fwhm, ypos)

        try:
            # remove correlation between flux and spot fwhm
            f = np.poly1d(np.polyfit(peak, fwhm, 1))
            meanF = np.mean(f(peak))
            fwhm_prime = fwhm - f(peak) + meanF

            # plot the corrected FWHM as a function of position and flux
            if self.artists.get('x_corr') is None:
                self.artists['x_corr'] = self.axs['x'].plot(
                    xpos, fwhm_prime, "b+", animated=True
                )[0]
                self.artists['y_corr'] = self.axs['y'].plot(
                    fwhm_prime, ypos, "b+", animated=True
                )[0]
            else:
                self.artists['x_corr'].set_data(xpos, fwhm_prime)
                self.artists['y_corr'].set_data(fwhm_prime, ypos)

            # plot the correlation between flux and spot fwhm
            idx = np.argsort(peak)
            if self.artists.get('corr') is None:
                self.artists['corr'] = self.axs['corr'].plot(
                    peak, fwhm, "r.", animated=True
                )[0]
                self.artists['corr_f'] = self.axs['corr'].plot(
                    peak[idx], f(peak[idx]), "b--", animated=True
                )[0]
            else:
                self.artists['corr'].set_data(peak, fwhm)
                self.artists['corr_f'].set_data(peak[idx], f(peak[idx]))

            # plot the image of the FWHM as a function of position
            xi = np.linspace(xpos.min(), xpos.max(), self.bins)
            yi = np.linspace(ypos.min(), ypos.max(), self.bins)
            xx, yy = np.meshgrid(xi, yi)
            func = interpolate.Rbf(xpos, ypos, fwhm_prime, function="quintic")
            image = func(xx, yy)
            if self.artists.get('im') is None:
                self.artists['im'] = self.axs['im'].pcolormesh(
                    xx, yy, image, cmap="magma", animated=True
                )
            else:
                self.artists['im'].set_array(image.ravel())

        except:
            # if the fitting or interpolation fails, just leave those plots empty
            pass

        # plot the positions of the spots in the data and reference images
        if self.artists.get('pos') is None:
            self.artists['pos'] = self.axs['pos'].plot(
                xpos, ypos, "r.", animated=True
            )[0]
            self.artists['pos_ref'] = self.axs['pos'].plot(
                xref, yref, "b+", alpha=0.5, animated=True
            )[0]
        else:
            self.artists['pos'].set_data(xpos, ypos)
            self.artists['pos_ref'].set_data(xref, yref)

        # add text box with the alignment information
        theta, offset = matchSpots(xpos, ypos, xref, yref)
        xoff, yoff = offset
        textstr = f"Θ = {theta:.2f}; xoff={xoff:.1f}, yoff={yoff:.1f}"
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
        if self.artists.get('pos_text') is None:
            self.artists['pos_text'] = self.axs['pos'].text(
                self.left,
                1.02,
                textstr,
                horizontalalignment="left",
                verticalalignment="bottom",
                bbox=props,
                fontsize=10,
                transform=self.axs['pos'].transAxes,
                animated=True
            )
        else:
            self.artists['pos_text'].set_text(textstr)

        # Finally, refresh the display
        self.refresh(redraw=True)
