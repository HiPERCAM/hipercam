#!/usr/bin/env python
from __future__ import division, print_function

import os
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interpolate
import scipy.spatial as sp
from astropy.stats import sigma_clip
from astropy.time import Time
from matplotlib.ticker import NullFormatter
from sklearn.cluster import KMeans
from trm import cline
from trm.cline import Cline

import hipercam as hcam
from hipercam import spooler
from hipercam.ccd import crop_calib_frame_to
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

        device  : (string) [hidden]
          Plot device. PGPLOT is used so this should be a PGPLOT-style name,
          e.g. '/xs', '1/xs' etc. At the moment only ones ending /xs are
          supported.

        cmap : str [hidden]
           The matplotlib colour map to use. "Greys" gives the usual greyscale.
           "none" will give whatever the current default is. Many other choices:
           "viridis", "jet", "hot", "Oranges", etc. Enter an invalid one and
           the program will fail but return a huge list of possibles in the
           process. Note that some will not work well with the colours used
           to plot annotations.

        width   : (float) [hidden]
           plot width (inches). Set = 0 to let the program choose.

        height  : (float) [hidden]
           plot height (inches). Set = 0 to let the program choose. BOTH width
           AND height must be non-zero to have any effect

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

        pause   : (float) [hidden]
           seconds to pause between frames (defaults to 0)

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
        cl.register("device", Cline.LOCAL, Cline.HIDE)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)
        cl.register("run", Cline.GLOBAL, Cline.PROMPT)
        cl.register("first", Cline.LOCAL, Cline.PROMPT)
        cl.register("twait", Cline.LOCAL, Cline.HIDE)
        cl.register("tmax", Cline.LOCAL, Cline.HIDE)
        cl.register("flist", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("rccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("pause", Cline.LOCAL, Cline.HIDE)
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
        width = cl.get_value("width", "plot width (inches)", 0.0)
        height = cl.get_value("height", "plot height (inches)", 0.0)

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

    # Initialise the plot managers to None, so we can create them on the first frame
    imanager = None
    fwhmmanager = None

    # Start the main loop
    total_time = 0.0
    with spooler.data_source(source, resource, first) as spool:
        for n, mccd in enumerate(spool):

            if server_or_local:
                # Handle the waiting game ...
                give_up, try_again, total_time = spooler.hang_about(
                    mccd, twait, tmax, total_time
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

            # get the right classes for the CCD names
            ccd_data = {'reference': ref_mccd[ref_cnam], 'data': mccd[data_cnam]}

            # calibration
            if bias is not None:
                ccd_data['data'] -= bias[data_cnam]
            if msub:
                for wind in ccd_data['data'].values():
                    wind -= wind.median()
                for wind in ccd_data['reference'].values():
                    wind -= wind.median()

            # We're ready to start plotting
            if imanager is None:
                # First time through, create the image plot manager
                imanager = ImageManager(
                    ccd_data=ccd_data,
                    iset=iset,
                    ilim=(plo, phi) if iset == "p" else (ilo, ihi),
                    region=(xlo, xhi, ylo, yhi),
                    cmap=cmap,
                    figsize=(width, height),
                )

                # start live display
                plt.show(block=False)
                plt.pause(0.1)

            # Update the image plot (done on every loop)
            imanager.update(ccd_data)

            # Now print progress indicators
            nframe = mccd.head.get("NFRAME", n)
            tstamp = Time(mccd.head["TIMSTAMP"], format="isot", precision=3)
            tok = "ok" if mccd.head.get("GOODTIME", True) else "nok"
            print(f"{nframe}, utc= {tstamp.iso} ({tok}), ", end="")

            # Get statistics to print out
            message = f"exp= {mccd.head['EXPTIME']:.4f}"
            for cnam in ccd_data:
                message += f", CCD {cnam}: "
                if ccd_data[cnam] is None:
                    # skipped due to nskip > 0, no data in this frame
                    message += "none"
                    continue

                # accumulate string of image scalings
                imgs = imanager.caxes[cnam].get_images()
                vmin, vmax = imgs[0].get_clim()
                message += f"{vmin:.1f} to {vmax:.1f}"

            print(message)

            # time for measurement of spots!
            try:
                xpos, ypos, fwhm, flux = measureSpots(ccd_data['data'], thresh)
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

            # also measure spots in reference image
            try:
                xpos, ypos, fwhm, rflux = measureSpots(ccd_data['reference'], thresh)
                lx, ly, lf, lp, bx, by, bf, bp = separateBigSmallSpots(
                    xpos, ypos, fwhm, rflux
                )
                if not small_spots:
                    # use big spots
                    xref, yref = bx, by
                else:
                    xref, yref = lx, ly
            except Exception as err:
                print(f"Failed to find spots in reference: {err}")
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
            if fwhmmanager is None:
                # First time through, create the FWHM plot manager
                fwhmmanager = FwhmManager(
                    fwhm_data,
                    fix_fwhm_scale=False,
                )
                plt.show(block=False)
                plt.pause(0.1)

            # Update the FWHM plot (done on every loop)
            fwhmmanager.update(fwhm_data)

            if pause > 0.0:
                time.sleep(pause)

    # At then end the plots will close automatically,
    # so block this by asking the user to confirm
    input("\nHit carriage return to close the plot(s): ")


class ImageManager:
    """Class to control the image animation"""

    def __init__(
        self,
        ccd_data,
        iset,
        ilim,
        region,
        cmap,
        figsize,
    ):
        """Initialises and manages the plot.

        Arguments:

        ccd_data : dict of {'reference': :class:`CCD`|None, 'data': :class:`CCD`|None}
            Contains the calibrated CCD, or None if the CCD is being skipped,
            for both the reference and data frames.
            Used to create the artists; later frames go through update_images.

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

        """
        # create empty figure and axes
        self.fig, self.caxes = hcam.mpl.create_figure(
            cnams = ['reference', 'data'],
            ncols=1,
            region=region,
            title="CCD display",
            figsize=figsize,
            animated=True,
        )

        # grab the background on every draw
        self._bg = None
        self.fig.canvas.mpl_connect("draw_event", self.on_draw)

        # inputs for the image plots
        self.iset = iset
        self.ilim = ilim
        self.region = region
        self.cmap = cmap

        # now actually create and update the artists
        self.img_artists = {'reference': None, 'data': None}
        self.update(ccd_data)

        # adjust the layout to prevent overlap
        # (can only do after the artists are created)
        self.fig.tight_layout()

    def on_draw(self, event):
        """Callback to register with 'draw_event'."""
        if event is not None and event.canvas != self.fig.canvas:
            raise RuntimeError
        self._bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        self._draw_animated()

    def _draw_animated(self):
        """Draw all of the animated artists."""
        for cnam in self.caxes:
            # Draw the axes before the content if you want to see anything!
            self.fig.draw_artist(self.caxes[cnam])

            if self.img_artists[cnam] is not None:
                # 'ccd' has values which are all lists of artists
                for artists in self.img_artists[cnam].values():
                    for artist in artists:
                        self.fig.draw_artist(artist)

    def refresh(self):
        """Redraws the animated artists without changing any data.

        Safe to call repeatedly while waiting for a frame to arrive,
        to keep the window responsive.
        """
        if self._bg is None:
            self.on_draw(None)
        else:
            # restore the background
            self.fig.canvas.restore_region(self._bg)
            self._draw_animated()
            self.fig.canvas.blit(self.fig.bbox)
        self.fig.canvas.flush_events()

    def update(self, ccd_data):
        """Creates or updates the artists for each CCD with new data."""
        new_artists = {}

        # update CCD artists
        for cnam in ccd_data:
            if ccd_data[cnam] is None:
                # no new frame, just pass old artists through
                new_artists[cnam] = self.img_artists[cnam]
            else:
                # create/update the animated artists
                axes = self.caxes[cnam]
                ccd = ccd_data[cnam]
                old_artists = self.img_artists[cnam]
                new_artists[cnam] = {}

                # get scale limits
                vmin, vmax = hcam.mpl.get_scale_limits(
                    ccd,
                    mode=self.iset,
                    limits=self.ilim,
                    region=self.region,
                )

                # plot the CCD
                new_artists[cnam] = hcam.mpl.pCcd(
                    axes=axes,
                    ccd=ccd,
                    vmin=vmin,
                    vmax=vmax,
                    cmap=self.cmap,
                    tlabel=cnam,
                    animated=True,
                    artists=old_artists,
                )

                # In this case we are setting up for the first time
                # need to re-draw to avoid irritating distorted image
                if old_artists is None:
                    self.fig.canvas.draw()

        # Finally, update the artists and refresh the display
        self.img_artists = new_artists
        self.refresh()


class FwhmManager:
    """Class to start and control the FWHM history animation"""

    def __init__(self, data, fix_fwhm_scale=False, bins=20):
        """Initialises the FWHM history plot."""
        # create empty figure
        self.fig = plt.figure()
        self.fig.canvas.manager.set_window_title("FWHM plot")

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

        self.axes = {}
        self.axes['im'] = self.fig.add_axes(rect_im, animated=True)
        self.axes['x'] = self.fig.add_axes(rect_xcut, animated=True)
        self.axes['y'] = self.fig.add_axes(rect_ycut, animated=True)
        self.axes['pos'] = self.fig.add_axes(
            rect_pos, sharex=self.axes['im'], sharey=self.axes['im'], animated=True
        )
        self.axes['corr'] = self.fig.add_axes(rect_corr, animated=True)

        # grab the background on every draw
        self._bg = None
        self.fig.canvas.mpl_connect("draw_event", self.on_draw)

        # params
        self.fix_fwhm_scale = fix_fwhm_scale
        self.bins = bins

        # plot formatting
        self.axes['y'].yaxis.set_major_formatter(NullFormatter())
        self.axes['x'].xaxis.set_major_formatter(NullFormatter())
        self.axes['x'].yaxis.set_label_position("right")
        self.axes['x'].yaxis.tick_right()
        self.axes['pos'].set_xlabel("x")
        self.axes['pos'].set_ylabel("y")
        self.axes['im'].set_xlabel("x")
        self.axes['x'].set_ylabel("FWHM")
        self.axes['y'].set_xlabel("FWHM")
        if self.fix_fwhm_scale:
            self.axes['y'].set_xlim(1.9, 2.8)
            self.axes['x'].set_ylim(1.9, 2.8)

        # plot artists
        self.artists = {}
        self.update(data)

    def on_draw(self, event):
        """Callback to register with 'draw_event'."""
        if event is not None and event.canvas != self.fig.canvas:
            raise RuntimeError
        self._bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        self._draw_animated()

    def _draw_animated(self):
        """Draw all of the animated artists."""
        # Draw the axes before the content if you want to see anything!
        for ax in self.axes.values():
            self.fig.draw_artist(ax)

            # We also need to manually rescale the axes when animating
            if not self.fix_fwhm_scale:
                ax.relim()
                ax.autoscale_view()

        # Now draw the artists
        for artist in self.artists.values():
            if artist is not None:
                self.fig.draw_artist(artist)

    def refresh(self):
        """Redraws the animated artists without changing any data.

        Safe to call repeatedly while waiting for a frame to arrive,
        to keep the window responsive.
        """
        if self._bg is None:
            self.on_draw(None)
        else:
            # restore the background
            self.fig.canvas.restore_region(self._bg)
            self._draw_animated()
            self.fig.canvas.blit(self.fig.bbox)
        self.fig.canvas.flush_events()

    def update(self, data):
        """Displays the latest FWHM measurements."""
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
            self.artists['x'] = self.axes['x'].plot(xpos, fwhm, "r.", animated=True)[0]
            self.artists['y'] = self.axes['y'].plot(fwhm, ypos, "r.", animated=True)[0]
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
                self.artists['x_corr'] = self.axes['x'].plot(
                    xpos, fwhm_prime, "b+", animated=True
                )[0]
                self.artists['y_corr'] = self.axes['y'].plot(
                    fwhm_prime, ypos, "b+", animated=True
                )[0]
            else:
                self.artists['x_corr'].set_data(xpos, fwhm_prime)
                self.artists['y_corr'].set_data(fwhm_prime, ypos)

            # plot the correlation between flux and spot fwhm
            idx = np.argsort(peak)
            if self.artists.get('corr') is None:
                self.artists['corr'] = self.axes['corr'].plot(
                    peak, fwhm, "r.", animated=True
                )[0]
                self.artists['corr_f'] = self.axes['corr'].plot(
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
                self.artists['im'] = self.axes['im'].pcolormesh(
                    xx, yy, image, cmap="magma", animated=True
                )
            else:
                self.artists['im'].set_array(image.ravel())

        except:
            # if the fitting or interpolation fails, just leave those plots empty
            pass

        # plot the positions of the spots in the data and reference images
        if self.artists.get('pos') is None:
            self.artists['pos'] = self.axes['pos'].plot(
                xpos, ypos, "r.", animated=True
            )[0]
            self.artists['pos_ref'] = self.axes['pos'].plot(
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
            self.artists['pos_text'] = self.axes['pos'].text(
                self.left,
                1.02,
                textstr,
                horizontalalignment="left",
                verticalalignment="bottom",
                bbox=props,
                fontsize=10,
                transform=self.axes['pos'].transAxes,
                animated=True
            )
        else:
            self.artists['pos_text'].set_text(textstr)

        self.refresh()
