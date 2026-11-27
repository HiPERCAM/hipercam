"""Benchmark image plotting using PGPLOT vs MPL."""
import os
import sys
import time

import numpy as np
from astropy.time import Time
from matplotlib import pyplot as plt
from trm import cline
from trm.pgplot import pgenv, pglab, pgpanl, pgpap, pgsch, pgsci, pgsubp

import hipercam as hcam
from hipercam import spooler
from hipercam.cline import Cline
from hipercam.instruments import available_instruments


class ImageManager(hcam.mpl.AnimationManager):
    """A minimal version of the mpl plotting manager."""

    def __init__(self,
        cnams,
        ncols,
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

        cnams : list of str
            The names of the CCDs to plot.

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
            cnams=cnams,
            ncols=ncols,
            region=region,
            title="CCD display",
            figsize=figsize,
            position=position,
        )

        # initialise the animation manager with the figure
        super().__init__(self.fig)

        # inputs for the image plots
        self.cnams = cnams
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

    def update(self, mccd):
        """Create or update the artists with new data."""
        # Update the artists for the CCDs
        self._update_images(mccd)

        # Refresh the display
        self.refresh(redraw=True)

    def _update_images(self, mccd):
        """Create or update the CCD artists."""
        # Keep track of if new artists are created
        new_artists = False

        # Update CCD artists
        for cnam in self.cnams:
            if not mccd[cnam].is_data():
                # no new frame, just keep any existing artists
                continue

            # get scale limits
            vmin, vmax = hcam.mpl.get_scale_limits(
                mccd[cnam],
                mode=self.iset,
                limits=self.ilim,
                region=self.region,
            )

            # plot the CCD
            # pCcd will update the existing artists if they exist,
            # otherwise it will create new ones
            if cnam not in self.artists:
                new_artists = True
            self.artists[cnam] = hcam.mpl.pCcd(
                ax=self.axs[cnam],
                ccd=mccd[cnam],
                vmin=vmin,
                vmax=vmax,
                cmap=self.cmap,
                tlabel=f"CCD {cnam}",
                animated=True,
                artists=self.artists[cnam] if cnam in self.artists else None,
            )

        if new_artists:
            # In this case we are setting up for the first time
            # We need to set the layout and re-draw to avoid irritating distorted images
            self.fig.tight_layout()
            self.fig.canvas.draw()


def setup_implot(device, iwidth, iheight, ccds, nx, xlo, xhi, ylo, yhi):
    """Perform initial setup of the image plot for reduction outputs."""
    device = hcam.pgp.Device(device)
    if iwidth > 0 and iheight > 0:
        pgpap(iwidth, iheight / iwidth)

    # set up panels and axes
    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    # slice up viewport
    pgsubp(nx, ny)

    # plot axes, labels, titles. Happens once only
    for cnam in ccds:
        pgsci(hcam.pgp.Params["axis.ci"])
        pgsch(hcam.pgp.Params["axis.number.ch"])
        pgenv(xlo, xhi, ylo, yhi, 1, 0)
        pgsci(hcam.pgp.Params["axis.label.ci"])
        pgsch(hcam.pgp.Params["axis.label.ch"])
        pglab("X", "Y", "CCD {:s}".format(cnam))

    return device


def update_implot(
    mccd,
    device,
    ccds,
    nx,
    iset,
    plo,
    phi,
    ilo,
    ihi,
    xlo,
    xhi,
    ylo,
    yhi,
):
    """Plot updater routine.

    Arguments::

        mccd : MCCD
            current debiassed, flat-fielded frame. In the case of frame
            groups, this will be the last one added.

        device : hipercam.pgp.Device
            device to use for the plot

        ccds : list
            list of CCD labels to be plotted

        nx : int
            Number of panels in X for image display

        iset : character
            'a', 'd', or 'p' to define image display levels

        plo : float [only matters if iset =='p']
            Lower percentile level

        phi : float [only matters if iset =='p']
            Upper percentile level

        ilo : float [only matters if iset =='d']
            Lower intensity level

        ihi : float [only matters if iset =='d']
            Upper intensity level

        xlo : int
            left X-limit for determination of percentiles

        xhi : int
            right X-limit for determination of percentiles

        ylo : int
            lower Y-limit for determination of percentiles

        yhi : int
            upper Y-limit for determination of percentiles

    """
    # select the plot device
    device.select()

    # display the CCDs chosen
    for nc, cnam in enumerate(ccds):
        ccd = mccd[cnam]

        if ccd.is_data():
            # this should be data as opposed to a blank frame
            # between data frames that occur with nskip > 0

            # set to the correct panel and then plot CCD
            ix = (nc % nx) + 1
            iy = nc // nx + 1
            pgpanl(ix, iy)
            hcam.pgp.pCcd(
                ccd, iset, plo, phi, ilo, ihi, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi
            )

def main(mode, args=None):
    command, args = cline.script_args()

    # get the inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("source", Cline.GLOBAL, Cline.HIDE)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
        cl.register("imwidth", Cline.LOCAL, Cline.HIDE)
        cl.register("imheight", Cline.LOCAL, Cline.HIDE)
        cl.register("imxpos", Cline.LOCAL, Cline.HIDE)
        cl.register("imypos", Cline.LOCAL, Cline.HIDE)
        cl.register("memory", Cline.LOCAL, Cline.HIDE)
        cl.register("run", Cline.GLOBAL, Cline.PROMPT)
        cl.register("first", Cline.LOCAL, Cline.PROMPT)
        cl.register("last", Cline.LOCAL, Cline.HIDE)
        cl.register("twait", Cline.LOCAL, Cline.HIDE)
        cl.register("tmax", Cline.LOCAL, Cline.HIDE)
        cl.register("flist", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("tupdate", Cline.LOCAL, Cline.HIDE)
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
        default_source = os.environ.get('HIPERCAM_DEFAULT_SOURCE','hl')
        source = cl.get_value(
            "source",
            "data source [hs, hl, us, ul, hf, or <instrument>:local/<instrument>:server]",
            default_source,
            lvals=available_instruments(),
        )

        # set some flags
        server_or_local = source.endswith("s") or source.endswith("l")

        if server_or_local:
            resource = cl.get_value("run", "run name", "run005")
            if source == "hs":
                # HiPERCAM doesn't like first=0, so allow -ve values
                first = cl.get_value("first", "first frame to plot", 1)
            else:
                first = cl.get_value("first", "first frame to plot", 1, 0)

            cl.set_default("last", 0)
            last = cl.get_value("last", "last frame to reduce", 0, 0)
            if last and last < first:
                print("Cannot set last < first unless last == 0")
                print("*** reduce aborted")
                exit(1)

            twait = cl.get_value(
                "twait", "time to wait for a new frame [secs]", 1.0, 0.0
            )
            tmax = cl.get_value(
                "tmax", "maximum time to wait for a new frame [secs]", 10.0, 0.0
            )

        else:
            resource = cl.get_value(
                "flist", "file list", cline.Fname("files.lis", hcam.LIST)
            )
            first = 1
            last = 0

        # Some settings for the image plots
        cmap = cl.get_value(
            "cmap", "colour map to use ['none' for mpl default]", "Greys"
        )

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

        # define the panel grid. first get the labels and maximum dimensions
        ccdinf = spooler.get_ccd_pars(source, resource)

        nxdef = cl.get_default("nx", 3)

        if len(ccdinf) > 1:
            ccd = cl.get_value("ccd", "CCD(s) to plot [0 for all]", "0")
            if ccd == "0":
                ccds = list(ccdinf.keys())
            else:
                ccds = ccd.split()
                check = set(ccdinf.keys())
                if not set(ccds) <= check:
                    raise hcam.HipercamError("At least one invalid CCD label supplied")

            if len(ccds) > 1:
                nxdef = min(len(ccds), nxdef)
                cl.set_default("nx", nxdef)
                nx = cl.get_value("nx", "number of panels in X", 3, 1)
            else:
                nx = 1
        else:
            nx = 1
            ccds = list(ccdinf.keys())

        tupdate = cl.get_value(
            "tupdate", "target time between plot updates [secs]", 0.2, 0.01
        )

        iset = cl.get_value(
            "iset",
            "set intensity a(utomatically)," " d(irectly) or with p(ercentiles)?",
            "a",
            lvals=["a", "d", "p"],
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
            nxtot, nytot, nxpad, nypad = ccdinf[cnam]
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

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff

    if mode == "none":
        mpl = False
        pgp = False
    elif mode == "mpl":
        mpl = True
        pgp = False
    elif mode == "pgp":
        mpl = False
        pgp = True
    else:
        raise ValueError("Invalid mode")

    if mpl:
        # Initialise the plot manager
        imanager = ImageManager(
            cnams=ccds,
            ncols=nx,
            iset=iset,
            ilim=(plo, phi) if iset == "p" else (ilo, ihi),
            region=(xlo, xhi, ylo, yhi),
            cmap=cmap,
            figsize=(imwidth, imheight),
            position=(imxpos, imypos),
            memory_args=(command, args) if memory else None,
        )
        imanager.show()

    if pgp:
        # Set up the plot devices and buffers
        imdevice = setup_implot("/xs", imwidth, imheight, ccds, nx, xlo, xhi, ylo, yhi)

    # timing setup
    start_time = time.time()
    loop_times = []

    # Start the main loop
    total_time = 0.0
    with (
        spooler.data_source(source, resource, first, full=False) as spool,
    ):
        for n, mccd in enumerate(spool):
            if server_or_local:
                # Handle the waiting game ...
                # Includes some awkward stuff involving updating on a cycle
                # faster than twait to make the plots more responsive,
                # if twait is long.
                updaters = []
                if mpl:
                    updaters.append(imanager.refresh)

                give_up, try_again, total_time = spooler.hang_about(
                    mccd, twait, tmax, total_time, updaters, tupdate
                )

                if give_up:
                    print("ftargets stopped")
                    break
                elif try_again:
                    continue

            # We have a frame!
            nframe = mccd.head.get("NFRAME", n + 1)
            if source != "hf" and last and nframe > last:
                # finite last frame number
                print("\nHave reduced up to the last frame set.")
                break

            # Now print progress indicators
            tstamp = Time(mccd.head["TIMSTAMP"], format="isot", precision=3)
            tok = "ok" if mccd.head.get("GOODTIME", True) else "nok"
            print(f"{nframe}, utc= {tstamp.iso} ({tok})")

            # We're ready to start plotting
            if mpl:
                imanager.update(mccd)

            if pgp:
                update_implot(
                    mccd,
                    imdevice,
                    ccds,
                    nx,
                    iset,
                    plo,
                    phi,
                    ilo,
                    ihi,
                    xlo,
                    xhi,
                    ylo,
                    yhi,
                )

            # timing
            loop_times.append(time.time())

    # store image details
    total_pix = (xhi - xlo + 1) * (yhi - ylo + 1) * len(ccds)
    binfac = mccd.get_num(0).get_num(0).xbin  # assume same binning for all CCDs

    return start_time, loop_times, total_pix, binfac


if __name__ == "__main__":

    # Test without either
    start_time_none, loop_times_none, total_pix, binfac = main("none")

    # Test PGPLOT
    start_time_pgp, loop_times_pgp, total_pix, binfac = main("pgp")

    # Test MPL
    start_time_mpl, loop_times_mpl, total_pix, binfac = main("mpl")

    # Print
    nframes = len(loop_times_none)  # should be the same for all modes
    print(f"\nFinished processing {nframes} frames")
    print(f"  Binning factor: {binfac}")
    print(f"  Binned pixels per frame: {total_pix:.2e}\n")

    for mode, start_time, loop_times in [
        ("none", start_time_none, loop_times_none),
        ("pgp", start_time_pgp, loop_times_pgp),
        ("mpl", start_time_mpl, loop_times_mpl),
    ]:
        print(f"mode: {mode}")

        total_time = loop_times[-1] - start_time
        print(f"  Total time:                  {total_time:.2f} seconds")
        print(f"    Frames per second:         {nframes / total_time:.2f}")

        # Look at only the second half of the frames, to avoid any startup effects
        second_half_time = loop_times[-1] - loop_times[len(loop_times) // 2]
        print(f"  Total time (second half):    {second_half_time:.2f} seconds")
        print(f"    Frames per second:         {nframes // 2 / second_half_time:.2f}")

        # Find the average time per loop in the second half
        loop_difs = np.diff(loop_times)[len(loop_times) // 2:]
        print(f"  Average loop time:           {np.mean(loop_difs):.3f}±{np.std(loop_difs):.3f} s")

    print("\ncomparison:")
    print("  PGPLOT is x{:.1f} times slower than no plotting".format(
        (loop_times_pgp[-1] - loop_times_pgp[len(loop_times_pgp) // 2])
        / (loop_times_none[-1] - loop_times_none[len(loop_times_none) // 2])
    ))
    print("     MPL is x{:.1f} times slower than PGPLOT".format(
        (loop_times_mpl[-1] - loop_times_mpl[len(loop_times_mpl) // 2])
        / (loop_times_pgp[-1] - loop_times_pgp[len(loop_times_pgp) // 2])
    ))
