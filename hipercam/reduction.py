"""
This contains code used in common by the scripts reduce.py and psf_reduce.py
"""
from collections import OrderedDict
import sys
import warnings
import numpy as np
from astropy.time import Time
import hipercam as hcam
from hipercam import utils, fitting

from trm.pgplot import (
    pgpap,
    pgsubp,
    pgsci,
    pgsch,
    pgenv,
    pglab,
    pgqvp,
    pgvstd,
    pgeras,
    pgerry,
    pgpt,
    pgpanl,
    pgbbuf,
    pgebuf,
    pgsvp,
    pgswin,
    pgbox,
    pgmove,
    pgdraw,
    pgpt1,
)

import matplotlib.pyplot as plt
from hipercam import mpl

NaN = float('NaN')

class Rfile(OrderedDict):
    """Class to read and interpret reduce setup files. Similar to
    configparser but a bit freer. Basically reads lots of lines like
    'ncpu = 1' with sections defined by lines like '[general]'.

    The script genred writes out such files.

    """

    @classmethod
    def read(cls, filename, readcal=True):
        """Builds an Rfile from a reduce file. A few checks for viability are
        applied.

        readcal -- can be used to turn off default reading of
                   calibration data (e.g. see genred)

        """

        rfile = cls()
        insection = False
        with open(filename) as fp:

            for line in fp:

                if not line.startswith("#") and line != "" and not line.isspace():

                    # strip trailing comments
                    comm = line.find("#")
                    if comm != -1:
                        line = line[:comm].strip()

                    if line.startswith("["):
                        # new section
                        section = line[1 : line.find("]")].strip()
                        sec = rfile[section] = OrderedDict()
                        insection = True

                    elif insection:
                        # another entry to current section
                        key = line[: line.find("=")].strip()
                        val = line[line.find("=") + 1 :].strip()

                        if key in sec:
                            if isinstance(sec[key], list):
                                sec[key].append(val)
                            else:
                                sec[key] = [sec[key], val]
                        else:
                            sec[key] = val

                    else:
                        raise hcam.HipercamError(
                            "found entry line before"
                            " any section = \n{:s}".format(line)
                        )

        # process it. this is a matter of checking entries and
        # in some cases converting them to correct or more
        # convenient forms

        #
        # general section
        #
        sect = rfile["general"]
        if sect["version"] != hcam.REDUCE_FILE_VERSION:
            # check the version
            raise hcam.HipercamError(
                f"Version mismatch: reduce file = {sect['version']}, reduce = {hcam.REDUCE_FILE_VERSION}"
                "If your reduce file date is earlier than the reduce version, you may want to try"
                f"updating it using the pipeline command 'rupdate {filename}'"
            )

        sect["scale"] = float(sect["scale"])
        if sect["scale"] <= 0:
            raise hcam.HipercamError("general.scale must be > 0")

        sect["lwidth"] = float(sect["lwidth"])
        if sect["lwidth"] < 0:
            raise hcam.HipercamError("general.lwidth must be >= 0")

        sect["lheight"] = float(sect["lheight"])
        if sect["lheight"] < 0:
            raise hcam.HipercamError("general.lheight must be >= 0")

        sect["iwidth"] = float(sect["iwidth"])
        if sect["iwidth"] < 0:
            raise hcam.HipercamError("general.iwidth must be >= 0")

        sect["iheight"] = float(sect["iheight"])
        if sect["iheight"] < 0:
            raise hcam.HipercamError("general.iheight must be >= 0")

        sect["toffset"] = int(sect["toffset"])

        toBool(rfile, "general", "skipbadt")

        # handle the count level warnings
        warns = sect.get("warn", [])
        if isinstance(warns, str):
            warns = [warns]
        sect['warns'] = warns

        # store a dictionary of warning levels. If
        # any CCD does not have any warning levels,
        # a single warning of this will be issued at
        # the start.
        rfile.warn = {}
        for warn in warns:
            cnam, nonlinear, saturation = warn.split()
            rfile.warn[cnam] = {
                "nonlinear": float(nonlinear),
                "saturation": float(saturation),
            }

        sect["ncpu"] = int(sect["ncpu"])
        if sect["ncpu"] < 1:
            raise hcam.HipercamError("general.ncpu must be >= 1")

        sect["ngroup"] = int(sect["ngroup"]) if sect["ncpu"] > 1 else 1
        if sect["ngroup"] < 1:
            raise hcam.HipercamError("general.ngroup must be >= 1")

        #
        # apertures section
        #
        apsec = rfile["apertures"]

        rfile.aper = hcam.MccdAper.read(
            utils.add_extension(apsec["aperfile"], hcam.APER)
        )
        if apsec["location"] != "fixed" and apsec["location"] != "variable":
            raise hcam.HipercamError(
                "aperture location must either be 'fixed' or 'variable'"
            )

        if apsec["fit_method"] == "moffat":
            rfile.method = "m"
        elif apsec["fit_method"] == "gaussian":
            rfile.method = "g"
        else:
            raise hcam.HipercamError(
                "apertures.fit_method = {:s} not recognised".format(apsec["fit_method"])
            )

        # type conversions
        toBool(rfile, "apertures", "fit_fwhm_fixed")
        toBool(rfile, "apertures", "search_smooth_fft")

        apsec["search_half_width"] = float(apsec["search_half_width"])
        apsec["search_smooth_fwhm"] = float(apsec["search_smooth_fwhm"])
        apsec["fit_fwhm"] = float(apsec["fit_fwhm"])
        apsec["fit_fwhm_min"] = float(apsec["fit_fwhm_min"])
        apsec["fit_ndiv"] = int(apsec["fit_ndiv"])
        apsec["fit_beta"] = float(apsec["fit_beta"])
        apsec["fit_beta_max"] = float(apsec["fit_beta_max"])
        apsec["fit_half_width"] = float(apsec["fit_half_width"])
        apsec["fit_thresh"] = float(apsec["fit_thresh"])
        apsec["fit_height_min_ref"] = float(apsec["fit_height_min_ref"])
        apsec["fit_height_min_nrf"] = float(apsec["fit_height_min_nrf"])
        apsec["fit_max_shift"] = float(apsec["fit_max_shift"])
        apsec["fit_alpha"] = float(apsec["fit_alpha"])
        apsec["fit_diff"] = float(apsec["fit_diff"])

        # run a few checks
        if apsec["fit_beta"] <= 0.0:
            raise hcam.HipercamError("apertures.fit_beta <= 0")
        if apsec["fit_thresh"] <= 1.0:
            raise hcam.HipercamError("apertures.fit_thresh <= 1")
        if apsec["fit_max_shift"] <= 0.0:
            raise hcam.HipercamError("apertures.fit_diff <= 0")
        if apsec["fit_alpha"] <= 0.0 or apsec["fit_alpha"] > 1:
            raise hcam.HipercamError("apertures.fit_alpha outside interval (0,1].")
        if apsec["fit_diff"] <= 0.0:
            raise hcam.HipercamError("apertures.fit_diff <= 0")

        #
        # calibration section
        #
        calsec = rfile["calibration"]

        toBool(rfile, "calibration", "crop")

        if readcal:
            if calsec["bias"] != "":
                rfile.bias = hcam.MCCD.read(
                    utils.add_extension(calsec["bias"], hcam.HCAM)
                )
            else:
                rfile.bias = None

            if calsec["dark"] != "":
                rfile.dark = hcam.MCCD.read(
                    utils.add_extension(calsec["dark"], hcam.HCAM)
                )
            else:
                rfile.dark = None

            if calsec["flat"] != "":
                rfile.flat = hcam.MCCD.read(
                    utils.add_extension(calsec["flat"], hcam.HCAM)
                )
            else:
                rfile.flat = None

            if calsec["fmap"] != "":
                rfile.fmap = hcam.MCCD.read(
                    utils.add_extension(calsec["fmap"], hcam.HCAM)
                )
                rfile.fpair = hcam.fringe.MccdFringePair.read(
                    utils.add_extension(calsec["fpair"], hcam.FRNG)
                )
                calsec["nhalf"] = int(calsec["nhalf"])
                calsec["rmin"] = float(calsec["rmin"])
                calsec["rmax"] = float(calsec["rmax"])
            else:
                rfile.fmap = None

        try:
            rfile.readout = float(calsec["readout"])
        except TypeError:
            if readcal:
                rfile.readout = hcam.MCCD.read(
                    utils.add_extension(calsec["readout"], hcam.HCAM)
                )

        try:
            rfile.gain = float(calsec["gain"])
        except TypeError:
            if readcal:
                rfile.gain = hcam.MCCD.read(
                    utils.add_extension(calsec["gain"], hcam.HCAM)
                )

        # Extraction section

        # Separate the extraction entries into lists, check and
        # convert some entries
        extsec = rfile["extraction"]

        for cnam in extsec:

            extsec[cnam] = lst = extsec[cnam].split()

            if lst[0] != "variable" and lst[0] != "fixed":
                raise hcam.HipercamError(
                    "first entry of extraction lines must either"
                    " be 'variable' or 'fixed'"
                )

            if lst[1] != "normal" and lst[1] != "optimal":
                raise hcam.HipercamError(
                    "second entry of extraction lines"
                    " must either be 'normal' or 'optimal'"
                )

            # type conversions
            for i in range(2, len(lst)):
                lst[i] = float(lst[i])

        #
        # sky section
        #
        skysec = rfile["sky"]

        if skysec["error"] == "variance":
            if skysec["method"] == "median":
                raise hcam.HipercamError(
                    "sky.error == variance requires sky.method == clipped"
                )
        elif skysec["error"] != "photon":
            raise hcam.HipercamError(
                "sky.error must be either 'variance' or 'photon'"
            )

        if skysec["method"] != "clipped" and skysec["method"] != "median":
            raise hcam.HipercamError(
                "sky.method must be either 'clipped' or 'median'"
            )

        skysec["thresh"] = float(skysec["thresh"])

        #
        # light curve plot section
        #

        sect = rfile["lcplot"]

        sect["xrange"] = float(sect["xrange"])
        sect["extend_x"] = float(sect["extend_x"])
        if sect["extend_x"] <= 0:
            raise hcam.HipercamError("lcplot.extend_x must be > 0")

        #
        # light curve panel section
        #

        sect = rfile["light"]

        plot = sect["plot"]
        if isinstance(plot, str):
            plot = [plot]

        # convert entries to the right type here and try colours to
        # PGPLOT colour indices.
        for n in range(len(plot)):
            cnam, tnm, cnm, off, fac, dcol, ecol = plot[n].split()

            # some consistency checks
            if cnam not in rfile.aper:
                raise hcam.HipercamError(
                    (
                        "CCD = {:s} has a plot line in 'light' but is"
                        " not present in the aperture file"
                    ).format(cnam)
                )

            if cnam not in extsec:
                raise hcam.HipercamError(
                    (
                        "CCD = {:s} has a plot line in 'light' but no"
                        " corresponding entry under 'extraction'"
                    ).format(cnam)
                )

            if tnm not in rfile.aper[cnam]:
                raise hcam.HipercamError(
                    (
                        "Aperture = {:s} [target], CCD = {:s} has a plot line"
                        " in 'light' but is not present in the aperture file"
                    ).format(tnm, cnam)
                )

            if cnm != "!" and cnm not in rfile.aper[cnam]:
                raise hcam.HipercamError(
                    (
                        "Aperture = {:s} [comparison], CCD = {:s} has a plot line"
                        " in 'light' but is not present in the aperture file"
                    ).format(tnm, cnam)
                )

            plot[n] = {
                "ccd": cnam,
                "targ": tnm,
                "comp": cnm,
                "off": float(off),
                "fac": float(fac),
                "dcol": ctrans(dcol),
                "ecol": ctrans(ecol),
            }
        sect["plot"] = plot

        toBool(rfile, "light", "linear")
        toBool(rfile, "light", "y_fixed")
        sect["y1"] = float(sect["y1"])
        sect["y2"] = float(sect["y2"])
        sect["extend_y"] = float(sect["extend_y"])
        if sect["extend_y"] <= 0:
            raise hcam.HipercamError("light.extend_y must be > 0")

        #
        # position panel section
        #

        rfile.position = "position" in rfile
        if rfile.position:
            sect = rfile["position"]

            plot = sect["plot"]
            if isinstance(plot, str):
                plot = [plot]

            # convert entries to the right type here and try colours to
            # PGPLOT colour indices.
            for n in range(len(plot)):
                cnam, tnm, dcol, ecol = plot[n].split()

                # some consistency checks
                if cnam not in rfile.aper:
                    raise hcam.HipercamError(
                        (
                            "CCD = {:s} has a plot line in 'position'"
                            " but is not present in the aperture file"
                        ).format(cnam)
                    )

                if cnam not in extsec:
                    raise hcam.HipercamError(
                        (
                            "CCD = {:s} has a plot line in 'position'"
                            " but is no corresponding entry under 'extraction'"
                        ).format(cnam)
                    )

                if tnm not in rfile.aper[cnam]:
                    raise hcam.HipercamError(
                        (
                            "Aperture = {:s}, CCD = {:s} has a plot line in 'position'"
                            " but is not present in the aperture file"
                        ).format(tnm, cnam)
                    )

                plot[n] = {
                    "ccd": cnam,
                    "targ": tnm,
                    "dcol": ctrans(dcol),
                    "ecol": ctrans(ecol),
                }
            sect["plot"] = plot

            sect["height"] = float(sect["height"])
            if sect["height"] <= 0:
                raise hcam.HipercamError("position.height must be > 0")

            toBool(rfile, "position", "x_fixed")
            sect["x_min"] = float(sect["x_min"])
            sect["x_max"] = float(sect["x_max"])
            if sect["x_max"] <= sect["x_min"]:
                raise hcam.HipercamError("position.x_min must be < position.x_max")

            toBool(rfile, "position", "y_fixed")
            sect["y_min"] = float(sect["y_min"])
            sect["y_max"] = float(sect["y_max"])
            if sect["y_max"] <= sect["y_min"]:
                raise hcam.HipercamError("position.y_min must be < position.y_max")

            sect["extend_y"] = float(sect["extend_y"])
            if sect["extend_y"] <= 0:
                raise hcam.HipercamError("position.extend_y must be > 0")

        #
        # transmission panel section
        #

        rfile.transmission = "transmission" in rfile
        if rfile.transmission:
            sect = rfile["transmission"]

            plot = sect["plot"]
            if isinstance(plot, str):
                plot = [plot]

            # convert entries to the right type here and try colours to
            # PGPLOT colour indices.
            for n in range(len(plot)):
                cnam, tnm, dcol, ecol = plot[n].split()

                # some consistency checks
                if cnam not in rfile.aper:
                    raise hcam.HipercamError(
                        (
                            "CCD = {:s} has a plot line in 'transmission' but is"
                            " not present in the aperture file"
                        ).format(cnam)
                    )

                if cnam not in extsec:
                    raise hcam.HipercamError(
                        (
                            "CCD = {:s} has a plot line in 'transmission' but no"
                            " corresponding entry under 'extraction'"
                        ).format(cnam)
                    )

                if tnm not in rfile.aper[cnam]:
                    raise hcam.HipercamError(
                        (
                            "Aperture = {:s}, CCD = {:s} has a plot line in 'transmission'"
                            " but is not present in the aperture file"
                        ).format(tnm, cnam)
                    )

                plot[n] = {
                    "ccd": cnam,
                    "targ": tnm,
                    "dcol": ctrans(dcol),
                    "ecol": ctrans(ecol),
                }
            sect["plot"] = plot

            sect["height"] = float(sect["height"])
            if sect["height"] <= 0:
                raise hcam.HipercamError("transmission.height must be > 0")

            sect["ymax"] = float(sect["ymax"])
            if sect["ymax"] < 100:
                raise hcam.HipercamError("transmission.ymax must be >= 100")

        #
        # seeing panel section
        #
        rfile.seeing = "seeing" in rfile
        if rfile.seeing:
            sect = rfile["seeing"]

            plot = sect["plot"]
            if isinstance(plot, str):
                plot = [plot]

            # convert entries to the right type here and try colours to
            # PGPLOT colour indices.
            for n in range(len(plot)):
                cnam, tnm, dcol, ecol = plot[n].split()

                # some consistency checks
                if cnam not in rfile.aper:
                    raise hcam.HipercamError(
                        (
                            "CCD = {:s} has a plot line in 'position' but is"
                            " not present in the aperture file"
                        ).format(cnam)
                    )

                if cnam not in extsec:
                    raise hcam.HipercamError(
                        (
                            "CCD = {:s} has a plot line in 'position' but no"
                            " corresponding entry under 'extraction'"
                        ).format(cnam)
                    )

                if tnm not in rfile.aper[cnam]:
                    raise hcam.HipercamError(
                        (
                            "Aperture = {:s}, CCD = {:s} has a plot line in 'position'"
                            " but is not present in the aperture file"
                        ).format(tnm, cnam)
                    )

                plot[n] = {
                    "ccd": cnam,
                    "targ": tnm,
                    "dcol": ctrans(dcol),
                    "ecol": ctrans(ecol),
                }
            sect["plot"] = plot

            toBool(rfile, "seeing", "y_fixed")

            sect["height"] = float(sect["height"])
            if sect["height"] <= 0:
                raise hcam.HipercamError("seeing.height must be > 0")

            sect["ymax"] = float(sect["ymax"])
            if sect["ymax"] <= 0:
                raise hcam.HipercamError("seeing.ymax must be > 0")

            sect["extend_y"] = float(sect["extend_y"])
            if sect["extend_y"] <= 0:
                raise hcam.HipercamError("seeing.extend_y must be > 0")

        # focal plane mask section
        toBool(rfile, "focal_mask", "demask")
        sect = rfile["focal_mask"]
        sect["dthresh"] = float(sect["dthresh"])

        # Monitor section
        monsec = rfile["monitor"]
        for apnam in monsec:

            # interpret the bitmasks.
            monsec[apnam] = [eval("hcam." + entry) for entry in monsec[apnam].split()]

        # We are finally done reading and checking the reduce script.
        # rfile[section][param] should from now on return something
        # useful and somewhat reliable.
        rfile.filename = filename
        return rfile

    def crop(self, mccd):
        """This uses a template file 'mccd' to try to get the calibration
        files into the same format.

        """
        if self.bias is not None:
            self.bias = self.bias.crop(mccd)

        if self.dark is not None:
            self.dark = self.dark.crop(mccd)

        if self.flat is not None:
            self.flat = self.flat.crop(mccd)

        if self.fmap is not None:
            self.fmap = self.fmap.crop(mccd)
            nhalf = self['calibration']['nhalf']
            self.fpair = self.fpair.crop(mccd, nhalf)

        if isinstance(self.readout, hcam.MCCD):
            self.readout = self.readout.crop(mccd)

        if isinstance(self.gain, hcam.MCCD):
            self.gain = self.gain.crop(mccd)


def setup_plots(rfile, ccds, nx, plot_lims, implot=True, lplot=True):
    """
    Perform initial setup of image and lightcurve plots for reduction outputs
    """
    imdev = lcdev = spanel = tpanel = xpanel = ypanel = lpanel = None
    if implot:
        xlo, xhi, ylo, yhi = plot_lims
        # image plot
        imdev = hcam.pgp.Device(rfile["general"]["idevice"])
        iwidth = rfile["general"]["iwidth"]
        iheight = rfile["general"]["iheight"]

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

    if lplot:
        # open the light curve plot
        lcdev = hcam.pgp.Device(rfile["general"]["ldevice"])
        lwidth = rfile["general"]["lwidth"]
        lheight = rfile["general"]["lheight"]
        if lwidth > 0 and lheight > 0:
            pgpap(lwidth, lheight / lwidth)

        # define and draw panels. 'total' is the total height of all
        # panel which will be used to work out the fraction occupied
        # by each one. There is always a light curve panel of unit
        # height.
        lheight = 1.0
        total = lheight

        if rfile.position:
            pheight = rfile["position"]["height"]
            total += pheight

        if rfile.transmission:
            theight = rfile["transmission"]["height"]
            total += theight

        if rfile.seeing:
            sheight = rfile["seeing"]["height"]
            total += sheight

        # define standard viewport. We set the character height to ensure
        # there is enough room around the edges.
        pgsch(max(hcam.pgp.Params["axis.label.ch"], hcam.pgp.Params["axis.number.ch"]))
        pgvstd()
        xv1, xv2, yv1, yv2 = pgqvp()
        x1, x2 = 0, rfile["lcplot"]["extend_x"]

        # scale = vertical height of LC panel in device coordinates
        scale = (yv2 - yv1) / total

        # Work from the bottom up to get the X-axis labelling right
        # since it has to be turned off for the upper panels.  at each
        # step yv1, yv2 is the Panel's vertical range
        xlabel = "Time [mins]"
        xopt = "bcnst"

        if rfile.seeing:
            # the seeing panel
            yv2 = yv1 + scale * sheight
            spanel = Panel(
                lcdev,
                xv1,
                xv2,
                yv1,
                yv2,
                xlabel,
                'FWHM (")',
                "",
                xopt,
                "bcnst",
                x1,
                x2,
                0,
                rfile["seeing"]["ymax"],
            )

            xlabel = ""
            xopt = "bcst"
            yv1 = yv2

        if rfile.transmission:
            # the transmission panel
            yv2 = yv1 + scale * theight
            tpanel = Panel(
                lcdev,
                xv1,
                xv2,
                yv1,
                yv2,
                xlabel,
                "% trans",
                "",
                xopt,
                "bcnst",
                x1,
                x2,
                0,
                rfile["transmission"]["ymax"],
            )

            xlabel = ""
            xopt = "bcst"
            yv1 = yv2

        if rfile.position:
            # the X,Y position panels. First Y
            yv2 = yv1 + scale * pheight / 2.0
            ypanel = Panel(
                lcdev,
                xv1,
                xv2,
                yv1,
                yv2,
                xlabel,
                "Y",
                "",
                xopt,
                "bcnst",
                x1,
                x2,
                rfile["position"]["y_min"],
                rfile["position"]["y_max"],
            )

            xlabel = ""
            xopt = "bcst"
            yv1 = yv2

            # then X
            yv2 = yv1 + scale * pheight / 2.0
            xpanel = Panel(
                lcdev,
                xv1,
                xv2,
                yv1,
                yv2,
                xlabel,
                "X",
                "",
                xopt,
                "bcnst",
                x1,
                x2,
                rfile["position"]["x_min"],
                rfile["position"]["x_max"],
            )

            xlabel = ""
            xopt = "bcst"
            yv1 = yv2

        # Light curve is not an option and is always at the top
        yv2 = yv1 + scale * lheight

        lpanel = Panel(
            lcdev,
            xv1,
            xv2,
            yv1,
            yv2,
            xlabel,
            "Flux" if rfile["light"]["linear"] else "Magnitudes",
            "",
            xopt,
            "bcnst",
            x1,
            x2,
            rfile["light"]["y1"],
            rfile["light"]["y2"],
        )

    return imdev, lcdev, spanel, tpanel, xpanel, ypanel, lpanel


def setup_plot_buffers(rfile):
    # create buffers to store the points plotted in each panel
    lbuffer = []
    for plot_config in rfile["light"]["plot"]:
        lbuffer.append(LightCurve(plot_config, rfile["light"]["linear"]))

    xbuffer, ybuffer = [], []
    if rfile.position:
        for plot_config in rfile["position"]["plot"]:
            xbuffer.append(Xposition(plot_config))
            ybuffer.append(Yposition(plot_config))

    tbuffer = []
    if rfile.transmission:
        for plot_config in rfile["transmission"]["plot"]:
            tbuffer.append(Transmission(plot_config))

    sbuffer = []
    if rfile.seeing:
        for plot_config in rfile["seeing"]["plot"]:
            sbuffer.append(Seeing(plot_config, rfile["general"]["scale"]))

    return lbuffer, xbuffer, ybuffer, tbuffer, sbuffer


# Some fancy stuff to keep tzero internal to update_plots
def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func

    return decorate


@static_vars(tzero=None)
def update_plots(
    results,
    rfile,
    implot,
    lplot,
    imdev,
    lcdev,
    pccd,
    ccds,
    msub,
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
    lpanel,
    xpanel,
    ypanel,
    tpanel,
    spanel,
    tkeep,
    lbuffer,
    xbuffer,
    ybuffer,
    tbuffer,
    sbuffer,
):
    """Plot updater routine.

    Arguments::

       results : list
          contains tuples with the results for groups of frames for
          each CCD. See LogWriter.write for more details on its
          structure.

       rfile : Rfile
          reduce file parameters

       implot : bool
          whether to plot images

       lplot : bool
          whether to plot light curves

       imdev : hipercam.pgp.Device
          represents the image plot device (if any)

       lcdev : hipercam.pgp.Device
          represents the light curve plot device (if any)

       pccd : MCCD [only matters if imdev]
          current debiassed, flat-fielded frame. In the case of frame
          groups, this will be the last one added.

       ccds : list [only matters if imdev]
          list of CCD labels to be plotted

       msub : bool [only matters if imdev]
          True to force subtraction of median value

       nx : int [only matters if imdev]
          Number of panels in X for image display

       iset : character [only matters if imdev]
          'a', 'd', or 'p' to define image display levels

       plo : float [only matters if imdev and iset =='p']
          Lower percentile level

       phi : float [only matters if imdev and iset =='p']
          Upper percentile level

       ilo : float [only matters if imdev and iset =='d']
          Lower intensity level

       ihi : float [only matters if imdev and iset =='d']
          Upper intensity level

       xlo : int
          left X-limit for determination of percentiles

       xhi : int
          right X-limit for determination of percentiles

       ylo : int
          lower Y-limit for determination of percentiles

       yhi : int
          upper Y-limit for determination of percentiles

       lpanel, xpanel, ypanel, tpanel, spanel, tkeep, lbuffer,
        xbuffer, ybuffer, tbuffer, sbuffer -- need doing.

       lbuffer : list of LightCurve

    """

    # now the plotting sections
    if implot:
        # image plot
        # select the image device
        imdev.select()

        # display the CCDs chosen
        message = "; "
        for nc, cnam in enumerate(ccds):
            ccd = pccd[cnam]

            if ccd.is_data():
                # this should be data as opposed to a blank frame
                # between data frames that occur with nskip > 0

                if msub:
                    # subtract median from each window
                    for wind in ccd.values():
                        wind -= wind.median()

                # set to the correct panel and then plot CCD
                ix = (nc % nx) + 1
                iy = nc // nx + 1
                pgpanl(ix, iy)
                vmin, vmax = hcam.pgp.pCcd(
                    ccd, iset, plo, phi, ilo, ihi, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi
                )

                # accumulate string of image scalings
                if nc:
                    message += ", ccd {:s}: {:.2f} to {:.2f}".format(cnam, vmin, vmax)
                else:
                    message += "ccd {:s}: {:.2f} to {:.2f}".format(cnam, vmin, vmax)

                if cnam in rfile.aper:
                    # find the result for this ccd name
                    res = next(res for (this_cnam, res) in results if this_cnam == cnam)
                    # we always use the last result in the parallel group
                    last_res = res[-1]
                    ccd_aper = last_res[2]
                    hcam.pgp.pCcdAper(ccd_aper)

        # end of CCD display loop
        print(message)

    if lplot:

        # plot the light curve

        if update_plots.tzero is None:

            # first set the tzero if it is not set yet

            skipbadt = rfile["general"]["skipbadt"]
            for cnam, res in results:
                for (
                    nframe,
                    store,
                    ccdaper,
                    reses,
                    mjdint,
                    mjdfrac,
                    mjdok,
                    expose,
                ) in res:
                    if mjdok or not skipbadt:
                        update_plots.tzero = mjdint + mjdfrac
                        break
                if update_plots.tzero is not None:
                    break

        if update_plots.tzero is None:
            # Must have a tzero to get going at all
            return

        tzero = update_plots.tzero

        # track the maximum time
        tmax = None

        replot, ltmax = plotLight(lpanel, tzero, results, rfile, tkeep, lbuffer)
        tmax = tmax if ltmax is None else ltmax if tmax is None else max(tmax, ltmax)

        if rfile.position:
            # plot the positions
            rep, ptmax = plotPosition(
                xpanel, ypanel, tzero, results, rfile, tkeep, xbuffer, ybuffer
            )
            replot |= rep
            tmax = (
                tmax if ptmax is None else ptmax if tmax is None else max(tmax, ptmax)
            )

        if rfile.transmission:
            # plot the transmission
            rep, ttmax = plotTrans(tpanel, tzero, results, rfile, tkeep, tbuffer)
            replot |= rep
            tmax = (
                tmax if ttmax is None else ttmax if tmax is None else max(tmax, ttmax)
            )

        if rfile.seeing:
            # plot the seeing
            rep, stmax = plotSeeing(spanel, tzero, results, rfile, tkeep, sbuffer)
            replot |= rep
            tmax = (
                tmax if stmax is None else stmax if tmax is None else max(tmax, stmax)
            )

        # check the time
        if tmax is not None and tmax > lpanel.x2:
            # extend range by multiple of extend_x
            x2 = lpanel.x2
            while tmax > x2:
                x2 += rfile["lcplot"]["extend_x"]

            # need to re-plot
            replot = True
            lpanel.x2 = x2

            if rfile.position:
                xpanel.x2 = x2
                ypanel.x2 = x2

            if rfile.transmission:
                tpanel.x2 = x2

            if rfile.seeing:
                spanel.x2 = x2

        if replot:

            # re-plot.

            if tkeep > 0:
                # adjust start time if tkeep in use
                tstart = max(0, lpanel.x2 - tkeep)

                lpanel.x1 = tstart

                if rfile.position:
                    xpanel.x1 = tstart
                    ypanel.x1 = tstart

                if rfile.transmission:
                    tpanel.x1 = tstart

                if rfile.seeing:
                    spanel.x1 = tstart

            # start buffering
            pgbbuf()

            # erase plot
            pgeras()

            # re-draw the light curve panel
            lpanel.plot()

            for lc in lbuffer:
                # convert the buffered data into float32 ndarrays
                t = np.array(lc.t, dtype=np.float32)
                f = np.array(lc.f, dtype=np.float32)
                fe = np.array(lc.fe, dtype=np.float32)
                symbs = np.array(lc.symb, dtype=np.int)
                asymbs = set(symbs)

                # Plot the error bars
                if lc.ecol is not None:
                    pgsci(lc.ecol)
                    pgerry(t, f - fe, f + fe, 0)

                # Plot the data
                pgsci(lc.dcol)
                pgsch(0.5)
                for symb in asymbs:
                    ok = symbs == symb
                    pgpt(t[ok], f[ok], symb)

            if rfile.position:
                # re-draw the position panels

                xpanel.plot()
                for xpos in xbuffer:
                    # convert the buffered data into float32 ndarrays
                    t = np.array(xpos.t, dtype=np.float32)
                    f = np.array(xpos.f, dtype=np.float32)
                    fe = np.array(xpos.fe, dtype=np.float32)
                    symbs = np.array(xpos.symb, dtype=np.int)
                    asymbs = set(symbs)

                    # Plot the error bars
                    if xpos.ecol is not None:
                        pgsci(xpos.ecol)
                        pgerry(t, f - fe, f + fe, 0)

                    # Plot the data
                    pgsci(xpos.dcol)
                    pgsch(0.5)
                    for symb in asymbs:
                        ok = symbs == symb
                        pgpt(t[ok], f[ok], symb)

                ypanel.plot()
                for ypos in ybuffer:
                    # convert the buffered data into float32 ndarrays
                    t = np.array(ypos.t, dtype=np.float32)
                    f = np.array(ypos.f, dtype=np.float32)
                    fe = np.array(ypos.fe, dtype=np.float32)
                    symbs = np.array(ypos.symb, dtype=np.int)
                    asymbs = set(symbs)

                    # Plot the error bars
                    if ypos.ecol is not None:
                        pgsci(ypos.ecol)
                        pgerry(t, f - fe, f + fe, 0)

                    # Plot the data
                    pgsci(ypos.dcol)
                    pgsch(0.5)
                    for symb in asymbs:
                        ok = symbs == symb
                        pgpt(t[ok], f[ok], symb)

            if rfile.transmission:
                # re-draw the transmission panel
                tpanel.plot()

                for trans in tbuffer:
                    if trans.fmax:
                        # convert the buffered data into float32
                        # ndarrays
                        t = np.array(trans.t, dtype=np.float32)
                        f = np.array(trans.f, dtype=np.float32)
                        fe = np.array(trans.fe, dtype=np.float32)
                        scale = np.float32(100.0 / trans.fmax)
                        f *= scale
                        fe *= scale
                        symbs = np.array(trans.symb, dtype=np.int)
                        asymbs = set(symbs)

                        # Plot the error bars
                        if trans.ecol is not None:
                            pgsci(trans.ecol)
                            pgerry(t, f - fe, f + fe, 0)

                        # Plot the data
                        pgsci(trans.dcol)
                        pgsch(0.5)
                        for symb in asymbs:
                            ok = symbs == symb
                            pgpt(t[ok], f[ok], symb)

            if rfile.seeing:
                # re-draw the seeing panel
                spanel.plot()

                for see in sbuffer:
                    # convert the buffered data into float32 ndarrays
                    t = np.array(see.t, dtype=np.float32)
                    f = np.array(see.f, dtype=np.float32)
                    fe = np.array(see.fe, dtype=np.float32)
                    symbs = np.array(see.symb, dtype=np.int)
                    asymbs = set(symbs)

                    # Plot the error bars
                    if see.ecol is not None:
                        pgsci(see.ecol)
                        pgerry(t, f - fe, f + fe, 0)

                    # Plot the data
                    pgsci(see.dcol)
                    pgsch(0.5)
                    for symb in asymbs:
                        ok = symbs == symb
                        pgpt(t[ok], f[ok], symb)

            # end buffering
            pgebuf()


def initial_checks(mccd, rfile):
    """
    Perform sanity checks of aperture settings, and process calibrations
    """
    # This is the very first OK data we have. There are a few
    # things we do now on the assumption that all the data frames
    # will have the same format. Begin by checking for some poor
    # settings
    xbin = ybin = 1
    ok = True

    for cnam in mccd:
        for wnam in mccd[cnam]:
            xbin = max(xbin, mccd[cnam][wnam].xbin)
            ybin = max(ybin, mccd[cnam][wnam].ybin)

    if rfile["apertures"]["fit_half_width"] < 3 * max(xbin, ybin):
        print(
            (
                "'** fit_half_width' < 3*max(xbin,ybin) ({:.1f} vs {:d}); profile fits"
                " will fail. Please increase its value in the reduce file."
            ).format(rfile["apertures"]["fit_half_width"], 3 * max(xbin, ybin)),
            file=sys.stderr,
        )
        ok = False

    elif rfile["apertures"]["fit_half_width"] < 5 * max(xbin, ybin):
        warnings.warn(
            (
                "'** fit_half_width' < 5*max(xbin,ybin) ({:.1f} vs {:d}) ==> small windows for"
                " profile fits. Advise increasing its value in the reduce file if possible."
            ).format(rfile["apertures"]["fit_half_width"], 5 * max(xbin, ybin))
        )

    if rfile["calibration"]["crop"]:
        # Trim the calibrations, on the assumption that all
        # data frames have the same format
        rfile.crop(mccd)

    # create a readout noise frame if none read in
    if isinstance(rfile.readout, hcam.MCCD):
        read = rfile.readout
    else:
        read = mccd.copy()
        read.set_const(rfile.readout)

    # create a gain frame if none read in
    if isinstance(rfile.gain, hcam.MCCD):
        gain = rfile.gain
    else:
        gain = mccd.copy()
        gain.set_const(rfile.gain)

    if rfile.flat is not None:
        # correct these frames by the flat, if there is one.
        # These corrections allow the use of the usual
        # read**2+data/gain as the variance, after the data
        # have been flat-fielded.
        read /= rfile.flat
        gain *= rfile.flat

    return (read, gain, ok)


class ProcessCCDs:
    """Function object replacing old "process_ccds" function. It
    handles some of the persistent objects & arguemts required by
    that routine more gracefully.
    """

    def __init__(self, rfile, read, gain, ccdproc, pool):
        """Arguments::

           rfile : Rfile
              reduce file settings

           read : MCCD
              CCD representing readout noise divided by the flat field

           gain : MCCD
              CCD representing the gain multiplied by the flat field

           ccdproc : function
              function that carries out the per-CCD processing needed for the
              parallelisation step.  See 'scripts/reduce.py' for its calling
              signature.

           pool : multiprocessing.Pool | None
              for parallel processing. If 'None', it will be done in serial
              mode.

        """

        self.rfile = rfile
        self.read = read
        self.gain = gain
        self.ccdproc = ccdproc
        self.pool = pool

        # mwins: a persistent storage container that is set up in the
        # first call to this object and retained for later access. It
        # is a dictionary of dictionaries such that mwin[cnam][apnam]
        # = the window name containing the aperture. It is assumed
        # that the apertures do not drift from one window to another,
        # so this is only created once.
        self.mwins = {}

        # store: another container to hold some of the results that need
        # propagating between one group of frames and the next
        self.store = {}

    def __call__(self, pccds, bccds, mccds, nframes):

        """Carries out the multi-processing reduction of a set of frames

        Arguments::

           pccds : list of MCCDs
              processed frames (debiassed, flat-fielded). The idea is
              to read in a few frames before sending them off to
              parallel processing in order to reduce the overheads
              that are incurred by processing each frame as it comes
              in (old method). All frames should have an identical set
              of CCDs in them.

           bccds : list of MCCDs
              purely debiassed versions, used for variance computation

           mccds : list of MCCDs
              their unprocessed counterparts, used to determine saturation.

        Returns with list of ccdproc results, one set for each CCD.
        See 'ccdproc' in reduce.py for definition of the return which
        should be the same for any version.

        """

        arglist, allres = [], []

        for cnam in pccds[0]:

            # get the apertures
            if (
                cnam not in self.rfile.aper
                or cnam not in self.rfile["extraction"]
                or len(self.rfile.aper[cnam]) == 0
            ):
                continue

            # build lists of specific CCDs from each frame as opposed to lists
            # of frames. Only accumulate ones marked as having data (i.e. skip
            # NSKIP frames).  this is where the parellisation split is made.
            pcds, bcds, mcds, nfrs = [], [], [], []
            for pccd, bccd, mccd, nframe in zip(pccds, bccds, mccds, nframes):
                if pccd[cnam].is_data():
                    pcds.append(pccd[cnam])
                    bcds.append(bccd[cnam])
                    mcds.append(mccd[cnam])
                    nfrs.append(nframe)

            if len(pcds) == 0:
                continue

            if cnam not in self.mwins:
                # first time through, work out an array of which window each
                # aperture lies in. we will assume this is fixed for the whole
                # run, i.e. that apertures do not drift from one window to
                # another. Set to None if no window found
                self.mwins[cnam] = {}
                for apnam, aper in self.rfile.aper[cnam].items():
                    for wnam, wind in mccd[cnam].items():
                        if wind.distance(aper.x, aper.y) > 0:
                            self.mwins[cnam][apnam] = wnam
                            break
                        else:
                            self.mwins[cnam][apnam] = None

            if cnam not in self.store:
                # initialisation
                self.store[cnam] = {"mfwhm": -1., "mbeta": -1.}

            if self.pool is None:
                # carry out processing serially, store results
                res = self.ccdproc(
                    cnam,
                    pcds,
                    bcds,
                    mcds,
                    nfrs,
                    self.read[cnam],
                    self.gain[cnam],
                    self.mwins[cnam],
                    self.rfile,
                    self.store[cnam],
                )
                allres.append(res)

            else:
                # store arguments lists for later parallel processing
                arglist.append(
                    (
                        cnam,
                        pcds,
                        bcds,
                        mcds,
                        nfrs,
                        self.read[cnam],
                        self.gain[cnam],
                        self.mwins[cnam],
                        self.rfile,
                        self.store[cnam],
                    )
                )

        if len(arglist):
            # delayed reduction now carried out in parallel
            allres = self.pool.starmap(self.ccdproc, arglist)

            # update store and rfile.aper with the final results for
            # each CCD. This is necessary because the alterations to
            # these made within the routine are otherwise lost
            for cnam, res in allres:
                if len(res):
                    (nframe, store, ccdaper, results, mjdint, mjdfrac, mjdok, expose) = res[-1]
                    self.store[cnam] = store
                    self.rfile.aper[cnam] = ccdaper

        # pass back the results
        return allres

def moveApers(cnam, ccd, bccd, read, gain, ccdwin, rfile, store):
    """Encapsulates aperture re-positioning. 'store' is a dictionary of results
    that will be used to start the fits from one frame to the next. It must
    start as {'mfwhm': -1., 'mbeta': -1.}. The values of these will be revised
    and extra stuff added in addition.

    It operates by first shifting any reference apertures, then non-linked
    apertures, and finally linked apertures. The overall sequence is quite
    complex and is described elsewhere.

    Arguments::

       cnam : string
           CCD label

       ccd : CCD
           the processed CCD.

       bccd : CCD
           debiassed only CCD.

       read : CCD
           readnoise

       gain : CCD
           gain, e-/ADU.

       ccdwin : dict
           the Window label corresponding to each Aperture

       rfile : Rfile
           reduce file configuration parameters

       store : dict
           used to store various parameters needed further down the
           line.  Initialise to {'mfwhm': -1., 'mbeta': -1.}. For each
           aperture in ccdaper, parameters `xe`, `ye`, `fwhm`,
           `fwhme`, `beta`, `betae`, `dx`, `dy` will be added. They
           represent: the uncertainty in the fitted X, Y position (ex,
           ey), the fitted FWHM and its uncertainty (fwhm, fwhme), the
           same for beta and its uncertainty (beta, betae), and
           finally the change in x and y position.

    Returns: True or False to indicate whether to move onto extraction or not.
    If False, extraction will be skipped.

    """

    ccdaper = rfile.aper[cnam]

    # short-hand that will used a lot
    apsec = rfile["apertures"]

    if apsec["location"] == "fixed":
        for apnam in ccdaper:
            store[apnam] = {
                "xe": NaN,
                "ye": NaN,
                "fwhm": NaN,
                "fwhme": NaN,
                "beta": NaN,
                "betae": NaN,
                "dx": 0.,
                "dy": 0.,
            }
        return True

    # some initialisations
    xsum, ysum = 0.0, 0.0
    wxsum, wysum = 0.0, 0.0
    ref = False
    fsum, wfsum = 0, 0
    bsum, wbsum = 0, 0
    shbox = apsec["search_half_width"]

    # first of all try to get a mean shift from the reference apertures.  we
    # move any of these apertures that are fitted OK. We work out weighted
    # mean FWHM and beta values once any other apertures are fitted.

    # next to store shifts
    shifts = []

    # next is a shift that will be updated after the first fit. This allows
    # the later reference stars to be not so good as the first.
    xshift, yshift = 0., 0.

    try:
        for apnam, aper in ccdaper.items():
            if aper.ref:
                ref = True

                # name of window for this aperture
                wnam = ccdwin[apnam]

                # extract the Window of the processed data, read noise and
                # gain frames
                wdata = ccd[wnam]
                wbias = bccd[wnam]
                wread = read[wnam]
                wgain = gain[wnam]

                # get sub-window around start position
                swdata = wdata.window(
                    aper.x + xshift - shbox, aper.x + xshift + shbox,
                    aper.y + yshift - shbox, aper.y + yshift + shbox
                )

                # carry out initial search
                x, y, peak = swdata.search(
                    apsec["search_smooth_fwhm"],
                    aper.x + xshift,
                    aper.y + yshift,
                    apsec["fit_height_min_ref"],
                    apsec["search_smooth_fft"],
                )

                # Now for a more refined fit. First extract fit Window
                fhbox = apsec["fit_half_width"]
                fwdata = wdata.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)
                fwbias = wbias.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)
                fwread = wread.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)
                fwgain = wgain.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)

                # initial estimate of background
                sky = np.percentile(fwdata.data, 50)

                # get some parameters from the previous run where possible
                fit_fwhm = store["mfwhm"] if store["mfwhm"] > 0.0 else apsec["fit_fwhm"]
                fit_beta = store["mbeta"] if store["mbeta"] > 0.0 else apsec["fit_beta"]

                # limit the initial value of beta because of tendency to
                # wander to high values and never come down.
                fit_beta = min(fit_beta, apsec["fit_beta_max"])

                # This is where the pure-debiassed data are used
                sigma = np.sqrt(fwread.data**2+np.maximum(0,fwbias.data)/fwgain.data)

                # refine the Aperture position by fitting the profile
                (
                    (sky, height, x, y, fwhm, beta),
                    (esky, eheight, ex, ey, efwhm, ebeta),
                    extras,
                ) = hcam.fitting.combFit(
                    fwdata, sigma,
                    rfile.method,
                    sky,
                    peak - sky,
                    x,
                    y,
                    fit_fwhm,
                    apsec["fit_fwhm_min"],
                    apsec["fit_fwhm_fixed"],
                    fit_beta,
                    apsec["fit_beta_max"],
                    False,
                    apsec["fit_thresh"],
                    apsec["fit_ndiv"],
                )

                # check that x, y are within both the search and fit boxes
                if swdata.distance(x, y) < 0.5:
                    raise hcam.HipercamError(
                        (
                            "Fitted position ({:.1f},{:.1f}) too close to or"
                            " beyond edge of search window = {!s}".format(
                                x, y, swdata.format(True)
                            )
                        )
                    )
                elif fwdata.distance(x, y) < 0.5:
                    raise hcam.HipercamError(
                        (
                            "Fitted position ({:.1f},{:.1f}) too close to or"
                            " beyond edge of fit window = {!s}".format(
                                x, y, swdata.format(True)
                            )
                        )
                    )

                if height >= apsec["fit_height_min_ref"]:

                    # The peak height check is probably not required at this
                    # point since the search routine applies it more
                    # stringently but I have left it in for safety.
                    dx = x - aper.x
                    wx = 1.0 / ex ** 2
                    wxsum += wx
                    xsum += wx * dx

                    dy = y - aper.y
                    wy = 1.0 / ey ** 2
                    wysum += wy
                    ysum += wy * dy

                    if len(shifts) == 0:
                        # store the first to help any subsequent ones
                        xshift, yshift = dx, dy

                    shifts.append((apnam, dx, dy))


                    # store stuff
                    store[apnam] = {
                        "xe": ex,
                        "ye": ey,
                        "fwhm": fwhm,
                        "fwhme": efwhm,
                        "beta": beta,
                        "betae": ebeta,
                        "dx": dx,
                        "dy": dy,
                    }

                    if efwhm > 0.0:
                        # average FWHM computation
                        wf = 1.0 / efwhm ** 2
                        fsum += wf * fwhm
                        wfsum += wf

                    if ebeta > 0.0:
                        # average beta computation
                        wb = 1.0 / ebeta ** 2
                        bsum += wb * beta
                        wbsum += wb

                else:
                    raise hcam.HipercamError(
                        (
                            f"CCD {cnam}, reference aperture {apnam}"
                            f", peak = {height:.1f}"
                            f" < {apsec['fit_height_min_ref']:.1f}"
                        )
                    )

    except hcam.HipercamError as err:
        # trap problems during the fits
        print(
            f"CCD {cnam}, reference aperture {apnam},"
            f" fit failed; extraction aborted. Error = {err}",
            file=sys.stderr,
        )

        for apnam, aper in ccdaper.items():
            store[apnam] = {
                "xe": NaN,
                "ye": NaN,
                "fwhm": NaN,
                "fwhme": NaN,
                "beta": NaN,
                "betae": NaN,
                "dx": 0.,
                "dy": 0.,
            }

        store["mfwhm"] = NaN
        store["mbeta"] = NaN
        return False

    # at this point we are done with measuring the positions of the reference
    # apertures, although their positions have yet to be fixed because another
    # test might have to be passed, namely the differential shift when there
    # is more than one reference star

    if ref:

        # if more than one reference aperture shift is measured, check all
        # shifts for consistency as a guard against cosmic rays and other
        # offsets
        if len(shifts) > 1:
            for n, (apn1, x1, y1) in enumerate(shifts[:-1]):
                for (apn2, x2, y2) in shifts[n + 1 :]:
                    diff = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

                    if np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) > apsec["fit_diff"]:
                        # have exceeded threshold differential shift
                        for apnam, aper in ccdaper.items():
                            store[apnam] = {
                                "xe": NaN,
                                "ye": NaN,
                                "fwhm": NaN,
                                "fwhme": NaN,
                                "beta": NaN,
                                "betae": NaN,
                                "dx": 0.,
                                "dy": 0.,
                            }

                        store["mfwhm"] = NaN
                        store["mbeta"] = NaN

                        # give some info on the problem
                        print(
                            (
                                "CCD {:s}: reference aperture differential "
                                "shift = {:.2f} exceeded limit = {:.2f}"
                            ).format(cnam, diff, apsec["fit_diff"]),
                            file=sys.stderr,
                        )
                        for apn, xs, ys in shifts:
                            print(
                                (
                                    "  Aperture {:s}: xshift, yshift = {:.2f}, {:.2f}"
                                ).format(apn, xs, ys),
                                file=sys.stderr,
                            )
                        return False

        if wxsum > 0.0 and wysum > 0.0:
            # things are OK if we get here. Work out mean shift
            xshift = xsum / wxsum
            yshift = ysum / wysum

            # Finally apply individual shifts to the reference targets this
            # means we do *not* change any positions if the differential shift
            # test fails.
            for apnam, aper in ccdaper.items():
                if aper.ref:
                    aper.x += store[apnam]["dx"]
                    aper.y += store[apnam]["dy"]

        else:

            # all reference fits have failed. Set all others to bad values
            # and return
            for apnam, aper in ccdaper.items():
                if not aper.ref:
                    store[apnam] = {
                        "xe": NaN,
                        "ye": NaN,
                        "fwhm": NaN,
                        "fwhme": NaN,
                        "beta": NaN,
                        "betae": NaN,
                        "dx": 0.,
                        "dy": 0.,
                    }
            store["mfwhm"] = NaN
            store["mbeta"] = NaN

            print(
                f"CCD {cnam}: no reference aperture "
                "fit was successful; extraction aborted",
                file=sys.stderr,
            )
            return False

    else:
        # no reference apertures. All individual
        xshift, yshift = 0.0, 0.0

    # now go over non-reference, non-linked apertures. Failed reference
    # apertures are shifted by the mean shift just determined
    for apnam, aper in ccdaper.items():

        if aper.ref and store[apnam]["xe"] <= 0.0:

            # Move failed reference fit to the mean shift
            aper.x += xshift
            aper.y += yshift
            store[apnam]["dx"] = xshift
            store[apnam]["dy"] = yshift

        elif not aper.ref and not aper.linked:

            # extract Window for data, rflat and flat
            wnam = ccdwin[apnam]
            wdata = ccd[wnam]
            wbias = bccd[wnam]
            wread = read[wnam]
            wgain = gain[wnam]

            try:

                # extract search sub-window around start position.
                swdata = wdata.window(
                    aper.x + xshift - shbox,
                    aper.x + xshift + shbox,
                    aper.y + yshift - shbox,
                    aper.y + yshift + shbox,
                )

                if ref:
                    # if we have reference stars, use the shift from the
                    # the reference stars rather than attempting a search
                    # around the target
                    xold = x = aper.x + xshift
                    yold = y = aper.y + yshift
                    peak = swdata.data[
                        int(round(swdata.y_pixel(y))), int(round(swdata.x_pixel(x)))
                    ]

                else:
                    # no reference star, carry out an explicit search
                    x, y, peak = swdata.search(
                        apsec["search_smooth_fwhm"],
                        aper.x + xshift,
                        aper.y + yshift,
                        apsec["fit_height_min_nrf"],
                        apsec["search_smooth_fft"],
                    )

                # now for a more refined fit. First extract fit Window
                fhbox = apsec["fit_half_width"]
                fwdata = wdata.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)
                fwbias = wbias.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)
                fwread = wread.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)
                fwgain = wgain.window(x - fhbox, x + fhbox, y - fhbox, y + fhbox)

                sky = np.percentile(fwdata.data, 50)

                # get some parameters from previous run where possible
                fit_fwhm = (
                    store[apnam]["fwhm"]
                    if apnam in store and store[apnam]["fwhme"] > 0.0
                    else apsec["fit_fwhm"]
                )

                fit_beta = (
                    store[apnam]["beta"]
                    if apnam in store and store[apnam]["betae"] > 0.0
                    else apsec["fit_beta"]
                )

                # limit the initial value of beta because of tendency
                # to wander to high values and never come down.
                fit_beta = min(fit_beta, apsec["fit_beta_max"])

                # This is where the pure-debiassed data are used
                sigma = np.sqrt(fwread.data**2+np.maximum(0,fwbias.data)/fwgain.data)

                # finally, attempt to fit the target profile
                (
                    (sky, height, x, y, fwhm, beta),
                    (esky, eheight, ex, ey, efwhm, ebeta),
                    extras,
                ) = hcam.fitting.combFit(
                    fwdata, sigma,
                    rfile.method,
                    sky,
                    peak - sky,
                    x,
                    y,
                    fit_fwhm,
                    apsec["fit_fwhm_min"],
                    apsec["fit_fwhm_fixed"],
                    fit_beta,
                    apsec["fit_beta_max"],
                    False,
                    apsec["fit_thresh"],
                    apsec["fit_ndiv"],
                )

                # check the position
                if swdata.distance(x, y) < 0.5:
                    error_message = (
                        "Fitted position ({:.1f},{:.1f}) too close to or"
                        " beyond edge of search window = {!s}"
                    ).format(x, y, swdata.format(True))

                    raise hcam.HipercamError(error_message)

                if fwdata.distance(x, y) < 0.5:
                    error_message = (
                        "Fitted position ({:.1f},{:.1f}) too close to or"
                        " beyond edge of fit window = {!s}".format(
                            x, y, swdata.format(True)
                        )
                    )

                    raise hcam.HipercamError(error_message)

                # check for overly large shifts in the case that we have
                # reference apertures
                if ref:
                    shift = np.sqrt((x - xold) ** 2 + (y - yold) ** 2)
                    if shift > apsec["fit_max_shift"]:
                        error_message = (
                            "Position of non-reference aperture"
                            " shifted by {:.1f} which exceeds "
                            "fit_max_shift = {:.1f}"
                        ).format(shift, apsec["fit_max_shift"])

                        raise hcam.HipercamError(error_message)

                if height >= apsec["fit_height_min_nrf"]:
                    # Height check is important here since no search has
                    # taken place if there are reference stars. If position has
                    # passed checks and the height is OK, then these are useful
                    # data to update parameters
                    store[apnam] = {
                        "xe": ex,
                        "ye": ey,
                        "fwhm": fwhm,
                        "fwhme": efwhm,
                        "beta": beta,
                        "betae": ebeta,
                        "dx": x - aper.x,
                        "dy": y - aper.y
                    }

                    if ref:
                        # apply a fraction 'fit_alpha' times the
                        # change is position relative to the expected
                        # position.
                        aper.x = xold + apsec["fit_alpha"] * (x - xold)
                        aper.y = yold + apsec["fit_alpha"] * (y - yold)
                    else:
                        aper.x = x
                        aper.y = y

                    if efwhm > 0.0:
                        # average FWHM computation
                        wf = 1.0 / efwhm ** 2
                        fsum += wf * fwhm
                        wfsum += wf

                    if ebeta > 0.0:
                        # average beta computation
                        wb = 1.0 / ebeta ** 2
                        bsum += wb * beta
                        wbsum += wb

                else:

                    error_message = (
                        "Non-reference aperture peak = {:.1f} < {:.1f}"
                        ).format(height, apsec["fit_height_min_nrf"])

                    raise hcam.HipercamError(error_message)

            except (hcam.HipercamError, IndexError) as err:
                print(
                    "CCD {:s}, aperture {:s}, fit failed (but extraction may still occur), error = {!s}".format(
                        cnam, apnam, err
                    ),
                    file=sys.stderr,
                )
                aper.x += xshift
                aper.y += yshift

                store[apnam] = {
                    "xe": NaN,
                    "ye": NaN,
                    "fwhm": NaN,
                    "fwhme": NaN,
                    "beta": NaN,
                    "betae": NaN,
                    "dx": xshift,
                    "dy": yshift,
                }

    # finally the linked ones
    for apnam, aper in ccdaper.items():

        if aper.linked:
            aper.x += store[aper.link]["dx"]
            aper.y += store[aper.link]["dy"]

            store[apnam] = {
                "xe": NaN,
                "ye": NaN,
                "fwhm": NaN,
                "fwhme": NaN,
                "beta": NaN,
                "betae": NaN,
                "dx": store[aper.link]["dx"],
                "dy": store[aper.link]["dy"],
            }

    # update mfwhm and mbeta if we can. If not, set them to -1 as a flag down
    # the line that there is no measured value of these parameters. If things
    # work, store them in the config file for retrieval if mfwhm and mbeta go
    # wrong later.
    if wfsum > 0.0:
        store["mfwhm"] = fsum / wfsum
        apsec["fit_fwhm"] = store["mfwhm"]
    else:
        store["mfwhm"] = -1

    if wbsum > 0.0:
        store["mbeta"] = bsum / wbsum
        apsec["fit_beta"] = store["mbeta"]
    else:
        store["mbeta"] = -1

    return True


class Panel:
    """
    Keeps track of the configuration of particular panels of plots so that
    they can be easily re-plotted and selected for additional plotting if
    need be.
    """

    def __init__(
        self,
        device,
        xv1,
        xv2,
        yv1,
        yv2,
        xlabel,
        ylabel,
        tlabel,
        xopt,
        yopt,
        x1,
        x2,
        y1,
        y2,
    ):
        """
        This takes all the arguments needs to set up some axes
        at an arbitrary location, using PGPLOT commands pgsvp, pgswin,
        pgbox and pglab. 'device' is the hipercam.pgp.Device to use
        for the plot. It only stores these values. 'plot' actually
        draws the axes.
        """
        self.device = device
        self.xv1 = xv1
        self.xv2 = xv2
        self.yv1 = yv1
        self.yv2 = yv2
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.tlabel = tlabel
        self.xopt = xopt
        self.yopt = yopt
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2

        # Plot it
        self.plot()

    def plot(self):
        """Plot the Panel. After running this you can plot to the Panel. If you plot
        to another Panel, you can return to this one using 'select'.

        """

        # select the device
        self.device.select()

        # draw the axes
        pgsci(hcam.pgp.Params["axis.ci"])
        pgsch(hcam.pgp.Params["axis.number.ch"])
        pgsvp(self.xv1, self.xv2, self.yv1, self.yv2)

        # avoid invalid limits warnings from PGPLOT
        if self.x1 == self.x2:
            x1, x2 = 0, 1
        else:
            x1, x2 = self.x1, self.x2
        if self.y1 == self.y2:
            y1, y2 = 0, 1
        else:
            y1, y2 = self.y1, self.y2
        pgswin(x1, x2, y1, y2)

        pgbox(self.xopt, 0, 0, self.yopt, 0, 0)

        # plot the labels
        pgsci(hcam.pgp.Params["axis.label.ci"])
        pgsch(hcam.pgp.Params["axis.label.ch"])
        pglab(self.xlabel, self.ylabel, self.tlabel)

    def select(self):
        """
        Selects this panel as the one to plot to. You can only use this if you
        have plotted the panel.
        """

        # select the device
        self.device.select()

        # set the physical scales and the viewport
        pgsvp(self.xv1, self.xv2, self.yv1, self.yv2)

        # avoid invalid limits warnings from PGPLOT
        if self.x1 == self.x2:
            x1, x2 = 0, 1.0
        else:
            x1, x2 = self.x1, self.x2
        if self.y1 == self.y2:
            y1, y2 = 0, 1.0
        else:
            y1, y2 = self.y1, self.y2
        pgswin(x1, x2, y1, y2)


def plotLight(panel, tzero, results, rfile, tkeep, lbuffer):
    """Plots one set of results in the light curve panel. Handles storage of
    points for future re-plots, computing new plot limits where needed.

    It returns a bool which if True means that there is a need to re-plot, and
    a maximum time of a plotted point for use in adjusting the X-axis.  The
    latter is returned as None if no plot is plotted.  This is deferred since
    there may be other panels to be adjusted as well since if a plot is
    cleared, everything has to be re-built.

    """

    # by default, don't re-plot
    replot = False

    skipbadt = rfile["general"]["skipbadt"]

    # shorthand
    sect = rfile["light"]

    # select the light curve panel
    panel.select()

    # Get the current y-range
    ymin, ymax = panel.y1, panel.y2
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    # add points to the plot and buffers, tracking the minimum and
    # maximum values
    tmax = fmin = fmax = None
    for lc in lbuffer:
        tmax, fmin, fmax = lc.add_point(results, tzero, tmax, fmin, fmax, skipbadt)
        lc.trim(tkeep)

    if (
        not sect["y_fixed"]
        and fmin is not None
        and ((ymin == ymax) or (fmin < ymin or fmax > ymax))
    ):

        # we are going to have to replot because we have moved
        # outside the y-limits of the panel. We extend a little bit
        # more than necessary according to extend_y in order to
        # reduce the amount of such re-plotting
        replot = True

        if ymin == ymax:
            # First time through just use data
            extend = sect["extend_y"] * (fmax - fmin)
            # for one CCD, fmin==fmax after first frame,
            # and we never get out of ymin==ymax and the scale
            # is never reset
            if extend == 0:
                extend = 0.001 * fmax if sect["linear"] else 0.001
            ymin = fmin - extend
            ymax = fmax + extend
        else:
            # subsequently use the plot range
            extend = sect["extend_y"] * (ymax - ymin)
            if fmin < ymin:
                ymin = fmin - extend
            if fmax > ymax:
                ymax = fmax + extend

        if sect["linear"]:
            panel.y1, panel.y2 = ymin, ymax
        else:
            panel.y1, panel.y2 = ymax, ymin

    return (replot, tmax)


def plotPosition(xpanel, ypanel, tzero, results, rfile, tkeep, xbuffer, ybuffer):
    """Plots one set of results in the seeing panel. Handles storage of
    points for future re-plots, computing new plot limits where needed.

    It returns a bool which if True means that there is a need to re-plot.
    This is deferred since there may be other panels to be adjusted as well
    since if a plot is cleared, everything has to be re-built.
    """

    # by default, don't re-plot
    replot = False

    skipbadt = rfile["general"]["skipbadt"]

    # shorthand
    sect = rfile["position"]

    # select the X panel
    xpanel.select()

    # add points to the plot and buffers
    tmax = xmin = xmax = None
    for xpos in xbuffer:
        tmax, xmin, xmax = xpos.add_point(results, tzero, tmax, xmin, xmax, skipbadt)
        xpos.trim(tkeep)

    if (
        xmin is not None
        and (xmin < xpanel.y1 or xmax > xpanel.y2)
        and not sect["x_fixed"]
    ):
        replot = True
        extend = sect["extend_y"] * (xpanel.y2 - xpanel.y1)
        if xmin < xpanel.y1:
            xpanel.y1 = xmin - extend

        if xmax > xpanel.y2:
            xpanel.y2 = xmax + extend

    # select the panel
    ypanel.select()

    # add points to the plot and buffers
    ymin = ymax = None
    for ypos in ybuffer:
        tmax, ymin, ymax = ypos.add_point(results, tzero, tmax, ymin, ymax, skipbadt)
        ypos.trim(tkeep)

    if (
        ymin is not None
        and (ymin < ypanel.y1 or ymax > ypanel.y2)
        and not sect["y_fixed"]
    ):
        replot = True
        extend = sect["extend_y"] * (ypanel.y2 - ypanel.y1)
        if ymin < ypanel.y1:
            ypanel.y1 = ymin - extend

        if ymax > ypanel.y2:
            ypanel.y2 = ymax + extend

    return (replot, tmax)


def plotTrans(panel, tzero, results, rfile, tkeep, tbuffer):
    """Plots one set of results in the transmission panel. Handles storage of
    points for future re-plots, computing new plot limits where needed.

    It returns a bool which if True means that there is a need to re-plot.
    This is deferred since there may be other panels to be adjusted as well
    since if a plot is cleared, everything has to be re-built.
    """

    # by default, don't re-plot
    replot = False

    skipbadt = rfile["general"]["skipbadt"]

    # select the light curve panel
    panel.select()

    # add points to the plot and buffers, resetting the maximum
    # transmission where necessary
    tmax = None
    for trans in tbuffer:
        tmax, fmax = trans.add_point(results, tzero, tmax, skipbadt)
        if fmax is not None and fmax > panel.y2:
            trans.fmax *= fmax / 100
            replot = True
        trans.trim(tkeep)

    return (replot, tmax)


def plotSeeing(panel, tzero, results, rfile, tkeep, sbuffer):
    """Plots one set of results in the seeing panel. Handles storage of
    points for future re-plots, computing new plot limits where needed.

    It returns a bool which if True means that there is a need to re-plot.
    This is deferred since there may be other panels to be adjusted as well
    since if a plot is cleared, everything has to be re-built.
    """

    # by default, don't re-plot
    replot = False

    skipbadt = rfile["general"]["skipbadt"]

    # shorthand
    sect = rfile["seeing"]

    # select the light curve panel
    panel.select()

    # add points to the plot and buffers, resetting the maximum
    # transmission where necessary
    tmax = fmax = None
    for see in sbuffer:
        tmax, fmax = see.add_point(results, tzero, tmax, fmax, skipbadt)
        see.trim(tkeep)

    if fmax is not None and fmax > panel.y2 and not sect["y_fixed"]:
        replot = True
        panel.y2 = (1 + sect["extend_y"]) * fmax

    return (replot, tmax)


class BaseBuffer:
    """
    Base class for buffer classes to define a few things in common.
    Container for light curves so they can be re-plotted as they come in.
    There should be one of these per plot line in the 'light' section.
    """

    def __init__(self, plot_config):
        self.cnam = plot_config["ccd"]
        self.targ = plot_config["targ"]
        self.dcol = plot_config["dcol"]
        self.ecol = plot_config["ecol"]
        self.t = []
        self.f = []
        self.fe = []
        self.symb = []

    def trim(self, tkeep):
        """
        Trims points more than tkeep minutes before the last point.
        Assumes times rise monotonically. Nothing done if tkeep <= 0,
        or if the first point does not exceed tkeep by at least 1.
        The idea is to avoid doing this too often. Returns the first time
        """
        if len(self.t) > 1 and tkeep > 0.0 and self.t[-1] > self.t[0] + tkeep + 1:
            tlast = self.t[-1]
            for ntrim, t in enumerate(self.t):
                if t > tlast - tkeep:
                    break
            if ntrim:
                self.t = self.t[ntrim:]
                self.f = self.f[ntrim:]
                self.fe = self.fe[ntrim:]
                self.symb = self.symb[ntrim:]

    def tstart(self):
        """
        Returns the start time of the buffer, 0 if one is not defined
        """
        if len(self.t) == 0:
            return 0.0
        else:
            return self.t[0]


class LightCurve(BaseBuffer):
    """
    Container for light curves so they can be re-plotted as they come in.
    There should be one of these per plot line in the 'light' section.
    """

    def __init__(self, plot_config, linear):
        super().__init__(plot_config)
        self.comp = plot_config["comp"]
        self.off = plot_config["off"]
        self.fac = plot_config["fac"]
        self.linear = linear

    def add_point(self, results, tzero, tmax, fmin, fmax, skipbadt):
        """Extracts the data to be plotted on the light curve plot for the
        given the results of a CCD group reduction as returned by a
        ProcessCCDs object. Assuming all is OK (errors > 0 for both
        comparison and target), it stores (t, f, fe) for possible
        re-plotting, plots the point and returns the value of 'f'
        plotted to help with re-scaling or None if nothing was
        plotted.

        Arguments::

           results : list of results for each CCD group
              see reduce.py for usage, ProcessCCDs for more insight

           tzero : float
              time offset zero point

           tmax : float
              the maximum time plotted so far (minutes). None if not yet set.

           fmin : float
              minimum flux plotted so far. None if not yet set.

           fmax : float
              maximum flux plotted so far. None if not yet set.

           skipbadt : bool
              controls whether points with bad times are skipped or included.
              Unfortunately, hipercam flags times as bad even when they are not.

        Returns (potentially) new values of (tmax, fmin, fmax), according to
        what was plotted.

        """

        for cnam, res in results:
            # loop over the CCD groups. We are only interested in the
            # one that matches the CCD defined for this buffer
            if cnam == self.cnam:
                break
        else:
            return (tmax, fmin, fmax)

        for nframe, store, ccdaper, reses, mjdint, mjdfrac, mjdok, expose in res:

            # loop over each frame of the group

            if mjdok or not skipbadt:

                targ = reses[self.targ]
                ft = targ["counts"]
                fte = targ["countse"]
                saturated = targ["flag"] & hcam.TARGET_SATURATED

                if fte == fte:

                    if self.comp != "!":
                        comp = reses[self.comp]
                        fc = comp["counts"]
                        fce = comp["countse"]
                        saturated |= comp["flag"] & hcam.TARGET_SATURATED

                        if fc == fc and fce == fce:
                            f = ft / fc
                            fe = np.sqrt(
                                (fte / fc) ** 2 + (ft * fce / fc ** 2) ** 2
                            )
                        else:
                            continue
                    else:
                        f = ft
                        fe = fte
                else:
                    continue

                if not self.linear:
                    if f <= 0.0:
                        continue

                    fe = 2.5 / np.log(10) * (fe / f)
                    f = -2.5 * np.log10(f)

                # apply scaling factor and offset
                f *= self.fac
                fe *= self.fac
                f += self.off

                # compute time in terms of minutes from the zeropoint
                tmins = hcam.DMINS * (mjdint + mjdfrac - tzero)

                # OK, we are done for this frame. Store new point
                self.t.append(tmins)
                self.f.append(f)
                self.fe.append(fe)
                if saturated:
                    # mark saturated data with cross
                    self.symb.append(5)
                else:
                    # blob if OK
                    self.symb.append(17)

                # Plot the point in minutes from start point
                pgsch(0.5)
                if self.ecol is not None:
                    pgsci(self.ecol)
                    pgmove(tmins, f - fe)
                    pgdraw(tmins, f + fe)

                pgsci(self.dcol)
                pgpt1(tmins, f, self.symb[-1])

                # track the limits
                tmax = tmins if tmax is None else max(tmins, tmax)
                fmin = f if fmin is None else min(f, fmin)
                fmax = f if fmax is None else max(f, fmax)

        return (tmax, fmin, fmax)


class Xposition(BaseBuffer):
    """Container for X measurements so they can be re-plotted as they come in.
    There should be one of these per plot line in the 'position' section.

    """

    def __init__(self, plot_config):
        super().__init__(plot_config)
        self.xzero = None

    def add_point(self, results, tzero, tmax, xmin, xmax, skipbadt):
        """
        Extracts the data to be plotted on the position plot.

        Arguments::

           results : list of results for each CCD group
              see reduce.py for usage, ProcessCCDs for more insight

           tzero : float
              time offset zero point

           tmax : float
              the maximum time plotted so far (minutes). None if not yet set.

           xmin : float
              minimum X plotted so far. None if not yet set.

           xmax : float
              maximum X plotted so far. None if not yet set.

        Returns (potentially) new values of (tmax, xmin, xmax), according to
        what was plotted.
        """

        for cnam, res in results:
            # loop over the CCD groups. We are only interested in the
            # one that matches the CCD defined for this buffer
            if cnam == self.cnam:
                break
        else:
            return (tmax, xmin, xmax)

        for nframe, store, ccdaper, reses, mjdint, mjdfrac, mjdok, expose in res:

            if mjdok or not skipbadt:

                # loop over each frame of the group

                targ = reses[self.targ]
                x = targ["x"]
                xe = targ["xe"]

                if xe != xe:
                    # skip junk
                    continue

                if self.xzero is None:
                    # initialise the start X position
                    self.xzero = x

                x -= self.xzero

                # compute time in terms of minutes from the zeropoint
                tmins = hcam.DMINS * (mjdint + mjdfrac - tzero)

                # Store new point
                self.t.append(tmins)
                self.f.append(x)
                self.fe.append(xe)
                if targ["flag"] & hcam.TARGET_SATURATED:
                    # mark saturated data with cross
                    self.symb.append(5)
                else:
                    # blob if OK
                    self.symb.append(17)

                # Plot the point in minutes from start point
                pgsch(0.5)
                if self.ecol is not None:
                    pgsci(self.ecol)
                    pgmove(tmins, x - xe)
                    pgdraw(tmins, x + xe)
                pgsci(self.dcol)
                pgpt1(tmins, x, self.symb[-1])

                # track the limits
                tmax = tmins if tmax is None else max(tmins, tmax)
                xmin = x if xmin is None else min(x, xmin)
                xmax = x if xmax is None else max(x, xmax)

        return (tmax, xmin, xmax)


class Yposition(BaseBuffer):
    """Container for Y measurements so they can be re-plotted as they come in.
    There should be one of these per plot line in the 'position' section.

    """

    def __init__(self, plot_config):
        super().__init__(plot_config)
        self.yzero = None

    def add_point(self, results, tzero, tmax, ymin, ymax, skipbadt):
        """
        Extracts the data to be plotted on the position plot.

        Arguments::

           results : list of results for each CCD group
              see reduce.py for usage, ProcessCCDs for more insight

           tzero : float
              time offset zero point

           tmax : float
              the maximum time plotted so far (minutes). None if not yet set.

           ymin : float
              minimum Y plotted so far. None if not yet set.

           ymax : float
              maximum Y plotted so far. None if not yet set.

        Returns (potentially) new values of (tmax, ymin, ymax), according to
        what was plotted.
        """

        for cnam, res in results:
            # loop over the CCD groups. We are only interested in the
            # one that matches the CCD defined for this buffer
            if cnam == self.cnam:
                break
        else:
            return (tmax, ymin, ymax)

        for nframe, store, ccdaper, reses, mjdint, mjdfrac, mjdok, expose in res:

            if mjdok or not skipbadt:

                # loop over each frame of the group

                targ = reses[self.targ]
                y = targ["y"]
                ye = targ["ye"]

                if ye != ye:
                    # skip junk
                    continue

                if self.yzero is None:
                    # initialise the start Y position
                    self.yzero = y

                y -= self.yzero

                # compute time in terms of minutes from the zeropoint
                tmins = hcam.DMINS * (mjdint + mjdfrac - tzero)

                # Store new point
                self.t.append(tmins)
                self.f.append(y)
                self.fe.append(ye)
                if targ["flag"] & hcam.TARGET_SATURATED:
                    # mark saturated data with cross
                    self.symb.append(5)
                else:
                    # blob if OK
                    self.symb.append(17)

                # Plot the point in minutes from start point
                pgsch(0.5)
                if self.ecol is not None:
                    pgsci(self.ecol)
                    pgmove(tmins, y - ye)
                    pgdraw(tmins, y + ye)
                pgsci(self.dcol)
                pgpt1(tmins, y, self.symb[-1])

                # track the limits
                tmax = tmins if tmax is None else max(tmins, tmax)
                ymin = y if ymin is None else min(y, ymin)
                ymax = y if ymax is None else max(y, ymax)

        return (tmax, ymin, ymax)


class Transmission(BaseBuffer):
    """
    Container for light curves so they can be re-plotted as they come in.
    There should be one of these per plot line in the 'light' section.
    """

    def __init__(self, plot_config):
        super().__init__(plot_config)
        self.fmax = None  # Maximum to scale the transmission

    def add_point(self, results, tzero, tmax, skipbadt):
        """
        Extracts the data to be plotted on the transmission plot for given the
        time and the results (for all CCDs, as returned by
        extractFlux. Assuming all is OK (errors > 0), it stores (t, f, fe) for
        possible re-plotting, plots the point and returns the value of 'f'
        plotted to help with re-scaling or None if nothing was plotted.  't'
        is the time in minutes since the start of the run
        """
        fmax = None
        for cnam, res in results:
            # loop over the CCD groups. We are only interested in the
            # one that matches the CCD defined for this buffer
            if cnam == self.cnam:
                break
        else:
            return (tmax, fmax)

        for nframe, store, ccdaper, reses, mjdint, mjdfrac, mjdok, expose in res:

            if mjdok or not skipbadt:

                # loop over each frame of the group
                targ = reses[self.targ]
                f = targ["counts"]
                fe = targ["countse"]

                if f <= 0 or fe != fe:
                    # skip junk
                    continue

                # compute time in terms of minutes from the zeropoint
                tmins = hcam.DMINS * (mjdint + mjdfrac - tzero)

                # Store new point
                self.t.append(tmins)
                self.f.append(f)
                self.fe.append(fe)
                if targ["flag"] & hcam.TARGET_SATURATED:
                    # mark saturated data with cross
                    self.symb.append(5)
                else:
                    # blob if OK
                    self.symb.append(17)

                if self.fmax is None:
                    # initialise the maximum flux
                    self.fmax = f

                f *= 100 / self.fmax
                fe *= 100 / self.fmax

                # Plot the point in minutes from start point
                pgsch(0.5)
                if self.ecol is not None:
                    pgsci(self.ecol)
                    pgmove(tmins, f - fe)
                    pgdraw(tmins, f + fe)
                pgsci(self.dcol)
                pgpt1(tmins, f, self.symb[-1])

                # track the limits
                tmax = tmins if tmax is None else max(tmins, tmax)
                fmax = f if fmax is None else max(f, fmax)

        return (tmax, fmax)


class Seeing(BaseBuffer):
    """
    Container for light curves so they can be re-plotted as they come in.
    There should be one of these per plot line in the 'light' section.
    """

    def __init__(self, plot_config, scale):
        super().__init__(plot_config)
        self.scale = scale

    def add_point(self, results, tzero, tmax, fmax, skipbadt):
        """
        Extracts the data to be plotted on the transmission plot for given the
        time and the results (for all CCDs, as returned by
        extractFlux. Assuming all is OK (errors > 0), it stores (t, f, fe) for
        possible re-plotting, plots the point and returns the value of 'f'
        plotted to help with re-scaling or None if nothing was plotted.  't'
        is the time in minutes since the start of the run
        """
        for cnam, res in results:
            # loop over the CCD groups. We are only interested in the
            # one that matches the CCD defined for this buffer
            if cnam == self.cnam:
                break
        else:
            return (tmax, fmax)

        for nframe, store, ccdaper, reses, mjdint, mjdfrac, mjdok, expose in res:

            if mjdok or not skipbadt:

                # loop over each frame of the group
                targ = reses[self.targ]
                f = targ["fwhm"]
                fe = targ["fwhme"]

                if f <= 0 or fe != fe:
                    # skip junk
                    continue

                f *= self.scale
                fe *= self.scale

                # compute time in terms of minutes from the zeropoint
                tmins = hcam.DMINS * (mjdint + mjdfrac - tzero)

                # Store new point
                self.t.append(tmins)
                self.f.append(f)
                self.fe.append(fe)
                if targ["flag"] & hcam.TARGET_SATURATED:
                    # mark saturated data with cross
                    self.symb.append(5)
                else:
                    # blob if OK
                    self.symb.append(17)

                # Plot the point in minutes from start point
                pgsch(0.5)
                if self.ecol is not None:
                    pgsci(self.ecol)
                    pgmove(tmins, f - fe)
                    pgdraw(tmins, f + fe)
                pgsci(self.dcol)
                pgpt1(tmins, f, self.symb[-1])

                # track the limits
                tmax = tmins if tmax is None else max(tmins, tmax)
                fmax = f if fmax is None else max(f, fmax)

        return (tmax, fmax)


def ctrans(cname):
    """
    Translates a colour name (cname) into a PGPLOT index
    which is the return value. Defaults to index 1 and prints
    a message if cname not recognised. Returns None if cname is None
    """

    if cname == "!":
        return None
    elif cname in hcam.CNAMS:
        return hcam.CNAMS[cname]
    else:
        rnames = list(hcam.CNAMS.keys())
        rnames.sort()

        warnings.warn(
            "Failed to recognize colour = '{:s}'; defaulting to black.\n"
            "Recognised colours are {:s}\n".format(
                cname, ", ".join(["'{:s}'".format(name) for name in rnames])
            )
        )
        return 1


def toBool(rfile, section, param):
    """Converts yes / no / true / false string responses into boolean
    True / False. This is used a few times in the code to read a
    reduce file.

    Arguments::

       rfile : Rfile
         the reduce file, an Odict of Odicts

      section : str
         the section name

      param : str
         the parameter

    Returns nothing; rfile modified on exit. A ValueError is
    raised if the initial value is neither 'yes', 'no', 'true'
    of 'false' (case-independent)

    """

    lstr = rfile[section][param].lower()

    if lstr == "yes" or lstr == "true":
        rfile[section][param] = True
    elif lstr == "no" or lstr == "false":
        rfile[section][param] = False
    else:
        raise hcam.HipercamError(
            f"{section}.{param} is a boolean; 'yes', 'no', 'true' and 'false'"
            " (case-insensitive) are the only supported values. It's current"
            f" (invalid) value = {rfile[section][param]}; please correct it."
        )


class LogWriter:
    """Context manager to handle opening logfiles, writing headers to logfiles
    and safe exit.

    """

    def __init__(self, filename, rfile, hipercam_version, plist):
        self.rfile = rfile
        self.filename = filename
        self.hipercam_version = hipercam_version
        self.plist = plist
        self.toffset = rfile["general"]["toffset"]

    def __enter__(self):
        # open filehandle, write header
        self.log = open(self.filename, "w")
        self.write_header()
        return self

    def __exit__(self, *args):
        self.log.close()

    def write_results(self, results):
        """
        Append results to open logfile.

        Arguments::

           results : list of 2-element tuples.
              each tuple is composed of (cnam, list) where cnam is a CCD label
              and list is a list of 8-element tuples each of which is composed
              of (nframe, store, ccdaper, results, mjdint, mjdfrac, mjdok, expose)
              storing all the relevant data for the CCD / exposure in question.
        """
        monitor = self.rfile["monitor"]

        def float2str(value, form, sval=None):
            """Formats a float as a string accounting for NaNs. value is the value,
            form is the format to be used for the string, e.g. '.2f'. If value matches
            the special value sval then also convert to NaN.
            """
            if value == sval or np.isnan(value):
                return 'nan'
            else:
                return f'{value:{form}}'

        alerts = []
        for cnam, res in results:
            # Loop over all CCDs

            for nframe, store, ccdaper, reses, mjdint, mjdfrac, mjdok, expose in res:
                # Loop over all frames for the specific CCD group in question

                # A spacer because the lines are long
                self.log.write("#\n")

                # The time is worked out so that it can especially accurate if a suitable
                # integer offset is applied. The parentheses are important here.
                mjd = (mjdint - self.toffset) + mjdfrac

                # get mean profile parameters
                mfwhm = store["mfwhm"]
                mbeta = store["mbeta"]

                # write generic data. have to take a little care because of NaN values
                self.log.write(
                    "{:s} {:d} {:.14f} {:b} {:.8g} {} {} ".format(
                        cnam, nframe, mjd, mjdok, expose,
                        float2str(mfwhm,'.2f',-1.),
                        float2str(mbeta,'.2f',-1.),
                    )
                )

                # now for data per aperture
                for apnam in ccdaper:
                    r = reses[apnam]
                    self.log.write(
                        "{} {} {} {} {} {} {} {} {} {} {} {} "
                        "{:d} {:d} {:d} {:d} ".format(
                            float2str(r["x"],'.4f'),
                            float2str(r["xe"],'.4f'),
                            float2str(r["y"],'.4f'),
                            float2str(r["ye"],'.4f'),
                            float2str(r["fwhm"],'.3f'),
                            float2str(r["fwhme"],'.3f',-1.),
                            float2str(r["beta"],'.3f'),
                            float2str(r["betae"],'.3f'),
                            float2str(r["counts"],'.1f'),
                            float2str(r["countse"],'.1f'),
                            float2str(r["sky"],'.2f'),
                            float2str(r["skye"],'.2f'),
                            r["nsky"],
                            r["nrej"],
                            r["cmax"],
                            r["flag"],
                        )
                    )

                    if apnam in monitor:
                        # accumulate any problems with particular targets
                        bitmasks = monitor[apnam]
                        flag = r["flag"]
                        messes = []
                        for bitmask in bitmasks:
                            if (flag & bitmask) and bitmask in hcam.FLAG_MESSAGES:
                                messes.append(hcam.FLAG_MESSAGES[bitmask])

                        if len(messes):
                            alerts.append(
                                " *** WARNING: Frame {:d}, CCD {:s}, aperture {:s}: {:s}".format(
                                    nframe, cnam, apnam, ", ".join(messes)
                                )
                            )

                # finish the line
                self.log.write("\n")

        # flush the buffer to ensure that we have complete lines if we hit ctrl-C
        self.log.flush()

        # Return with the accumulated list of alerts.
        return alerts

    def write_header(self):

        # first, a general description
        self.log.write(
            """#
# This is a logfile produced by the HiPERCAM pipeline command 'reduce'. It consists
# of one line per reduced CCD per exposure. Each line contains all the information
# from all apertures defined for the CCD. The column names are defined just before
# the data below. The logfile was produced using version number:
#
#  {hipercam_version}
#
# of the HiPERCAM reduction software, and was generated using the following
# command-line inputs to 'reduce':
#
""".format(
                hipercam_version=self.hipercam_version
            )
        )

        # second, list the command-line inputs to the logfile
        for line in self.plist:
            self.log.write("# {:s}".format(line))

        # third, list the reduce file
        self.log.write(
            """#
# and here is a minimal version of the reduce file used ['rfile' above] with
# all between-line comments removed for compactness:
#
"""
        )

        # skip these as they only affect the on the fly plots,
        # not the final values
        skip_sections = ("lcplot", "light", "transmission", "seeing")
        skip = False
        with open(self.rfile.filename) as fred:
            for line in fred:
                if not line.startswith("#") and not line.isspace():
                    if line.startswith("["):
                        sect = line[1 : line.find("]")].strip()
                        if sect in skip_sections:
                            skip = True
                            continue
                        else:
                            skip = False
                        self.log.write("#\n#   {:s}".format(line))
                    elif not skip:
                        self.log.write("#   {:s}".format(line))

        # fourth, write the apertures
        self.log.write(
            """#
# Next here is the aperture file used in JSON-style format that (without
# the initial comment hashes) is readable by setaper and reduce:
#
"""
        )
        # convert aperture file to JSON-style string, split line by line,
        # pre-pend comment and indentation, write out to the logfile.
        lines = [
            "#   {:s}\n".format(line) for line in self.rfile.aper.toString().split("\n")
        ]
        for line in lines:
            self.log.write(line)

        # fifth the column names for each CCD which has any
        # apertures
        self.log.write(
            """#
# Now follow column name definitions for each CCD. These include all apertures
# of the CCD. Since there are 15 items stored per aperture and each column name
# is built from the item name followed by an underscore and finally the aperture
# name, thus these definition can be very long and won't be very readable. Each
# line starts with the '<CCD label> = '. The various items have the following
# meanings. First the generic ones at the start of the line:
#
#    CCD    : CCD label
#    nframe : integer frame number
#    MJD    : MJD at the centre of the exposure
#    MJDok  : flag to say whether the MJD is thought reliable
#    Exptim : exposure time, seconds
#    mfwhm  : the mean FWHM used to determine the aperture scale (-1 if none)
#    mbeta  : the mean Moffat beta exponent (-1 if none)
#
# Then a set of 16 items that is repeated for each aperture and will
# have '_<apnam>' added to them. The 'errors' are RMS uncertainties, and
# may often be set to -1 if no direct measurement is made. This can happen
# through problems that occur or because an aperture is linked for instance.
#
#    x       : X-position
#    xe      : X-position error
#    y       : Y-position
#    ye      : Y-position error
#    fwhm    : FWHM
#    fwhme   : FWHM error
#    beta    : Moffat exponent
#    betae   : Moffat exponent error
#    counts   : sky-subtracted counts in aperture
#    countse  : error in       "          "
#    sky     : sky level, counts per pixel
#    skye    : sky level error, counts per pixel
#    nsky    : number of contributing sky pixels (those not rejected)
#    nrej    : number of sky pixels rejected
#    cmax    : maximum raw counts in aperture
#    flag    : status flag, 0 = all OK.
#
# Start of column name definitions:
#
"""
        )
        toffset = self.rfile["general"]["toffset"]
        for cnam, ccdaper in self.rfile.aper.items():
            if len(ccdaper) == 0:
                # nothing will be written for CCDs without
                # apertures
                continue

            if toffset == 0:
                cnames = "# {:s} = CCD nframe MJD MJDok Exptim mfwhm mbeta ".format(
                    cnam
                )
            else:
                cnames = "# {:s} = CCD nframe MJD-{:d} MJDok Exptim mfwhm mbeta ".format(
                    cnam, toffset
                )

            for apnam in ccdaper:
                cnames += (
                    "x_{0:s} xe_{0:s} y_{0:s} ye_{0:s} "
                    "fwhm_{0:s} fwhme_{0:s} beta_{0:s} betae_{0:s} "
                    "counts_{0:s} countse_{0:s} sky_{0:s} skye_{0:s} "
                    "nsky_{0:s} nrej_{0:s} cmax_{0:s} flag_{0:s} ".format(apnam)
                )
            self.log.write(cnames + "\n")

        # now the datatypes for building into structured arrays
        self.log.write(
            """#
# End of column name definitions
#
# Now follow a similar series of datatypes which are designed to be used to
# build a numpy.dtype for each CCD to allow them to read into numpy structured
# arrays.
#
# Start of data type definitions:
#
"""
        )
        atypes = "f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 i4 i4 i4 u4 "
        for cnam, ccdaper in self.rfile.aper.items():
            if len(ccdaper) == 0:
                # nothing will be written for CCDs without
                # apertures
                continue

            self.log.write(
                "# {:s} = s i4 f8 ? f4 f4 f4 {:s}\n".format(cnam, len(ccdaper) * atypes)
            )

        self.log.write(
            """#
# End of data type definitions
#
"""
        )

    # that's it for the headers!
