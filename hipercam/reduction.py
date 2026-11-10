"""
This contains code used in common by the script reduce.py
"""

from collections import OrderedDict
import sys
import warnings

import numpy as np

import hipercam as hcam
from hipercam import fitting, utils
from hipercam.ccd import crop_calib_frame_to

NaN = float("NaN")


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
                            "found entry line before any section = \n{:s}".format(line)
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
            v = sect["version"]
            rfv = hcam.REDUCE_FILE_VERSION
            raise hcam.HipercamError(
                f"Version mismatch: reduce file = {v}, reduce = {rfv}. "
                "If your reduce file date is earlier than the reduce version, you may want to try "
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
        sect["warns"] = warns

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
            raise hcam.HipercamError("sky.error must be either 'variance' or 'photon'")

        if skysec["method"] != "clipped" and skysec["method"] != "median":
            raise hcam.HipercamError("sky.method must be either 'clipped' or 'median'")

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
            self.bias = crop_calib_frame_to(mccd, self.bias, "bias")

        if self.dark is not None:
            self.dark = crop_calib_frame_to(mccd, self.dark, "dark")

        if self.flat is not None:
            self.flat = crop_calib_frame_to(mccd, self.flat, "flat")

        if self.fmap is not None:
            self.fmap = crop_calib_frame_to(mccd, self.fmap, "fmap")
            nhalf = self["calibration"]["nhalf"]
            self.fpair = self.fpair.crop(mccd, nhalf)

        if isinstance(self.readout, hcam.MCCD):
            self.readout = crop_calib_frame_to(mccd, self.readout, "readout")

        if isinstance(self.gain, hcam.MCCD):
            self.gain = crop_calib_frame_to(mccd, self.gain, "gain")


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
                self.store[cnam] = {"mfwhm": -1.0, "mbeta": -1.0}

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
                    (
                        nframe,
                        store,
                        ccdaper,
                        results,
                        mjdint,
                        mjdfrac,
                        mjdok,
                        expose,
                    ) = res[-1]
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
                "dx": 0.0,
                "dy": 0.0,
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
    xshift, yshift = 0.0, 0.0

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
                    aper.x + xshift - shbox,
                    aper.x + xshift + shbox,
                    aper.y + yshift - shbox,
                    aper.y + yshift + shbox,
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
                sigma = np.sqrt(
                    fwread.data**2 + np.maximum(0, fwbias.data) / fwgain.data
                )

                # refine the Aperture position by fitting the profile
                (
                    (sky, height, x, y, fwhm, beta),
                    (esky, eheight, ex, ey, efwhm, ebeta),
                    extras,
                ) = hcam.fitting.combFit(
                    fwdata,
                    sigma,
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

                    # reminder: compo apertures shift in opposite sense
                    # store dx, dy as the shift corrected for this flip
                    dx = x - aper.x if not aper.compo else aper.x - x
                    wx = 1.0 / ex**2
                    wxsum += wx
                    xsum += wx * dx

                    dy = y - aper.y if not aper.compo else aper.y - y
                    wy = 1.0 / ey**2
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
                        wf = 1.0 / efwhm**2
                        fsum += wf * fwhm
                        wfsum += wf

                    if ebeta > 0.0:
                        # average beta computation
                        wb = 1.0 / ebeta**2
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
                "dx": 0.0,
                "dy": 0.0,
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
                for apn2, x2, y2 in shifts[n + 1 :]:
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
                                "dx": 0.0,
                                "dy": 0.0,
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
                    if not aper.compo:
                        aper.x += store[apnam]["dx"]
                        aper.y += store[apnam]["dy"]
                    else:
                        aper.x -= store[apnam]["dx"]
                        aper.y -= store[apnam]["dy"]

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
                        "dx": 0.0,
                        "dy": 0.0,
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
            if not aper.compo:
                aper.x += xshift
                aper.y += yshift
            else:
                aper.x -= xshift
                aper.y -= yshift
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
                    if not aper.compo:
                        xold = x = aper.x + xshift
                        yold = y = aper.y + yshift
                    else:
                        xold = x = aper.x - xshift
                        yold = y = aper.y - yshift
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
                sigma = np.sqrt(
                    fwread.data**2 + np.maximum(0, fwbias.data) / fwgain.data
                )

                # finally, attempt to fit the target profile
                (
                    (sky, height, x, y, fwhm, beta),
                    (esky, eheight, ex, ey, efwhm, ebeta),
                    extras,
                ) = hcam.fitting.combFit(
                    fwdata,
                    sigma,
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
                    if aper.compo:
                        dx = aper.x - x
                        dy = aper.y - y
                    else:
                        dx = x - aper.x
                        dy = y - aper.y

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
                        wf = 1.0 / efwhm**2
                        fsum += wf * fwhm
                        wfsum += wf

                    if ebeta > 0.0:
                        # average beta computation
                        wb = 1.0 / ebeta**2
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
                if not aper.compo:
                    aper.x += xshift
                    aper.y += yshift
                else:
                    aper.x -= xshift
                    aper.y -= yshift

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
            # dx, dy in store always saved in the sense of non-injected apertures
            # even if linked aperture is injected, 'store[aper.link]["dx"]' needs
            # inverting for injected apertures.
            if not aper.compo:
                aper.x += store[aper.link]["dx"]
                aper.y += store[aper.link]["dy"]
            else:
                aper.x -= store[aper.link]["dx"]
                aper.y -= store[aper.link]["dy"]

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


def ccdproc(
    extraction_fn, cnam, ccds, bccds, rccds, nframes, read, gain, ccdwin, rfile, store
):
    """Processing steps for a sequential set of images from the same
    CCD. This is designed for parallelising the processing across CCDs
    of multi-arm cameras like ULTRACAM and HiPERCAM using
    multiprocessing. To be called *after* checking that any processing
    is needed.

    Arguments::

       extraction_fn : function
          callable to process a single CCD. For the signature of this function
          see the docstring to extractFlux below.

       cnam : string
          name of CCD, for information purposes (e.g. 'red', '3', etc)

       ccds : List of CCDs
          the CCDs for processing which should have been debiassed, flat
          fielded and multiplied by the gain to get into electrons.

       bccds : List of CCDs
          the CCDs for processing which should just have been debiassed

       rccds : List of CCDs
          unprocessed CCDs, one-to-one correspondence with 'ccds', used
          to measure saturation

       nframes : List of ints
          frame numbers for each CCD

       read : CCD
          readnoise frame divided by the flatfield

       gain : CCD
          gain frame dmultiplied by the flatfield

       ccdwin : dict
          label of the Window enclosing each aperture

       rfile : Rfile object
          reduction control parameters. rfile.aper used to store the aperture
          parameters.

       store : dict
          dictionary of results

    Returns: (cnam, list[res]) where 'res' represents a tuple
    of results for each input CCD and contains the following:

    (nframe, store, ccdaper, results, mjdint, mjdfrac, mjdok, expose)

    """

    # At this point 'ccds' contains a list of CCD each of which
    # contains all the Windows of a CCD, 'ccdaper' all of its
    # apertures, 'ccdwin' the label of the Window enclosing each
    # aperture, 'rfile' contains control parameters, 'rflat' contains
    # the readout noise in electrons and divided by the flat as a CCD,
    # 'store' is a dictionary initially with jus 'mfwhm' and 'mbeta'
    # set = -1, but will pick up extra stuff from moveApers for use by
    # extractFlux along with revised values of mfwhm and mbeta which
    # are used to initialise profile fits next time.

    res = []
    for ccd, bccd, rccd, nframe in zip(ccds, bccds, rccds, nframes):
        # Loop through the CCDs supplied

        # move the apertures. Also fits the profile to the reference stars
        moveApers(cnam, ccd, bccd, read, gain, ccdwin, rfile, store)

        # extract flux from all apertures of each CCD. Return with the CCD
        # name, the store dictionary, ccdaper and then the results from
        # extractFlux for compatibility with multiprocessing.
        results = extraction_fn(cnam, ccd, bccd, rccd, read, gain, ccdwin, rfile, store)

        # Save the essentials
        res.append(
            (
                nframe,
                store,
                rfile.aper[cnam],
                results,
                ccd.head["MJDINT"],
                ccd.head["MJDFRAC"],
                ccd.head.get("GOODTIME", True),
                ccd.head.get("EXPTIME", 1.0),
            )
        )

    return (cnam, res)


def extractFlux(cnam, ccd, bccd, rccd, read, gain, ccdwin, rfile, store):
    """This extracts the flux of all apertures of a given CCD.

    The steps are (1) aperture resizing, (2) sky background estimation, (3)
    flux extraction. The apertures are assumed to be correctly positioned.

    It returns the results as a dictionary keyed on the aperture label. Each
    entry returns a list:

    [x, ex, y, ey, fwhm, efwhm, beta, ebeta, counts, countse, sky, esky,
    nsky, nrej, flag]

    flag = bitmask. See hipercam.core to see all the options which are
    referred to by name in the code e.g. ALL_OK. The various flags can
    signal that there no sky pixels (NO_SKY), the sky aperture was off
    the edge of the window (SKY_AT_EDGE), etc.

    This code::

       >> bset = flag & TARGET_SATURATED

    determines whether the data saturation flag is set for example.

    Arguments::

       cnam : string
          CCD identifier label

       ccd : CCD
           the fully processed CCD (debiassed, dark-corrected, flat-fielded, fringe-corrected)

       bccd : CCD
           debiassed-only CCD. Used in variance computation.

       rccd : CCD
          corresponding raw CCD, used to work out whether data are
          saturated in target aperture.

       read : CCD
           readnoise divided by the flat-field

       gain : CCD
           gain multiplied by the flat field

       ccdwin : dictionary of strings
           the Window label corresponding to each Aperture

       rfile : Rfile
           reduce file configuration parameters

       store : dict of dicts
           see moveApers for what this contains.

    """

    # initialise flag
    flag = hcam.ALL_OK

    ccdaper = rfile.aper[cnam]

    # get the control parameters
    (
        resize,
        extype,
        r1fac,
        r1min,
        r1max,
        r2fac,
        r2min,
        r2max,
        r3fac,
        r3min,
        r3max,
    ) = rfile["extraction"][cnam]

    results = {}
    mfwhm = store["mfwhm"]

    if resize == "variable" or extype == "optimal":
        if mfwhm <= 0:
            # return early here as there is nothing we can do.
            print(
                (
                    " *** WARNING: CCD {:s}: no measured FWHM to re-size"
                    " apertures or carry out optimal extraction; no"
                    " extraction possible"
                ).format(cnam)
            )
            # set flag to indicate no FWHM
            flag = hcam.NO_FWHM

            for apnam, aper in ccdaper.items():
                info = store[apnam]
                results[apnam] = {
                    "x": aper.x,
                    "xe": info["xe"],
                    "y": aper.y,
                    "ye": info["ye"],
                    "fwhm": info["fwhm"],
                    "fwhme": info["fwhme"],
                    "beta": info["beta"],
                    "betae": info["betae"],
                    "counts": NaN,
                    "countse": NaN,
                    "sky": NaN,
                    "skye": NaN,
                    "nsky": 0,
                    "nrej": 0,
                    "flag": flag,
                    "cmax": 0,
                }
            return results

        else:
            # Re-size the apertures
            for aper in ccdaper.values():
                aper.rtarg = max(r1min, min(r1max, r1fac * mfwhm))
                aper.rsky1 = max(r2min, min(r2max, r2fac * mfwhm))
                aper.rsky2 = max(r3min, min(r3max, r3fac * mfwhm))

    elif resize == "fixed":
        # just apply the max and min limits
        for aper in ccdaper.values():
            aper.rtarg = max(r1min, min(r1max, aper.rtarg))
            aper.rsky1 = max(r2min, min(r2max, aper.rsky1))
            aper.rsky2 = max(r3min, min(r3max, aper.rsky2))

    else:
        raise hcam.HipercamError(
            "CCD {:s}: 'variable' and 'fixed' are the only"
            " aperture resizing options".format(cnam)
        )

    # apertures have been positioned in moveApers and now re-sized. Finally
    # we can extract something.
    for apnam, aper in ccdaper.items():
        # initialise flag
        flag = hcam.ALL_OK

        # extract Windows relevant for this aperture
        wnam = ccdwin[apnam]

        wdata = ccd[wnam]
        wbias = bccd[wnam]
        wread = read[wnam]
        wgain = gain[wnam]
        wraw = rccd[wnam]

        # extract sub-windows that include all of the pixels that could
        # conceivably affect the aperture. We have to check that 'extra'
        # apertures do not go beyond rsky2 which would normally be expected to
        # be the default outer radius
        rmax = aper.rsky2
        for xoff, yoff in aper.extra:
            rmax = max(rmax, np.sqrt(xoff**2 + yoff**2) + aper.rtarg)

        # this is the region of interest
        x1, x2, y1, y2 = (
            aper.x - aper.rsky2 - wdata.xbin,
            aper.x + aper.rsky2 + wdata.xbin,
            aper.y - aper.rsky2 - wdata.ybin,
            aper.y + aper.rsky2 + wdata.ybin,
        )

        try:
            # extract sub-Windows
            swdata = wdata.window(x1, x2, y1, y2)
            swbias = wbias.window(x1, x2, y1, y2)
            swread = wread.window(x1, x2, y1, y2)
            swgain = wgain.window(x1, x2, y1, y2)
            swraw = wraw.window(x1, x2, y1, y2)

            # some checks for possible problems. bitmask flags will be set if
            # they are encountered.
            xlo, xhi, ylo, yhi = swdata.extent()
            if (
                xlo > aper.x - aper.rsky2
                or xhi < aper.x + aper.rsky2
                or ylo > aper.y - aper.rsky2
                or yhi < aper.y + aper.rsky2
            ):
                # the sky aperture overlaps the edge of the window
                flag |= hcam.SKY_AT_EDGE

            if (
                xlo > aper.x - aper.rtarg
                or xhi < aper.x + aper.rtarg
                or ylo > aper.y - aper.rtarg
                or yhi < aper.y + aper.rtarg
            ):
                # the target aperture overlaps the edge of the window
                flag |= hcam.TARGET_AT_EDGE

            for xoff, yoff in aper.extra:
                rout = np.sqrt(xoff**2 + yoff**2) + aper.rtarg
                if (
                    xlo > aper.x - rout
                    or xhi < aper.x + rout
                    or ylo > aper.y - rout
                    or yhi < aper.y + rout
                ):
                    # an extra target aperture overlaps the edge of the window
                    flag |= hcam.TARGET_AT_EDGE

            # compute X, Y arrays over the sub-window relative to the centre
            # of the aperture and the distance squared from the centre (Rsq)
            # to save a little effort.
            x = swdata.x(np.arange(swdata.nx)) - aper.x
            y = swdata.y(np.arange(swdata.ny)) - aper.y
            X, Y = np.meshgrid(x, y)
            Rsq = X**2 + Y**2

            # squared aperture radii for comparison
            R1sq, R2sq, R3sq = aper.rtarg**2, aper.rsky1**2, aper.rsky2**2

            # sky selection, accounting for masks and extra (which we assume
            # acts like a sky mask as well)
            sok = (Rsq > R2sq) & (Rsq < R3sq)
            for xoff, yoff, radius in aper.mask:
                sok &= (X - xoff) ** 2 + (Y - yoff) ** 2 > radius**2
            for xoff, yoff in aper.extra:
                sok &= (X - xoff) ** 2 + (Y - yoff) ** 2 > R1sq

            # sky data
            dsky = swdata.data[sok]

            if len(dsky):
                # we have some sky!

                if rfile["sky"]["method"] == "clipped":
                    # clipped mean. Take average, compute RMS,
                    # reject pixels > thresh*rms from the mean.
                    # repeat until no new pixels are rejected.

                    thresh = rfile["sky"]["thresh"]
                    ok = np.ones_like(dsky, dtype=bool)
                    nrej = 1
                    while nrej:
                        slevel = dsky[ok].mean()
                        srms = dsky[ok].std()
                        nold = len(dsky[ok])
                        ok = ok & (np.abs(dsky - slevel) < thresh * srms)
                        nrej = nold - len(dsky[ok])
                        if len(dsky[ok]) == 0:
                            # no sky. will still return the flux in
                            # the aperture but set flag and the sky
                            # uncertainty to -1
                            flag |= hcam.NO_SKY
                            slevel = 0
                            serror = -1
                            nrej = 0
                            break

                    nsky = len(dsky[ok])
                    if nsky:
                        # serror -- error in the sky estimate.
                        serror = srms / np.sqrt(nsky)

                else:
                    # 'median' goes with 'photon'
                    slevel = dsky.median()
                    nsky = len(dsky)
                    nrej = 0

                    # read*gain/flat and flat over sky region
                    dread = swread.data[sok]
                    dgain = swgain.data[sok]

                    serror = np.sqrt(
                        (dread**2 + np.max(0, dsky) / dgain).sum() / nsky**2
                    )

            else:
                # no sky. will still return the flux in the aperture but set
                # flag and the sky uncertainty to -1
                flag |= hcam.NO_SKY
                slevel = NaN
                serror = NaN
                nsky = 0
                nrej = 0

            # size of a pixel which is used to taper pixels as they approach
            # the edge of the aperture to reduce pixellation noise
            size = np.sqrt(wdata.xbin * wdata.ybin)

            # target selection, accounting for extra apertures and allowing
            # pixels to contribute if their centres are as far as size/2 beyond
            # the edge of the circle (but with a tapered weight)
            dok = Rsq < (aper.rtarg + size / 2.0) ** 2

            if not dok.any():
                # check there are some valid pixels
                flag |= hcam.NO_DATA
                raise hcam.HipercamError("no valid pixels in aperture")

            # check for saturation and nonlinearity
            cmax = int(swraw.data[dok].max())
            if cnam in rfile.warn:
                if cmax >= rfile.warn[cnam]["saturation"]:
                    flag |= hcam.TARGET_SATURATED

                if cmax >= rfile.warn[cnam]["nonlinear"]:
                    flag |= hcam.TARGET_NONLINEAR

            else:
                warnings.warn(
                    "CCD {:s} has no nonlinearity or saturation levels set".format(cnam)
                )

            # Pixellation amelioration:
            #
            # The weight of a pixel is set to 1 at the most and then linearly
            # declines as it approaches the edge of the aperture. The scale over
            # which it declines is set by 'size', the geometric mean of the
            # binning factors. A pixel with its centre exactly on the edge
            # gets a weight of 0.5.
            wgt = np.minimum(
                1, np.maximum(0, (aper.rtarg + size / 2.0 - np.sqrt(Rsq)) / size)
            )
            for xoff, yoff in aper.extra:
                rsq = (X - xoff) ** 2 + (Y - yoff) ** 2
                dok |= rsq < (aper.rtarg + size / 2.0) ** 2
                wg = np.minimum(
                    1, np.maximum(0, (aper.rtarg + size / 2.0 - np.sqrt(rsq)) / size)
                )
                wgt = np.maximum(wgt, wg)

            # the values needed to extract the flux.
            dtarg = swdata.data[dok]
            dbtarg = swbias.data[dok]
            dread = swread.data[dok]
            dgain = swgain.data[dok]
            wtarg = wgt[dok]

            # 'override' to indicate we want to override the readout noise.
            if nsky and rfile["sky"]["error"] == "variance":
                # from sky variance
                rd = srms
                override = True
            else:
                rd = dread
                override = False

            # count above sky
            diff = dtarg - slevel

            if extype == "normal" or extype == "optimal":
                if extype == "optimal":
                    # optimal extraction. Need the profile
                    warnings.warn(
                        "Transmission plot is not reliable with optimal extraction"
                    )

                    mbeta = store["mbeta"]
                    if mbeta > 0.0:
                        prof = fitting.moffat(
                            X[dok],
                            Y[dok],
                            0.0,
                            1.0,
                            0.0,
                            0.0,
                            mfwhm,
                            mbeta,
                            wdata.xbin,
                            wdata.ybin,
                            rfile["apertures"]["fit_ndiv"],
                        )
                    else:
                        prof = fitting.gaussian(
                            X[dok],
                            Y[dok],
                            0.0,
                            1.0,
                            0.0,
                            0.0,
                            mfwhm,
                            wdata.xbin,
                            wdata.ybin,
                            rfile["apertures"]["fit_ndiv"],
                        )

                    # multiply weights by the profile
                    wtarg *= prof

                # now extract
                counts = (wtarg * diff).sum()

                if override:
                    # in this case, the "readout noise" includes the component
                    # due to the sky background so we use the sky-subtracted
                    # counts above sky for the object contribution.
                    var = (wtarg**2 * (rd**2 + np.maximum(0, diff) / dgain)).sum()
                else:
                    # in this case we are using the true readout noise and we
                    # just use the data (which should be debiassed) without
                    # removal of the sky.
                    var = (wtarg**2 * (rd**2 + np.maximum(0, dbtarg) / dgain)).sum()

                if serror > 0:
                    # add in factor due to uncertainty in sky estimate
                    var += (wtarg.sum() * serror) ** 2

                countse = np.sqrt(var)

            else:
                raise hcam.HipercamError(
                    "extraction type = {:s} not recognised".format(extype)
                )

            info = store[apnam]

            results[apnam] = {
                "x": aper.x,
                "xe": info["xe"],
                "y": aper.y,
                "ye": info["ye"],
                "fwhm": info["fwhm"],
                "fwhme": info["fwhme"],
                "beta": info["beta"],
                "betae": info["betae"],
                "counts": counts,
                "countse": countse,
                "sky": slevel,
                "skye": serror,
                "nsky": nsky,
                "nrej": nrej,
                "flag": flag,
                "cmax": cmax,
            }

        except hcam.HipercamError as err:
            print(
                f"CCD {cnam}, reference aperture {apnam},"
                f" fit failed; extraction aborted. Error = {err}",
                file=sys.stderr,
            )
            info = store[apnam]
            flag |= hcam.NO_EXTRACTION

            results[apnam] = {
                "x": aper.x,
                "xe": info["xe"],
                "y": aper.y,
                "ye": info["ye"],
                "fwhm": info["fwhm"],
                "fwhme": info["fwhme"],
                "beta": info["beta"],
                "betae": info["betae"],
                "counts": NaN,
                "countse": NaN,
                "sky": NaN,
                "skye": NaN,
                "nsky": 0,
                "nrej": 0,
                "flag": flag,
                "cmax": 0,
            }

    # finally, we are done
    return results


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
                return "nan"
            else:
                return f"{value:{form}}"

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
                        cnam,
                        nframe,
                        mjd,
                        mjdok,
                        expose,
                        float2str(mfwhm, ".2f", -1.0),
                        float2str(mbeta, ".2f", -1.0),
                    )
                )

                # now for data per aperture
                for apnam in ccdaper:
                    r = reses[apnam]
                    self.log.write(
                        "{} {} {} {} {} {} {} {} {} {} {} {} "
                        "{:d} {:d} {:d} {:d} ".format(
                            float2str(r["x"], ".4f"),
                            float2str(r["xe"], ".4f"),
                            float2str(r["y"], ".4f"),
                            float2str(r["ye"], ".4f"),
                            float2str(r["fwhm"], ".3f"),
                            float2str(r["fwhme"], ".3f", -1.0),
                            float2str(r["beta"], ".3f"),
                            float2str(r["betae"], ".3f"),
                            float2str(r["counts"], ".1f"),
                            float2str(r["countse"], ".1f"),
                            float2str(r["sky"], ".2f"),
                            float2str(r["skye"], ".2f"),
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
        self.log.write("""#
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
""".format(hipercam_version=self.hipercam_version))

        # second, list the command-line inputs to the logfile
        for line in self.plist:
            self.log.write("# {:s}".format(line))

        # third, list the reduce file
        self.log.write("""#
# and here is a minimal version of the reduce file used ['rfile' above] with
# all between-line comments removed for compactness:
#
""")

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
        self.log.write("""#
# Next here is the aperture file used in JSON-style format that (without
# the initial comment hashes) is readable by setaper and reduce:
#
""")
        # convert aperture file to JSON-style string, split line by line,
        # pre-pend comment and indentation, write out to the logfile.
        lines = [
            "#   {:s}\n".format(line) for line in self.rfile.aper.toString().split("\n")
        ]
        for line in lines:
            self.log.write(line)

        # fifth the column names for each CCD which has any
        # apertures
        self.log.write("""#
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
""")
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
                cnames = (
                    "# {:s} = CCD nframe MJD-{:d} MJDok Exptim mfwhm mbeta ".format(
                        cnam, toffset
                    )
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
        self.log.write("""#
# End of column name definitions
#
# Now follow a similar series of datatypes which are designed to be used to
# build a numpy.dtype for each CCD to allow them to read into numpy structured
# arrays.
#
# Start of data type definitions:
#
""")
        atypes = "f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 i4 i4 i4 u4 "
        for cnam, ccdaper in self.rfile.aper.items():
            if len(ccdaper) == 0:
                # nothing will be written for CCDs without
                # apertures
                continue

            self.log.write(
                "# {:s} = s i4 f8 ? f4 f4 f4 {:s}\n".format(cnam, len(ccdaper) * atypes)
            )

        self.log.write("""#
# End of data type definitions
#
""")

    # that's it for the headers!
