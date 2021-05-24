"""Command line script to grab images"""

import sys
import os
import time
import tempfile

import numpy as np

import hipercam as hcam
from hipercam import cline, utils, spooler, fringe
from hipercam.cline import Cline

__all__ = [
    "grab",
]

###########################################################
#
# grab -- downloads a series of images from a raw data file
#
###########################################################


def grab(args=None):
    """``grab [source] run [temp] (ndigit) first last trim [twait tmax]
    ([ncol nrow]) bias dark flat fmap (fpair [nhalf rmin rmax verbose])
    [dtype]``

    This downloads a sequence of images from a raw data file and writes them
    out to a series CCD / MCCD files.

    Parameters:

       source : string [hidden]
           Data source, four options:

              | 'hs' : HiPERCAM server
              | 'hl' : local HiPERCAM FITS file
              | 'us' : ULTRACAM server
              | 'ul' : local ULTRACAM .xml/.dat files

           The standard start-off default for ``source`` can be set
           using the environment variable
           HIPERCAM_DEFAULT_SOURCE. e.g. in bash :code:`export
           HIPERCAM_DEFAULT_SOURCE="us"` would ensure it always
           started with the ULTRACAM server by default. If
           unspecified, it defaults to 'hl'.

       run : string
           run name to access

       temp : bool [hidden, defaults to False]
           True to indicate that the frames should be written to temporary
           files with automatically-generated names in the expectation of
           deleting them later. This also writes out a file listing the names.
           The aim of this is to allow grab to be used as a feeder for other
           routines in larger scripts.  If temp == True, grab will return with
           the name of the list of hcm files. Look at 'makebias' for an
           example that uses this feature.

       ndigit : int [if not temp]
           Files created will be written as 'run005_0013.fits' etc. `ndigit`
           is the number of digits used for the frame number (4 in this
           case). Any pre-existing files of the same name will be
           over-written.

       first : int
           First frame to access

       last : int
           Last frame to access, 0 for the lot

       trim : bool
           True to trim columns and/or rows off the edges of windows nearest
           the readout. Particularly useful for ULTRACAM windowed data where
           the first few rows and columns can contain bad data.

       ncol : int [if trim, hidden]
           Number of columns to remove (on left of left-hand window, and right
           of right-hand windows)

       nrow : int [if trim, hidden]
           Number of rows to remove (bottom of windows)

       twait : float [hidden]
           time to wait between attempts to find a new exposure, seconds.

       tmax : float [hidden]
           maximum time to wait between attempts to find a new exposure,
           seconds.

       bias : str
           Name of bias frame to subtract, 'none' to ignore.

       dark : str
           Name of dark frame, 'none' to ignore.

       flat : str
           Name of flat field to divide by, 'none' to ignore. Should normally
           only be used in conjunction with a bias, although it does allow you
           to specify a flat even if you haven't specified a bias.

       fmap : str
           Name of fringe map (see e.g. `makefringe`), 'none' to ignore.

       fpair : str [if fmap is not 'none']
           Name of fringe pair file (see e.g. `setfringe`). Required if
           a fringe map has been specified.

       nhalf : int [if fmap is not 'none', hidden]
           When calculating the differences for fringe measurement,
           a region extending +/-nhalf binned pixels will be used when
           measuring the amplitudes. Basically helps the stats.

       rmin : float [if fmap is not 'none', hidden]
           Minimum individual ratio to accept prior to calculating the overall
           median in order to reduce the effect of outliers. Although all ratios
           should be positive, you might want to set this a little below zero
           to allow for some statistical fluctuation.

       rmax : float [if fmap is not 'none', hidden]
           Maximum individual ratio to accept prior to calculating the overall
           median in order to reduce the effect of outliers. Probably typically
           < 1 if fringe map was created from longer exposure data.

       verbose : bool
           True to print lots of details of fringe ratios

       dtype : string [hidden, defaults to 'f32']
           Data type on output. Options:

            | 'f32' : output as 32-bit floats (default)

            | 'f64' : output as 64-bit floats.

            | 'u16' : output as 16-bit unsigned integers. A warning will be
                      issued if loss of precision occurs; an exception will
                      be raised if the data are outside the range 0 to 65535.

    .. Note::

       |grab| is used by several other scripts such as |averun| so take great
       care when changing anything to do with its input parameters.

    """

    command, args = utils.script_args(args)

    # get inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("source", Cline.GLOBAL, Cline.HIDE)
        cl.register("run", Cline.GLOBAL, Cline.PROMPT)
        cl.register("temp", Cline.GLOBAL, Cline.HIDE)
        cl.register("ndigit", Cline.LOCAL, Cline.PROMPT)
        cl.register("first", Cline.LOCAL, Cline.PROMPT)
        cl.register("last", Cline.LOCAL, Cline.PROMPT)
        cl.register("trim", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ncol", Cline.GLOBAL, Cline.HIDE)
        cl.register("nrow", Cline.GLOBAL, Cline.HIDE)
        cl.register("twait", Cline.LOCAL, Cline.HIDE)
        cl.register("tmax", Cline.LOCAL, Cline.HIDE)
        cl.register("bias", Cline.GLOBAL, Cline.PROMPT)
        cl.register("flat", Cline.GLOBAL, Cline.PROMPT)
        cl.register("dark", Cline.GLOBAL, Cline.PROMPT)
        cl.register("fmap", Cline.GLOBAL, Cline.PROMPT)
        cl.register("fpair", Cline.GLOBAL, Cline.PROMPT)
        cl.register("nhalf", Cline.GLOBAL, Cline.HIDE)
        cl.register("rmin", Cline.GLOBAL, Cline.HIDE)
        cl.register("rmax", Cline.GLOBAL, Cline.HIDE)
        cl.register("verbose", Cline.LOCAL, Cline.HIDE)
        cl.register("dtype", Cline.LOCAL, Cline.HIDE)

        # get inputs
        default_source = os.environ.get('HIPERCAM_DEFAULT_SOURCE','hl')
        source = cl.get_value(
            "source",
            "data source [hs, hl, us, ul]",
            default_source,
            lvals=("hs", "hl", "us", "ul"),
        )

        # OK, more inputs
        resource = cl.get_value("run", "run name", "run005")

        cl.set_default("temp", False)
        temp = cl.get_value(
            "temp", "save to temporary automatically-generated file names?", True
        )

        if not temp:
            ndigit = cl.get_value(
                "ndigit", "number of digits in frame identifier", 3, 0
            )

        first = cl.get_value("first", "first frame to grab", 1, 0)
        last = cl.get_value("last", "last frame to grab", 0)
        if last < first and last != 0:
            sys.stderr.write("last must be >= first or 0")
            sys.exit(1)

        trim = cl.get_value("trim", "do you want to trim edges of windows?", True)
        if trim:
            ncol = cl.get_value("ncol", "number of columns to trim from windows", 0)
            nrow = cl.get_value("nrow", "number of rows to trim from windows", 0)

        twait = cl.get_value("twait", "time to wait for a new frame [secs]", 1.0, 0.0)
        tmax = cl.get_value(
            "tmax", "maximum time to wait for a new frame [secs]", 10.0, 0.0
        )

        bias = cl.get_value(
            "bias",
            "bias frame ['none' to ignore]",
            cline.Fname("bias", hcam.HCAM),
            ignore="none",
        )
        if bias is not None:
            # read the bias frame
            bias = hcam.MCCD.read(bias)

        # dark
        dark = cl.get_value(
            "dark", "dark frame ['none' to ignore]",
            cline.Fname("dark", hcam.HCAM), ignore="none"
        )
        if dark is not None:
            # read the dark frame
            dark = hcam.MCCD.read(dark)

        # flat
        flat = cl.get_value(
            "flat", "flat field frame ['none' to ignore]",
            cline.Fname("flat", hcam.HCAM), ignore="none"
        )
        if flat is not None:
            # read the flat field frame
            flat = hcam.MCCD.read(flat)

        # fringe file (if any)
        fmap = cl.get_value(
            "fmap", "fringe map ['none' to ignore]",
            cline.Fname("fmap", hcam.HCAM), ignore="none",
        )
        if fmap is not None:
            # read the fringe frame and prompt other parameters
            fmap = hcam.MCCD.read(fmap)

            fpair = cl.get_value(
                "fpair", "fringe pair file",
                cline.Fname("fringe", hcam.FRNG)
            )
            fpair = fringe.MccdFringePair.read(fpair)

            nhalf = cl.get_value(
                "nhalf", "half-size of fringe measurement regions",
                2, 0
            )
            rmin = cl.get_value(
                "rmin", "minimum fringe pair ratio", -0.2
            )
            rmax = cl.get_value(
                "rmax", "maximum fringe pair ratio", 1.0
            )
            verbose = cl.get_value(
                "verbose", "verbose output", False
            )

        cl.set_default("dtype", "f32")
        dtype = cl.get_value(
            "dtype", "data type [f32, f64, u16]", "f32", lvals=("f32", "f64", "u16")
        )

    # Now the actual work.

    # strip off extensions
    if resource.endswith(hcam.HRAW):
        resource = resource[: resource.find(hcam.HRAW)]

    # initialisations
    total_time = 0  # time waiting for new frame
    nframe = first
    root = os.path.basename(resource)
    bframe = None

    # Finally, we can go
    if temp:
        fnames = []
        tdir = utils.temp_dir()

    with spooler.data_source(source, resource, first) as spool:

        try:

            for mccd in spool:

                # Handle the waiting game ...
                give_up, try_again, total_time = spooler.hang_about(
                    mccd, twait, tmax, total_time
                )

                if give_up:
                    print("grab stopped")
                    break
                elif try_again:
                    continue

                # Trim the frames: ULTRACAM windowed data has bad
                # columns and rows on the sides of windows closest to
                # the readout which can badly affect reduction. This
                # option strips them.
                if trim:
                    hcam.ccd.trim_ultracam(mccd, ncol, nrow)

                if nframe == first:
                    # First time through, need to manipulate calibration data
                    if bias is not None:
                        bias = bias.crop(mccd)
                        bexpose = bias.head.get("EXPTIME", 0.0)
                    else:
                        bexpose = 0.
                    if flat is not None:
                        flat = flat.crop(mccd)
                    if dark is not None:
                        dark = dark.crop(mccd)
                    if fmap is not None:
                        fmap = fmap.crop(mccd)
                        fpair = fpair.crop(mccd, nhalf)

                # now any time through, apply calibrations
                if bias is not None:
                    mccd -= bias

                if dark is not None:
                    # subtract dark, CCD by CCD
                    dexpose = dark.head["EXPTIME"]
                    for cnam in mccd:
                        ccd = mccd[cnam]
                        cexpose = ccd.head["EXPTIME"]
                        scale = (cexpose - bexpose) / dexpose
                        ccd -= scale * dark[cnam]

                if flat is not None:
                    mccd /= flat

                if fmap is not None:
                    # apply CCD by CCD
                    for cnam in fmap:
                        if cnam in fpair:
                            ccd = mccd[cnam]
                            if verbose:
                                print(f' CCD {cnam}')

                            fscale = fpair[cnam].scale(
                                ccd, fmap[cnam], nhalf, rmin, rmax,
                                verbose=verbose
                            )

                            if verbose:
                                print(f'  Median scale factor = {fscale}')

                            ccd -= fscale*fmap[cnam]

                if dtype == "u16":
                    mccd.uint16()
                elif dtype == "f32":
                    mccd.float32()
                elif dtype == "f64":
                    mccd.float64()

                # write to disk
                if temp:
                    # generate name automatically
                    fd, fname = tempfile.mkstemp(suffix=hcam.HCAM, dir=tdir)
                    mccd.write(fname, True)
                    os.close(fd)
                    fnames.append(fname)
                else:
                    fname = "{:s}_{:0{:d}}{:s}".format(root, nframe, ndigit, hcam.HCAM)
                    mccd.write(fname, True)

                print("Written frame {:d} to {:s}".format(nframe, fname))

                # update the frame number
                nframe += 1
                if last and nframe > last:
                    break

        except KeyboardInterrupt:
            # trap ctrl-C so we can delete temporary files if temp
            if temp:
                for fname in fnames:
                    os.remove(fname)
                print("\ntemporary files deleted")
                print("grab aborted")
            else:
                print("\ngrab aborted")
            sys.exit(1)

    if temp:
        if len(fnames) == 0:
            raise hcam.HipercamError(
                'no files were grabbed; please check input parameters, especially "first"'
            )
        # write the file names to a list
        fd, fname = tempfile.mkstemp(suffix=hcam.LIST, dir=tdir)
        with open(fname, "w") as fout:
            for fnam in fnames:
                fout.write(fnam + "\n")
        os.close(fd)
        print("temporary file names written to {:s}".format(fname))

        # return the name of the file list
        return fname
