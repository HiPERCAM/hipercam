import sys
import os
import tempfile

import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft

import hipercam as hcam
from hipercam import cline, utils, spooler, fringe
from hipercam.cline import Cline

__all__ = [
    "makefringe",
]

####################################################
#
# makefringe -- makes flat fields from a set of frames
#
####################################################


def makefringe(args=None):
    """``makefringe [source] (run first last [twait tmax] | flist) (bias
    flat dark) fpair ([nhalf]) ccd fwhm [clobber] output``

    Averages a set of images to make a frame for defringing (referred
    to elsewhere as a "fringe map").

    At long wavelengths, CCDs suffer what is known as "fringing",
    which in terms of structure looks something like the coloured
    patterns you see when there is a thin layer of oil on water. Both
    patterns are caused by interference.  In the case of CCDs this is
    caused by sky emission lines and is variable in strength, and
    additive in nature.

    De-fringing in |hiper| is implemented according to the method
    suggested by Snodgrass & Carry (2013Msngr.152...14S). The idea is
    to create a set of pairs of points marking the peak and troughs of
    fringes. The difference in intensity between these in the data to
    be corrected is compared with the difference in a reference
    "fringe map" by taking their ratio. This is done many times so
    that the median can be taken to eliminate ratios affected by
    celestial targets or cosmic rays. the ratios of which will be used
    to estimate the level of fringing compared to a reference
    frame. ``makefringe`` is used to make the reference fringe map. It
    works as follows: given an input list of files (or optionally a
    single run), it reads them all in, optionally debiases,
    dark-subtracts and flat-fields them, calculates a median count
    level of each one which is subtracted from each CCD
    individually. The pixel-by-pixel median of all frames is then
    calculated. It is also possible to apply the fringe pair ratio
    measurement to scale each frame before combining which allows
    frames of similar pattern but varying amplitude to be
    combined. Finally, the output can be smoothed because fringes are
    usually a medium to large scale pattern.

    Parameters:

        source : string [hidden]
           Data source, five options:

               | 'hs' : HiPERCAM server
               | 'hl' : local HiPERCAM FITS file
               | 'us' : ULTRACAM server
               | 'ul' : local ULTRACAM .xml/.dat files
               | 'hf' : list of HiPERCAM hcm FITS-format files

           'hf' is used to look at sets of frames generated by 'grab' or
           converted from foreign data formats.

        run : string [if source ends 's' or 'l']
           run number to access, e.g. 'run034'

        flist : string [if source ends 'f']
           name of file list. Assumed that bias, flat-fielding etc have been applied to
           these already.

        first : int [if source ends 's' or 'l']
           exposure number to start from. 1 = first frame ('0' is
           not supported).

        last : int [if source ends 's' or 'l']
           last exposure number must be >= first or 0 for the whole lot.

        twait : float [if source ends 's' or 'l'; hidden]
           time to wait between attempts to find a new exposure, seconds.

        tmax : float [if source ends 's' or 'l'; hidden]
           maximum time to wait between attempts to find a new exposure,
           seconds.

        bias : str [if source ends 's' or 'l']
           Name of bias frame to subtract, 'none' to ignore.

        flat : str [if source ends 's' or 'l']
           Name of flat field, 'none' to ignore.

        dark : str [if source ends 's' or 'l']
           Name of dark frame to subtract, 'none' to ignore. Note that
           it is assumed all CCDs have the same exposure time when making
           a dark correction.

        fpair : str
           FringePair file (see setfringe), or 'none' to ignore. If specified
           if will be used to scale the frames prior to combining.

        nhalf : int [if fpair is not 'none']
           When calculating the differences for fringe measurement,
           a region extending +/-nhalf binned pixels will be used when
           measuring the amplitudes. Basically helps the stats.

        ccd : str
           CCD(s) to process, '0' for all, '1 3' for '1' and '3' only, etc.

        fwhm : float
           FWHM (unbinned pixels) of gaussian smoothing to apply to the output.
           0 to ignore. Should be << smallest scale between fringes. Probably
           no more than 4.

        clobber : bool [hidden]
           clobber any pre-existing output files

        output  : str
           output file. Set by default to match the last part of "run"
           (but it will have a different extension so they won't
           clash)

      .. Note::

         This routine writes the files returned by 'grab' to
         automatically generated files, typically in .hipercam/tmp, to
         avoid polluting the working directory. These are removed at
         the end, but may not be if you ctrl-C. You should check
         .hipercam/tmp for redundant files every so often

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("source", Cline.GLOBAL, Cline.HIDE)
        cl.register("run", Cline.GLOBAL, Cline.PROMPT)
        cl.register("first", Cline.LOCAL, Cline.PROMPT)
        cl.register("last", Cline.LOCAL, Cline.PROMPT)
        cl.register("twait", Cline.LOCAL, Cline.HIDE)
        cl.register("tmax", Cline.LOCAL, Cline.HIDE)
        cl.register("flist", Cline.LOCAL, Cline.PROMPT)
        cl.register("bias", Cline.LOCAL, Cline.PROMPT)
        cl.register("flat", Cline.LOCAL, Cline.PROMPT)
        cl.register("dark", Cline.LOCAL, Cline.PROMPT)
        cl.register("fpair", Cline.LOCAL, Cline.PROMPT)
        cl.register("nhalf", Cline.LOCAL, Cline.HIDE)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("fwhm", Cline.LOCAL, Cline.PROMPT)
        cl.register("clobber", Cline.LOCAL, Cline.HIDE)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value(
            "source",
            "data source [hs, hl, us, ul, hf]",
            "hl",
            lvals=("hs", "hl", "us", "ul", "hf"),
        )

        # set a flag
        server_or_local = source.endswith("s") or source.endswith("l")

        if server_or_local:
            resource = cl.get_value("run", "run name", "run005")
            root = os.path.basename(resource)
            cl.set_default('output', cline.Fname(root, hcam.HCAM))
            first = cl.get_value("first", "first frame to average", 1, 1)
            last = cl.get_value("last", "last frame to average (0 for all)", first, 0)
            twait = cl.get_value(
                "twait", "time to wait for a new frame [secs]", 1.0, 0.0
            )
            tmax = cl.get_value(
                "tmax", "maximum time to wait for a new frame [secs]", 10.0, 0.0
            )

            # bias frame (if any)
            bias = cl.get_value(
                "bias",
                "bias frame ['none' to ignore]",
                cline.Fname("bias", hcam.HCAM),
                ignore="none",
            )

            # flat field (if any)
            flat = cl.get_value(
                "flat",
                "flat field frame ['none' to ignore]",
                cline.Fname("flat", hcam.HCAM),
                ignore="none",
            )

            # dark frame (if any)
            dark = cl.get_value(
                "dark",
                "dark frame ['none' to ignore]",
                cline.Fname("dark", hcam.HCAM),
                ignore="none",
            )

        else:
            resource = cl.get_value(
                "flist", "file list", cline.Fname("files.lis", hcam.LIST)
            )
            first = 1

        fpair = cl.get_value(
            "fpair", "fringe pair file",
            cline.Fname("fringe", hcam.FRNG),
            ignore="none"
        )
        if fpair is not None:
            nhalf = cl.get_value(
                "nhalf", "half-size of fringe measurement region",
                2, 0
            )

        ccdinf = spooler.get_ccd_pars(source, resource)

        if len(ccdinf) > 1:
            ccd = cl.get_value("ccd", "CCD(s) to process [0 for all]", "0")
            if ccd == "0":
                ccds = list(ccdinf.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(ccdinf.keys())

        fwhm = cl.get_value(
            "fwhm", "FWHM of gaussian smoothing to apply to final result [unbinned pixels, 0 to ignore]", 4., 0.
        )
        clobber = cl.get_value(
            "clobber", "clobber any pre-existing files on output", False
        )

        output = cl.get_value(
            "output",
            "median output fringe frame",
            cline.Fname(
                "hcam", hcam.HCAM, cline.Fname.NEW if clobber else cline.Fname.NOCLOBBER
            ),
        )

    # inputs done with.

    try:

        # big try / except section here to trap ctrl-C to allow the temporary
        # files to be deleted. First make a directory for the temporary files

        if server_or_local:
            print("\nCalling 'grab' ...")

            args = [
                None,
                "prompt",
                source,
                resource,
                "yes",
                str(first),
                str(last),
                "no",
                str(twait),
                str(tmax),
                "none" if bias is None else bias,
                "none" if dark is None else dark,
                "none" if flat is None else flat,
                "none", "f32",
            ]
            resource = hcam.scripts.grab(args)

        # at this point 'resource' is a list of files, no matter the input
        # method.

        # Read all the files to determine mean levels (after bias
        # subtraction) save the bias-subtracted, flat-fielded,
        # mean-level subtracted and normalised results to temporary
        # files
        print("Reading all files in to determine their mean levels")
        bframe, fframe, dframe, fpairobj = None, None, None, None
        medians = {}
        for cnam in ccds:
            medians[cnam] = {}

        # We might have a load of temporaries from grab, but we are about to
        # make some more to save the bias-subtracted normalised versions.
        tdir = utils.temp_dir()

        fnames = []

        mtype = []
        fref = {}
        with spooler.HcamListSpool(resource) as spool:

            for mccd in spool:

                if fpair is not None:
                    # prepare fringepairs
                    if fpairobj is None:
                        fpairobj = fringe.MccdFringePair.read(fpair)
                        fpairobj = fpairobj.crop(mccd, nhalf)

                # here we determine the median levels, subtract them
                # then (optionally) normalise by the ratio of fringe
                # amplitudes

                # generate the name to save to automatically
                fd, fname = tempfile.mkstemp(suffix=hcam.HCAM, dir=tdir)

                for cnam in ccds:
                    # its unlikely that fringe frames would be taken with skips, but
                    # you never know. Eliminate them from consideration now.
                    ccd = mccd[cnam]
                    if ccd.is_data():
                        cmedian = mccd[cnam].median()
                        medians[cnam][fname] = cmedian
                        mccd[cnam] -= cmedian
                        if cnam not in fref:
                            # save reference
                            fref[cnam] = mccd[cnam]

                        elif fpair is not None:
                            fscale = fpairobj[cnam].scale(
                                mccd[cnam], fref[cnam], nhalf
                            )
                            mccd[cnam] *= fscale
                            print(f'{cnam}, fscale = {fscale:.3f}')

                # write to disk, save the name, close the filehandle
                mccd.write(fname)
                fnames.append(fname)
                os.close(fd)

                # a bit of progress info
                print(f"Saved {', '.join(mtype)} frame to {fname}")

        # now we go through CCD by CCD, using the first as a template
        # for the window names in which we will also store the results.
        template = hcam.MCCD.read(fnames[0])

        # Now process each file CCD by CCD to reduce the memory
        # footprint
        for cnam in ccds:

            # Read files into memory, insisting that they
            # all have the same set of CCDs
            print(f"\nLoading all CCDs labelled '{cnam}'")

            accds = []
            with spooler.HcamListSpool(fnames, cnam) as spool:

                mean = None
                for ccd in spool:
                    if ccd.is_data():
                        accds.append(ccd)

            if len(accds) == 0:
                raise hcam.HipercamError(
                    f"Found no valid examples of CCD {cnam}"
                    f" in {fnames}"
                )
            else:
                print(f"Loaded {len(accds)} CCDs")

            # Median combine
            for wnam, wind in template[cnam].items():

                # build list of all data arrays
                arrs = [ccd[wnam].data for ccd in accds]

                # convert to 3D numpy array
                arr3d = np.stack(arrs)

                wind.data = np.median(arr3d, axis=0)

            # Add history
            template[cnam].head.add_history(
                f"Median combine of {len(accds)} images"
            )
            print("Computed and stored their pixel-by-pixel median")

        # Remove any CCDs not included to avoid impression of having done
        # something to them
        dcnams = []
        for cnam in template.keys():
            if cnam not in ccds:
                dcnams.append(cnam)
        for cnam in dcnams:
            del template[cnam]

        if fwhm > 0:
            # apply gaussian smoothing
            print(f'Smoothing with FWHM = {fwhm}')
            sigma = fwhm / np.sqrt(8 * np.log(2))
            kern = Gaussian2DKernel(sigma)
            for ccd in template.values():
                for wind in ccd.values():
                    wind.data = convolve_fft(wind.data, kern, "fill", wind.median())

        # write out
        template.write(output, clobber)
        print(f"\nFinal result written to {output}")

    except KeyboardInterrupt:
        print("\nmakefringe aborted")

    if server_or_local:
        # grab has created a load of temporaries, including the file
        # list 'resource'
        with open(resource) as fin:
            for fname in fin:
                os.remove(fname.strip())
        os.remove(resource)

    # and another load were created during the later running of the script
    for fname in fnames:
        os.remove(fname)

    print("temporary files have been removed.")
