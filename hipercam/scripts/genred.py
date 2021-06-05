import sys
import os
from time import gmtime, strftime

import hipercam as hcam
from hipercam import cline, utils, reduction
from hipercam.cline import Cline

# get hipercam version to write into the reduce file
from pkg_resources import get_distribution, DistributionNotFound

try:
    hipercam_version = get_distribution("hipercam").version
except DistributionNotFound:
    hipercam_version = "not found"

__all__ = [
    "genred",
]

################################################
#
# genred -- generates a reduce file
#
################################################


def genred(args=None):
    """``genred apfile rfile bias dark flat fmap fpair seeing (binfac)
    template (inst (nccd (ccd) nonlin sat scale readout gain))``

    Generates a reduce file as needed by |reduce| or |psf_reduce|. You
    give it the name of an aperture file, calibration frames and a few
    other parameters and it will write out a reduce file which you can
    then refine by hand. The aim is to try to get a reduce file that
    is self-consistent and works to get observing started as soon as
    possible. The parameters are not necessarily the best.

    It is assumed that the target is called '1', the main comparison '2'.

    To avoid endless prompts, genred does not prompt for all
    parameters but accepts a "template" file which can just be from a
    previous run of genred where such parameters can be altered by
    editing with a standard text editor before running genred. genred
    recognises some standard instrument telescope combinations which
    it uses to set parameters such as the pixel scale and readout
    noise, but if you choose other, you will be prompted for these
    details. The main parameters it does prompt for are calibration
    files, and it allows you to define a "seeing" which then controls
    the setting of various profile fit parameters. This can be ignored
    by setting = 0, in which case any template parameters will be passed
    unchanged.

    Parameters:

        apfile : str
           the input aperture file created using |setaper| (default
           extension .ape). This will be read for the targets. The
           main target will be assumed to have been called '1', the
           main comparison '2'. If there is a '3' it will be plotted
           relative to '2'; all others will be ignored for plotting
           purposes. Target '2' will be used to define the position
           and transmission plots for one CCD only [user
           definable]. Target '1' will be used for the seeing plot
           unless it is linked when target '2' will be used instead.

        rfile : str
           the output reduce file created by this script. The main
           target will be assumed to have been called '1', the main
           comparison '2'. If there is a '3' it will be plotted
           relative to '2'; all others will be ignored for plotting
           purposes.

        bias : str
           Name of bias frame; 'none' to ignore.

        dark : str
           Name of dark frame; 'none' to ignore.

        flat : str
           Name of flat field frame; 'none' to ignore.

        fmap : str
           Name of fringe frame; 'none' to ignore.

        fpair : str [if fmap not 'none']
           File defining pairs to measure fringe amplitudes.

        seeing : float
           estimate of seeing which will be used to define several of
           the profile fitting parameters. Enter 0 to ignore and use
           defaults from the template or genred instead.

        binfac : int [if seeing > 0]
           binning factor. e.g. 4 if using 4x4. Needed to optimise some
           profile parameters.

        template : str
           Reduce file to use as a template for any parameters not set
           by genred. This allows one to tweak settings by running
           genred without a template, then modifying the resultant
           file and using it as a template. 'none' to ignore, and then
           these values will be set by genred by default. Some will be
           modified anyway (e.g. the calibration file names), and if
           seeing > 0, then several of the profile fitting paramaters
           will be adapted to the value given along with the
           instrument parameters.  genred will try to keep as many of
           the 'extraction', 'position', 'transmission', 'light' etc
           lines associated with particular CCDs in this case without
           modification, and will pass readout, saturation levels etc
           unchanged from the template.

        inst : str [if template == 'none']
           the instrument. 'hipercam-gtc', 'ultracam-ntt', 'ultraspec-tnt',
           'ultracam-wht', 'ultracam-vlt',  or 'other'. Sets pixel scale and
           standard readnoise params.

        nccd : int [if inst == 'other']
           number of CCDs.

        ccd : int [if nccd > 1]
           the CCD to use for position plots

        nonlin : list(float) [if inst == 'other']
           values for warning of non linearity

        sat : list(float) [if inst == 'other']
           values for warning of saturation

        scale : float [if inst == 'other']
           image scale in arcsec/pixel

        readout : float | string [if inst == 'other']
           readout noise, RMS ADU. Can either be a single value or an
           hcm file.

        gain : float | string [if inst == 'other']
           gain, electrons per ADU. Can either be a single value or an
           hcm file.

    """

    command, args = utils.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("apfile", Cline.LOCAL, Cline.PROMPT)
        cl.register("rfile", Cline.GLOBAL, Cline.PROMPT)
        cl.register("comment", Cline.LOCAL, Cline.PROMPT)
        cl.register("bias", Cline.LOCAL, Cline.PROMPT)
        cl.register("dark", Cline.LOCAL, Cline.PROMPT)
        cl.register("flat", Cline.LOCAL, Cline.PROMPT)
        cl.register("fmap", Cline.LOCAL, Cline.PROMPT)
        cl.register("fpair", Cline.LOCAL, Cline.PROMPT)
        cl.register("seeing", Cline.LOCAL, Cline.PROMPT)
        cl.register("binfac", Cline.LOCAL, Cline.PROMPT)
        cl.register("inst", Cline.LOCAL, Cline.PROMPT)
        cl.register("nccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("nonlin", Cline.LOCAL, Cline.PROMPT)
        cl.register("sat", Cline.LOCAL, Cline.PROMPT)
        cl.register("scale", Cline.LOCAL, Cline.PROMPT)
        cl.register("template", Cline.LOCAL, Cline.PROMPT)

        # get inputs

        # the aperture file
        apfile = cl.get_value(
            "apfile", "aperture input file", cline.Fname("aper.ape", hcam.APER)
        )
        # Read the aperture file
        aper = hcam.MccdAper.read(apfile)

        # the reduce file
        root = os.path.splitext(os.path.basename(apfile))[0]
        cl.set_default("rfile", cline.Fname(root, hcam.RED, cline.Fname.NEW))
        rfile = cl.get_value(
            "rfile", "reduce output file",
            cline.Fname("reduce.red", hcam.RED, cline.Fname.NEW),
        )

        # calibrations

        # bias frame
        bias = cl.get_value(
            "bias",
            "bias frame ['none' to ignore]",
            cline.Fname("bias", hcam.HCAM),
            ignore="none",
        )

        # dark frame
        dark = cl.get_value(
            "dark",
            "dark field frame ['none' to ignore]",
            cline.Fname("dark", hcam.HCAM),
            ignore="none",
        )

        # flat field frame
        flat = cl.get_value(
            "flat",
            "flat field frame ['none' to ignore]",
            cline.Fname("flat", hcam.HCAM),
            ignore="none",
        )

        # fringe frame
        fmap = cl.get_value(
            "fmap",
            "fringe frame ['none' to ignore]",
            cline.Fname("fringe", hcam.HCAM),
            ignore="none",
        )
        if fmap is not None:
            fpair = cl.get_value(
                "fpair", "fringe pair file",
                cline.Fname("fpair", hcam.FRNG),
            )
        else:
            fpair = None

        seeing = cl.get_value(
            "seeing",
            "approximate seeing [arcsec, 0 to ignore]",
            1.0, 0.
        )
        if seeing > 0.:
            binfac = cl.get_value(
                "binfac", "binning factor",
                1, 1
            )

        # template
        template = cl.get_value(
            "template",
            "template reduce file ['none' to ignore]",
            cline.Fname("template", hcam.RED),
            ignore="none",
        )

        if template is None:

            instrument = cl.get_value(
                "inst",
                "the instrument-telescope",
                "hipercam-gtc",
                lvals=[
                    "hipercam-gtc", "ultracam-ntt",
                    "ultracam-wht", "ultracam-vlt",
                    "ultraspec-tnt", "other"
                ],
            )

            if instrument.startswith("hipercam"):
                warn_levels = """
warn = 1 50000 64000
warn = 2 50000 64000
warn = 3 50000 64000
warn = 4 50000 64000
warn = 5 50000 64000
                """
                maxcpu = 5
                skipbadt = False
                scale = 0.081
                readout = '4.2'
                gain = '1.1'
                nccd = 5
                ccd = '2'

            elif instrument.startswith("ultracam"):
                warn_levels = """
warn = 1 28000 64000
warn = 2 28000 64000
warn = 3 50000 64000
            """
                maxcpu = 3
                skipbadt = True
                readout = '2.8'
                gain = '1.1'
                nccd = 3
                ccd = '2'
                if instrument.endswith("ntt"):
                    scale = 0.35
                elif instrument.endswith("wht"):
                    scale = 0.30
                elif instrument.endswith("vlt"):
                    scale = 0.15

            elif instrument.startswith("ultraspec"):
                warn_levels = """
warn = 1 60000 64000
            """
                maxcpu = 1
                skipbadt = True
                scale = 0.45
                readout = '5.0'
                gain = '1.1'
                nccd = 1
                ccd = '1'

            elif instrument == "other":

                # unrecognised instrument need extra stuff
                nccd = cl.get_value("nccd", "how many CCDs?", 3, 1)
                if nccd > 1:
                    ccd = cl.get_value("ccd", "which CCD for the position plot?", 1, 1, nccd)
                    ccd = str(ccd)
                else:
                    ccd = '1'

                # non-linearity warning levels
                nonlin = cl.get_default("nonlin")
                if nonlin is not None and len(nonlin) != nccd:
                    cl.set_default("nonlin", nccd * (50000.,))
                nonlins = cl.get_value(
                    "nonlin",
                    "levels to indicate non-linear behaviour, 1 per CCD",
                    nccd*(50000.,)
                )

                # saturation warning levels
                sat = cl.get_default("sat")
                if sat is not None and len(sat) != nccd:
                    cl.set_default("sat", nccd * (62000.,))
                sats = cl.get_value(
                    "sat",
                    "levels to indicate saturation, 1 per CCD",
                    nccd*(62000.,)
                )
                warn_levels = ""
                for n, (nonlin,sat) in zip(nonlins,sats):
                    warn_levels += \
                        f"warn = {n+1} {nonlin} {sat}\n"
                maxcpu = nccd

                scale = cl.get_value(
                    "scale", "image scale [arcsec/unbinned pixel]",
                    0.3, 0.001
                )

                readout = cl.get_value(
                    "readout", "readout noise, RMS ADU (float or file name)", "4.5"
                )

                gain = cl.get_value(
                    "gain", "gain, electrons/ADU, (float or file name)", "1.1"
                )

            else:
                raise hcam.HipercamError(
                    "Programming error: instrument unrecognised"
                )


            if ccd not in aper:
                # make sure 'ccd' is in aper, even if our
                # favoured one isn't
                for cnam in aper:
                    ccd = cnam
                    break

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff


    if template is not None:

        # read the template to define many unprompted values
        # shortcut names defines
        vals = reduction.Rfile.read(template, False)
        gsec = vals['general']
        instrument = gsec.get('instrument')
        scale = gsec.get('scale')
        warn_levels = ''.join(gsec['warn'])

        asec = vals['apertures']
        psfsec = vals['psf_photom']
        sksec = vals['sky']
        csec = vals['calibration']
        lcsec = vals['lcplot']
        lsec = vals['light']
        psec = vals.get('position',None)
        tsec = vals.get('transmission',None)
        ssec = vals.get('seeing',None)

    else:

        # define default values for unprompted params.
        if seeing == 0.:
            seeing = 1.5
            binfac = 1

        # same name as read into template
        vals = {}

        # general section
        gsec = vals['general'] = {}
        gsec['instrument'] = instrument
        gsec['scale'] = scale
        gsec['ldevice'] = '1/xs'
        gsec['lwidth'] = 0.
        gsec['lheight'] = 0.
        gsec['idevice'] = '2/xs'
        gsec['iwidth'] = 0.
        gsec['iheight'] = 0.
        gsec['toffset'] = 0
        gsec['skipbadt'] = 'yes'
        gsec['ncpu'] = 1
        gsec['ngroup'] = 1

        # aperture section
        asec = vals['apertures'] = {}
        asec['location'] = 'variable'
        asec['search_smooth_fft'] = 'no'
        asec['fit_method'] = 'moffat'
        asec['fit_beta'] = 5.0
        asec['fit_beta_max'] = 20.0
        asec['fit_ndiv'] = 0
        asec['fit_fwhm_fixed'] = 'no'
        asec['fit_thresh'] = 8.
        asec['fit_height_min_ref'] = 40.
        asec['fit_height_min_nrf'] = 10.
        asec['fit_alpha'] = 0.1

        # sky
        sksec = vals['sky'] = {}
        sksec['error'] = 'variance'
        sksec['method'] = 'clipped'
        sksec['thresh'] = 3.0

        # calibration
        csec = vals['calibration'] = {}
        csec['crop'] = 'yes'
        csec['nhalf'] = 3
        csec['rmin'] = -1.
        csec['rmax'] = 2.
        csec['readout'] = readout
        csec['gain'] = gain

        # PSF photom
        psfsec = vals['psf_photom'] = {}
        psfsec['gfac'] = 3.0
        psfsec['fit_half_width'] = 15.0
        psfsec['positions'] = 'fixed'

        # light curve
        lcsec = vals['lcplot'] = {}
        lcsec['xrange'] = 0
        lcsec['extend_x'] = 10.0

        # light
        lsec = vals['light'] = {}
        lsec['linear'] = 'no'
        lsec['y_fixed'] = 'no'
        lsec['y1'] = 0
        lsec['y2'] = 0
        lsec['extend_y'] = 0.1

        # position
        psec = vals['position'] = {}
        psec['height'] = 0.5
        psec['x_fixed'] = 'no'
        psec['x_min'] = -5.0
        psec['x_max'] = +5.0
        psec['y_fixed'] = 'no'
        psec['y_min'] = -5.0
        psec['y_max'] = +5.0
        psec['extend_y'] = 0.2

        # transmission
        tsec = vals['transmission'] = {}
        tsec['height'] = 0.5
        tsec['ymax'] = 110.

        # seeing
        ssec = vals['seeing'] = {}
        ssec['height'] = 0.5
        ssec['ymax'] = 2.999
        ssec['y_fixed'] = 'no'
        ssec['extend_y'] = 0.2

        fsec = vals['focal_mask'] = {}
        fsec['demask'] = 'no'
        fsec['dthresh'] = 3.0

    # apply seeing/instrument-related fixes...
    # *always* gets run in the no template case
    if seeing > 0.:
        asec['search_half_width'] = max(5., 4*seeing/scale)
        asec['search_smooth_fwhm'] = max(2.,seeing/scale/binfac)
        asec['fwhm_min'] = 2*binfac
        asec['fwhm'] = max(2*binfac, seeing/scale)
        asec['fwhm_max'] = 5.*seeing/scale
        asec['fit_max_shift'] = seeing/scale/3.
        asec['fit_half_width'] = max(5., 3*seeing/scale)
        asec['fit_diff'] = max(1., seeing/scale/4.)

        rfac = 1.8
        ramin = 3*binfac
        ramax = max(ramin,2.5*rfac*seeing/scale)
        sinner = ramax
        souter = 1.5*sinner

    # standard colours for CCDs
    if instrument.startswith("hipercam"):
        CCD_COLS = {
            "1": "blue",
            "2": "green",
            "3": "orange",
            "4": "red",
            "5": "darkred",
        }

    elif instrument.startswith("ultracam"):
        CCD_COLS = {
            "1": "red",
            "2": "green",
            "3": "blue",
        }

    elif instrument.startswith("ultraspec"):
        CCD_COLS = {
            "1": "green",
        }

    elif instrument == "other":
        # 'other'
        CCD_COLS = {
            "1": "green",
            "2": "blue",
            "3": "red",
            "4": "orange",
            "5": "purple",
        }
    else:
        raise hcam.HipercamError(
            "Programming error: instrument unrecognised"
        )

    # Extraction lines
    extraction = ""
    if template is None:
        # generate one for each CCD
        for cnam in aper:
            extraction += (
                f"{cnam} = variable normal"
                f" {rfac:.2f} {ramin:.1f} {ramax:.1f}"
                f" 2.5 {sinner:.1f} {sinner:.1f}"
                f" 3.0 {souter:.1f} {souter:.1f}\n"
            )
    else:
        # pass through unchanged where possible if
        # a template being used but check for consistency
        esec = template['extraction']
        for cnam, lst in esec.items():
            if cnam in aper:
                # unpack items
                ltype,etype,rfac,ramin,ramax,sifac,simin,simax,sofac,somin,somax = lst
                extraction += (
                    f"{cnam} = variable normal"
                    f" {rfac:.2f} {ramin:.1f} {ramax:.1f}"
                    f" {sifac:.2f} {simin:.1f} {simax:.1f}"
                    f" {sofac:.2f} {somin:.1f} {somax:.1f}\n"
                )
            else:
                warnings.warn(
                    f'CCD {nam} has an extraction line in {template} but no apertures in {apfile} and will be skipped'
                )

        # warn if there are apertures for a CCD but no extraction line
        for cnam in aper:
            if cnam not in esec:
                warnings.warn(
                    f'CCD {nam} has apertures defined in {apfile} but no extraction line in {template}'
                )

    # Generate the light curve plot lines
    light_plot = ""
    no_light = True

    if template is None:
        # Here we assume target in aperture 1, main comparison in aperture 2
        for cnam in aper:
            ccdaper = aper[cnam]
            if "1" in ccdaper and "2" in ccdaper:
                # favour 2 as the comparison
                light_plot += (
                    f'plot = {cnam} 1 2 0 1 {CCD_COLS.get(cnam,"black")} !  '
                    " # ccd, targ, comp, off, fac, dcol, ecol\n"
                )
                no_light = False

            elif "1" in ccdaper and "3" in ccdaper:
                # but try 3 if need be
                light_plot += (
                    f'plot = {cnam} 1 3 0 1 {CCD_COLS.get(cnam,"black")} !  '
                    " # ccd, targ, comp, off, fac, dcol, ecol\n"
                )
                no_light = False

            elif "1":
                # just plot 1 as last resort
                light_plot += (
                    f'plot = {cnam} 1 ! 0 1 {CCD_COLS.get(cnam,"black")} !  '
                    " # ccd, targ, comp, off, fac, dcol, ecol\n"
                )
                no_light = False

    else:

        # template case: pass as many of the lines as we can but
        # run checks of apertures
        plots = template['light']['plot']
        for pl in plots:
            cnam = pl['cnam']
            targ = pl['targ']
            comp = pl['comp']
            if cnam in aper and targ in aper[cnam] and (comp == '!' or comp in aper[cnam]):
                light_plot += (
                    f"plot = {cnam} {targ} {comp} {pl['off']} {pl['fac']} {pl['dcol']} {pl['ecol']}"
                    " # ccd, targ, comp, off, fac, dcol, ecol\n"
                )
                no_light = False
            else:
                warnings.warn(
                    "Light curve plot line:\n"
                    f"plot = {cnam} {targ} {comp} {pl['off']} {pl['fac']} {pl['dcol']} {pl['ecol']}\n"
                    f"from {template} is incompatible with {apfile} and has been skipped"
                )

    if no_light:
        raise hcam.HipercamError(
            "Found no targets for light curve plots in any CCD; cannot make light "
            "curve plot. Needs at least aperture '1' defined and preferably '2' as well,"
            "for at least one CCD"
        )

    # Generate the position plot lines. Favour the comparisons first, then the target.
    position_plot = ""
    no_position = True
    if template is None:
        ccdaper = aper[ccd]
        for targ in ('2','3','1'):
            if targ in ccdaper:
                position_plot += (
                    f'plot = {ccd} {targ} '
                    f'{CCD_COLS.get(cnam,"black")} ! '
                    '# ccd, targ, dcol, ecol\n'
                )
                no_position = False
                break

        if no_position:
            warnings.warn(
                f"Targets 2, 1, or 3 not found in CCD = {ccd} in {apfile}"
            )

    elif psec is not None:
        plots = psec['plot']
        for pl in plots:
            cnam, targ, dcol, ecol = pl
            if cnam in aper and targ in aper[cnam]:
                position_plot += (
                    f'plot = {cnam} {apnam} {dcol} {ecol} # ccd, targ, dcol, ecol\n'
                )
                no_position = False
            else:
                warnings.warn(
                    f"Position target {targ} / CCD {cnam} not found in {apfile}; will skip"
                )

    if no_position:
        warnings.warn(f"No position plot will be made")

    # Generate the transmission plot lines, again comparisons first,
    # target only if desperate.
    transmission_plot = ""
    no_transmission = True
    if template is None:
        for cnam in aper:
            ccdaper = aper[cnam]
            for targ in ('2','3','1'):
                if targ in ccdaper:
                    transmission_plot += (
                        f'plot = {cnam} {targ} '
                        f'{CCD_COLS.get(cnam,"black")} ! '
                        f"# ccd, targ, dcol, ecol\n"
                    )
                    no_transmission = False
                    break

        if no_transmission:
            warning.warn(
                f"Targets 2, 3, or 1 not found in any CCD within {apfile}; no transmission plot"
            )

    elif tsec is not None:
        plots = tsec['plot']
        for pl in plots:
            cnam, targ, dcol, ecol = pl
            if cnam in aper and targ in aper[cnam]:
                transmission_plot += (
                    f'plot = {cnam} {targ} {dcol} {ecol} # ccd, targ, dcol, ecol\n'
                )
                no_transmission = False
            else:
                warnings.warn(
                    f"Transmission target {targ} / CCD {cnam} not found in {apfile}; will skip"
                )

    if no_transmission:
        warnings.warn(f"No transmission plot will be made")

    # Generate the seeing plot lines. Prioritise the target in this case.
    seeing_plot = ""
    no_seeing = True
    if template is None:
        for cnam in aper:
            ccdaper = aper[cnam]
            for targ in ('1','2','3'):
                if targ in ccdaper and not ccdaper[targ].linked:
                    seeing_plot += (
                        f"plot = {cnam} {targ}"
                        f' {CCD_COLS.get(cnam,"black")} ! '
                        "# ccd, targ, dcol, ecol\n"
                    )
                    no_seeing = False
                    break

        if no_seeing:
            raise hcam.HipercamError(
                "Targets 1, 2 and 3 not found in any CCD"
                " (or they are linked); cannot make seeing plot"
            )

    elif ssec is not None:
        plots = tsec['plot']
        for pl in plots:
            cnam, targ, dcol, ecol = pl
            if cnam in aper and targ in aper[cnam] and not aper[cnam][targ].linked:
                seeing_plot += (
                    f'plot = {cnam} {targ} {dcol} {ecol} # ccd, targ, dcol, ecol\n'
                )
                no_seeing = False
            else:
                warnings.warn(
                    f"Seeing target {targ} / CCD {cnam} not found or linked in {apfile}; will skip"
                )

    if no_seeing:
        warnings.warn(f"No seeing plot will be made")

    # monitor targets (all of them by default)
    if template is None:
        targs = set()
        for cnam in aper:
            ccdaper = aper[cnam]
            for targ in ccdaper:
                targs.add(targ)
        monitor = ""
        for targ in sorted(targs):
            monitor += (
                f"{targ} = NO_EXTRACTION TARGET_SATURATED TARGET_AT_EDGE"
                " TARGET_NONLINEAR NO_SKY NO_FWHM NO_DATA SKY_AT_EDGE\n"
            )
    else:
        monitor = ""
        msec = template['monitor']
        rlook = dict([(k,v) for v,k in hcam.FLAGS])
        for targ, masks in msec:
            smasks = ' '.join([rlook[mask] for mask in masks])
            monitor += f'{targ} = {smasks}'

    # time stamp
    tstamp = strftime("%d %b %Y %H:%M:%S (UTC)", gmtime())

    # Finally, write out the reduce file ...
    with open(rfile, "w") as fout:
        # write out file
        fout.write(
            f"""#
# This is a HiPERCAM "reduce file" which defines the operation of the
# reduce script. It was written by the HiPERCAM pipeline command
# 'genred'.  It consists of a series of sections each of which
# contains a number of parameters. The file is self-documenting on the
# meaning of these parameters. The idea is that these are to large
# extent unchanging and it would be annoying to be prompted every time
# for them, but it also acts as a record of how reduction was carried
# out and is fed into the log file produce by 'reduce'.

# File written on {tstamp}
#
# HiPERCAM pipeline version: {hcam.REDUCE_FILE_VERSION}
#
# Start with some general items that tend not to change
# much. 'version' is the version of the reduce file format which
# changes more slowly than the software does. It must match the same
# parameter in 'reduce' for reduction to proceed. This is
# automatically the case at the time of creation, but old versions of
# the reduce file may become incompatible with later versions of
# reduce. Either they will require updating to be used, or the
# software version can be rolled back to give a compatible version of
# reduce using 'git'. The script 'rupdate', which attempts automatic
# update, may be worth trying if you need to update. It attempts to
# make the minimum changes needed to an old reduce file to run with
# later version dates.

[general]
version = {hcam.REDUCE_FILE_VERSION} # must be compatible with the version in reduce

instrument = {gsec['instrument']} # instrument
scale = {gsec['scale']} # pixel scale, arcsec/pixel

ldevice = {gsec['ldevice']} # PGPLOT plot device for light curve plots
lwidth = {gsec['lwidth']} # light curve plot width, inches, 0 to let program choose
lheight = {gsec['lheight']} # light curve plot height, inches

idevice = {gsec['idevice']} # PGPLOT plot device for image plots [if implot True]
iwidth = {gsec['iwidth']} # image curve plot width, inches, 0 to let program choose
iheight = {gsec['iheight']} # image curve plot height, inches

toffset = {gsec['toffset']} # integer offset to subtract from the MJD

# skip points with bad times in plots. HiPERCAM has a problem in not
# correctly indicating bad times so one does not usually want to
# skip "bad time" points, whereas one should for ULTRACAM and ULTRASPEC.
skipbadt = {gsec['skipbadt']}

# series of count levels at which warnings will be triggered for (a)
# non linearity and (b) saturation. Each line starts 'warn =', and is
# then followed by the CCD label, the non-linearity level and the
# saturation level

{warn_levels}

# The aperture reposition and extraction stages can be run in separate
# CPUs in parallel for each CCD potentially offering speed
# advantages. 'ncpu' is the number of CPUs to use for this. The
# maximum useful and best number to use is the number of CCDs in the
# instrument, e.g. 5 for HiPERCAM. You probably also want to leave at
# least one CPU to do other stuff, but if you have more than 2 CPUs,
# this parameter may help speed things. If you use this option (ncpu >
# 1), then there is also an advantage in terms of reducing
# parallelisation overheads in reading frames a few at a time before
# processing. This is controlled using 'ngroup'. i.e. with ngroup=10,
# 10 full frames are read before being processed. This parameter is
# ignored if ncpu==1

# I have had mixed results with this, which could depend a lot on the
# installation. On my home desktop is works better with ncpu=1 because
# it already seems to parallelise and it only gets worse if I split things
# up. But no harm trying.

ncpu = {gsec['ncpu']}
ngroup = {gsec['ngroup']}

# The next section '[apertures]' defines how the apertures are
# re-positioned from frame to frame. Apertures are re-positioned
# through a combination of a search near a start location followed by
# a 2D profile fit. Several parameters below are associated with this
# process and setting these right can be the key to a successful
# reduction.  If there are reference apertures, they are located first
# to give a mean shift. This is used to avoid the initial search for
# any non-reference apertures which has the advantage of reducing the
# chance of problems. The search is carried out by first extracting a
# sub-window centred on the last good position of a target. This is
# then smoothed by a gaussian (width 'search_smooth_fwhm'), and the
# closest peak to the last valid position higher than
# 'fit_height_min_ref' above background (median of the square box) is
# taken as the initial position for later profile fits. The smoothing
# serves to make the process more robust against cosmic rays. The
# width of the search box ('search_half_width') depends on how good
# the telescope guiding is. It should be large enough to cope with the
# largest likely shift in position between any two consecutive
# frames. Well-chosen reference targets, which should be isolated and
# bright, can help this process a great deal. The threshold is applied
# to the *smoothed* image. This means that it can be significantly
# lower than simply the raw peak height. e.g. a target might have a
# typical peak height around 100, in seeing of 4 pixels FWHM. If you
# smooth by 10 pixels, the peak height will drop to
# 100*4**2/(4**2+10**2) = 14 counts. It will be much more stable as a
# result, but you should then probably choose a threshold of 7 when
# you might have thought 50 was appropriate. The smoothing itself can
# be carried out by direct convolution or by an FFT-based method. The
# end-result is the same either way but for large values of
# 'search_smooth_fwhm', i.e. >> 1, FFTs may offer an advantage
# speed-wise. But the only way to tell is by explicity running with
# 'search_smooth_fft' switched from 'no' to 'yes'.

# The boxes for the fits ('fit_half_width') need to be large enough to
# include the target and a bit of sky to ensure that the FWHM is
# accurately measured, remembering that seeing can flare of course. If
# your target was defocussed, a gaussian or Moffat function will be a
# poor fit and you may be better keeping the FWHM fixed at a large
# value comparable to the widths of your defoccused images (and use
# the gaussian option in such cases). If the apertures are chosen to
# be fixed, there will be no search or fit carried out in which case
# you must choose 'fixed' as well when it comes the extraction since
# otherwise it needs a FWHM. 'fixed' is a last resort and you will
# very likely need to use large aperture radii in the extraction
# section.

# An obscure parameter is 'fit_ndiv'. If this is made > 0, the fit
# routine attempts to allow for pixellation by evaluating the profile
# at multiple points within each pixel of the fit. First it will
# evaluate the profile for every unbinned pixel within a binned pixel
# if the pixels are binned; second, it will evaluate the profile over
# an ndiv by ndiv square grid within each unbinned pixel. Thus ndiv=1
# means one evaluation at the centre of each unbinned pixel,
# etc. Obviously this will slow things, but it could help if your
# images are under-sampled. I would always start with fit_ndiv=0, and
# only raise it if the measured FWHM seem to be close to or below two
# binned pixels.

# If you use reference targets (you should if possible), the initial
# positions for the non-reference targets should be good. You can then
# guard further against problems using the parameter 'fit_max_shift'
# to reject positions for the non-reference targets that shift too far
# from the initial guess. 'fit_alpha' is another parameter that
# applies only in this case. If reference apertures are being used,
# the expected locations of non-reference apertures can be predicted
# with some confidence. In this case when the non-reference aperture's
# position is measured, its position will be adjusted by 'fit_alpha'
# times the measured change in its position. Its value is bounded by 0
# < fit_alpha <= 1. "1" just means use the full measured change from
# the current frame to update the position. Anything < 1 builds in a
# bit of past history. The hope is that this could make the aperture
# positioning, especially for faint targets, more robust to cosmic
# rays and other issues.  Of course it will correlate the positions
# from frame to frame. fit_alpha = 0.1 for instance will lead to a
# correlation length ~ 10 frames.

# If you use > 1 reference targets, then the parameter 'fit_diff'
# comes into play.  Multiple reference targets should move together
# and give very consistent shifts. If they don't, then a problem may
# have occurred, e.g. one or more have been affected by a meteor trail
# for instance. The maximum acceptable differential shift is defined
# by 'fit_diff'. If exceeded, then the entire extraction will be
# aborted and positions held fixed.

# To get and idea of the right values of some of these parameters, in
# particular the 'search_half_width', the height thresholds,
# 'fit_max_shift' and 'fit_diff', the easiest approach is probably to
# run a reduction with loose values and see how it goes.

[apertures]
aperfile = {apfile} # file of software apertures for each CCD
location = {asec['location']} # aperture locations: 'fixed' or 'variable'

search_half_width = {asec['search_half_width']:.1f} # for initial search for objects around previous position, unbinned pixels
search_smooth_fwhm = {asec['search_smooth_fwhm']:.1f} # smoothing FWHM, binned pixels
search_smooth_fft = {asec['search_smooth_fft']} # use FFTs for smoothing, 'yes' or 'no'.

fit_method = {asec['fit_method']} # gaussian or moffat
fit_beta = {asec['fit_beta']:.1f} # Moffat exponent
fit_beta_max = {asec['fit_beta_max']:.1f} # max Moffat expt
fit_fwhm = {asec['fwhm']:.1f} # FWHM, unbinned pixels
fit_fwhm_min = {asec['fwhm_min']:.1f} # Minimum FWHM, unbinned pixels [MUST be < fit_fwhm!]
fit_fwhm_max = {asec['fwhm_max']:.1f} # Maximum FWHM, unbinned pixels
fit_ndiv = {asec['fit_ndiv']} # sub-pixellation factor
fit_fwhm_fixed = {asec['fit_fwhm_fixed']} # Might want to set = 'yes' for defocussed images
fit_half_width = {asec['fit_half_width']:.1f} # for fit, unbinned pixels
fit_thresh = {asec['fit_thresh']:.2f} # RMS rejection threshold for fits
fit_height_min_ref = {asec['fit_height_min_ref']:.1f} # minimum height to accept a fit, reference aperture
fit_height_min_nrf = {asec['fit_height_min_nrf']:.1f} # minimum height to accept a fit, non-reference aperture
fit_max_shift = {asec['fit_max_shift']:.1f} # max. non-ref. shift, unbinned pixels.
fit_alpha = {asec['fit_alpha']:.2f} # Fraction of non-reference aperture shift to apply
fit_diff = {asec['fit_diff']:.2f} # Maximum differential shift of multiple reference apertures, unbinned

# The next lines define how the apertures will be re-sized and how the
# flux will be extracted from the aperture. There is one line per CCD
# with format:

# <CCD label> = <resize> <extract method> [scale min max] [scale min max]
#               [scale min max]

# where: <CCD label> is the CCD label; <resize> is either 'variable'
# or 'fixed' and sets how the aperture size will be determined -- if
# variable it will be scaled relative to the FWHM, so profile fitting
# will be attempted; <extract method> is either 'normal' or 'optimal'
# to say how the flux will be extracted -- 'normal' means a straight
# sum of sky subtracted flux over the aperture, 'optimal' use Tim
# Naylor's profile weighting, and requires profile fits to
# work. Finally there follow a series of numbers in three triplets,
# each of which is a scale factor relative to the FWHM for the
# aperture radius if the 'variable' option was chosen, then a minimum
# and a maximum aperture radius in unbinned pixels.  The three triples
# correspond to the target aperture radius, the inner sky radius and
# finally the outer sky radius. The mininum and maximum also apply if
# you choose 'fixed' apertures and can be used to override whatever
# value comes from the aperture file. A common approach is set them
# equal to each other to give a fixed value, especially for the sky
# where one does not necessarily want the radii to vary.  For PSF
# photometry, all these settings have no effect, but this section can
# still be used to determine which CCDs have fluxes extracted.

[extraction]
{extraction}

# The next lines are specific to the PSF photometry option. 'gfac' is
# used to label the sources according to groups, such that stars
# closer than 'gfac' times the FWHM are labelled in the same
# group. Each group has a PSF model fit independently. The reason
# behind the construction of groups is to reduce the dimensionality of
# the fitting procedure. Usually you want closely seperated stars to
# be fit simultaneously, but too large a value will mean fitting a
# model with many free parameters, which can fail to converge. The
# size of the box over which data is collected for fitting is set by
# 'fit_half_width'. Finally, 'positions' determines whether the star's
# positions should be considered variable in the PSF fitting. If this
# is set to fixed, the positions are held at the locations found in
# the aperture repositioning step, otherwise the positions are refined
# during PSF fitting. This step can fail for PSF photometry of faint
# sources.

[psf_photom]
gfac = {psfsec['gfac']:.1f}  # multiple of the FWHM to use in grouping objects
fit_half_width = {psfsec['fit_half_width']:.1f}  # size of window used to collect the data to do the fitting
positions = {psfsec['positions']}   # 'fixed' or 'variable'

# Next lines determine how the sky background level is
# calculated. Note you can only set error = variance if method =
# 'clipped'. 'median' should usually be avoided as it can cause
# noticable steps in light curves. It's here as a comparator.

[sky]
method = {sksec['method']} # 'clipped' | 'median'
error  = {sksec['error']} # 'variance' | 'photon': first uses actual variance of sky
thresh = {sksec['thresh']:.1f} # threshold in terms of RMS for 'clipped'

# Calibration frames and constants

# If you specify "!" for the readout, an attempt to estimate it from
# +/- 1 sigma percentiles will be made. This could help if you have no
# bias (and hence variance calculation from the count level will be
# wrong)

# The fringe parameters are the most complex. To de-fringe you need
# first a fringe map, e.g. as returned by makefringe. Then you need a
# file of FringePairs as returned by 'setfringe'. These are pairs of
# points placed at peaks and troughs of fringes. They will be used to
# measure a series of differences in each frame of interest and the
# fringe map from which a series of ratios is formed. The differences
# are measured from a group of pixels extending +/-nhlaf around the
# central pixel of each point of a FringePair. Thus nhalf=2 would give
# a 5x5 array (unbinned pixels). Finally to reduce the effect of bad
# values the individual ratios can be pruned if they lie outside a
# range rmin to rmax prior to their overall median being calculated.
# rmin should be small, but probably slightly less than zero to allow
# for a bit of statistical fluctuation, depending on your setup. One
# would normally expect rmax < 1 if the fringe map came from frames with
# a long exposure compared to the data.

[calibration]
crop = {csec['crop']} # Crop calibrations to match the data
bias = {'' if bias is None else bias} # Bias frame, blank to ignore
flat = {'' if flat is None else flat} # Flat field frame, blank to ignore
dark = {'' if dark is None else dark} # Dark frame, blank to ignore
fmap = {'' if fmap is None else fmap} # Fringe map frame, blank to ignore
fpair = {'' if fpair is None else fpair} # File of fringe pairs, ignored if fmap blank

nhalf = {csec['nhalf']} # Half-size of region used for fringe ratios, binned pix, ignored if fmap blank=""
rmin = {csec['rmin']:.2f} # Mininum acceptable individual fmap ratio, ignored if fmap blank
rmax = {csec['rmax']:.2f} # Maximum acceptable individual fmap ratio, ignored if fmap blank
readout = {csec['readout']} # RMS ADU. Float or string name of a file or "!" to estimate on the fly
gain = {csec['gain']} # Gain, electrons/ADU. Float or string name of a file

# The light curve plot which consists of light curves, X & Y
# poistions, the transmission and seeing. All but the light curves can
# be switched off by commenting them out (in full). First a couple of
# general parameters.

[lcplot]
xrange  = {lcsec['xrange']} # maximum range in X to plot (minutes), <= 0 for everything
extend_x = {lcsec['extend_x']} # amount by which to extend xrange, minutes.

# The light curve panel (must be present). Mostly obvious, then a
# series of lines, each starting 'plot' which specify one light curve
# to be plotted giving CCD, target, comparison ('!' if you don't want
# a comparison), an additive offset, a multiplicative scaling factor
# and then a colour for the data and a colour for the error bar There
# will always be a light curve plot, whereas later elements are
# optional, therefore the light curve panel is defined to have unit
# height and all others are scaled relative to this.

[light]
linear  = {lsec['linear']} # linear vertical scale (else magnitudes): 'yes' or 'no'
y_fixed = {lsec['y_fixed']} # keep a fixed vertical range or not: 'yes' or 'no'
y1 = {lsec['y1']} # initial lower y value
y2 = {lsec['y2']} # initial upper y value. y1=y2 for auto scaling
extend_y = {lsec['extend_y']} # fraction of plot height to extend when rescaling

# line or lines defining the targets to plot
{light_plot}
"""
        )

        # position plot (optional)
        if no_position:
            fout.write(
                f"""# The X,Y position panel. Uncomment to show

#[position]
#height = 0.5 # height relative to light curve plot
#x_fixed = no # keep X-position vertical range fixed
#x_min = -5 # lower limit for X-position
#x_max = +5 # upper limit for X-position
#y_fixed = no # keep Y-position vertical range fixed
#y_min = -5 # lower limit for Y-position
#y_max = +5 # upper limit for Y-position
#extend_y = 0.2 # Vertical extension fraction if limits exceeded

## line or lines defining the targets to plot
# plot = 1 1 darkred ! # ccd, targ, dcol, ecol
"""
            )

        else:
            fout.write(
                f"""# The X,Y position panel. Comment *every* line
# if you don't want it.

[position]
height = {psec['height']} # height relative to light curve plot
x_fixed = {psec['x_fixed']} # keep X-position vertical range fixed
x_min = {psec['x_min']} # lower limit for X-position
x_max = {psec['x_max']} # upper limit for X-position
y_fixed = {psec['y_fixed']} # keep Y-position vertical range fixed
y_min = {psec['y_min']} # lower limit for Y-position
y_max = {psec['y_max']} # upper limit for Y-position
extend_y = {psec['extend_y']} # Vertical extension fraction if limits exceeded

# line or lines defining the targets to plot
{position_plot}
"""
            )

        # Transmission plot, optional
        if no_transmission:
            fout.write(
                f"""# The transmission panel. Uncomment if you want it
# Simply plots the flux in whatever apertures are chosen,
# scaling them by their maximum (hence one can sometimes
# find that what you thought was 100% transmission was
# actually only 50% revealed as the cloud clears).

#[transmission]
#height = 0.5 # height relative to the light curve plot
#ymax = 110 # Maximum transmission to plot (>= 100 to slow replotting)

## line or lines defining the targets to plot
# plot = 1 2 green ! # ccd, targ, dcol, ecol
"""
            )
        else:
            fout.write(
                f"""# The transmission panel. If you comment it out,
# comment *every* line.
# Simply plots the flux in whatever apertures are chosen,
# scaling them by their maximum (hence one can sometimes
# find that what you thought was 100% transmission was
# actually only 50% revealed as the cloud clears).

[transmission]
height = {tsec['height']} # height relative to the light curve plot
ymax = {tsec['ymax']} # Maximum transmission to plot (>= 100 to slow replotting)

# line or lines defining the targets to plot
{transmission_plot}
"""
            )

        if no_seeing:
            fout.write(
                f"""# The seeing panel. Uncomment if you want it.
# You can have multiple plot lines. Don't choose linked
# targets as their FWHMs are not measured.

#[seeing]
#height = 0.5 # height relative to the light curve plot
#ymax = 2.999 # Initial maximum seeing
#y_fixed = no # fix the seeing scale (or not)
#extend_y = 0.2 # Y extension fraction if out of range and not fixed

## line or lines defining the targets to plot
# plot = 1 2 green ! # ccd, targ, dcol, ecol
"""
            )

        else:
            fout.write(
                f"""# The seeing panel. Can be commented out if you don't want one but make
# sure to comment it out completely, section name and all
# parameters. You can have multiple plot lines. Don't choose linked
# targets as their FWHMs are not measured.

[seeing]
height = {ssec['height']} # height relative to the light curve plot
ymax = {ssec['ymax']} # Initial maximum seeing
y_fixed = {ssec['y_fixed']} # fix the seeing scale (or not)
extend_y = {ssec['extend_y']} # Y extension fraction if out of range and not fixed

# line or lines defining the targets to plot
{seeing_plot}
"""
            )

        # Finish up
        fout.write(
            f"""
# This option attempts to correct for a badly-positioned focal plane mask
# in drift mode which combined with a high background can lead to steps in
# illumination in the Y direction. This tries to subtract the median in the
# X-direction of each window. 'dthresh' is a threshold used to reject X
# pixels prior to taking the median. The purpose is to prevent the stars
# from distorting the median. Take care with this option which is experimental.

[focal_mask]
demask = {fsec['demask']}
dthresh = {fsec['dthresh']}

# Monitor section. This section allows you to monitor particular
# targets for problems. If they occur, then messages will be printed
# to the terminal during reduce. The messages are determined by the
# bitmask flag set during the extraction of each
# target. Possibilities:

#  NO_FWHM           : no FWHM measured
#  NO_SKY            : no sky pixels at all
#  SKY_AT_EDGE       : sky aperture off edge of window
#  TARGET_AT_EDGE    : target aperture off edge of window
#  TARGET_SATURATED  : at least one pixel in target above saturation level
#  TARGET_NONLINEAR  : at least one pixel in target above nonlinear level
#  NO_EXTRACTION     : no extraction possible
#  NO_DATA           : no valid pixels in aperture

# For a target you want to monitor, type its label, '=', then the
# bitmask patterns you want to be flagged up if they are set. This is
# designed mainly for observing, as there is less you can do once the
# data have been taken, but it still may prove useful.

[monitor]
{monitor}
"""
        )

    print("Reduce file written to {:s}".format(rfile))
