import sys
import os
from time import gmtime, strftime

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

# get hipercam version to write into the reduce file
from pkg_resources import get_distribution, DistributionNotFound
try:
    hipercam_version = get_distribution('hipercam').version
except DistributionNotFound:
    hipercam_version = 'not found'

__all__ = ['genred',]

################################################
#
# genred -- generates a reduce file
#
################################################

def genred(args=None):
    """``genred apfile rfile comment bias flat dark linear [ccd
    smooth_fwhm fwhm fwhm_min rfac rmin rmax sinner souter]``

    Generates a reduce file as needed by `reduce`. You give it the
    name of an aperture file and a few other parameters and it will
    write out a reduce file which you can then refine by hand. A few
    simplifying assumptions are made, e.g. that the target is called '1',
    see below for more. This script effectively defines the format of
    reduce files. The script generates a self-consistent reduce file.
    e.g. if there are apertures in CCD 5, it does not attempt to plot
    any corresponding light curves.

    To avoid excessive prompting, `genred` has many hidden parameters
    that you might want to set up at the start of a run but leave fixed
    thereafter. Specify 'prompt' on the command line to see all of these.
    They are chosen to be the parameters most likely to vary with telescope
    or conditions; many others are left at default values and require
    editing to change.

    Parameters:

        apfile   : string
           the input aperture file created using `setaper` (default extension
           .ape). This will be read for the targets. The main target will be
           assumed to have been called '1', the main comparison '2'. If there
           is a '3' it will be plotted relative to '2'; all others will be
           ignored for plotting purposes. Target '2' will be used to define
           the position and transmission plots for one CCD only [user
           definable]. Target '1' will be used for the seeing plot unless it
           is linked when target '2' will be used instead.

        rfile    : string
           the output reduce file created using `setaper`. This will be read
           for the targets. The main target will be assumed to have been
           called '1', the main comparison '2'. If there is a '3' it will be
           plotted relative to '2'; all others will be ignored for plotting
           purposes.

        comment : string
           comment to add near the top of the reduce file. Obvious things to say
           are the target name and the name of the observer for instance.

        bias    : string
           Name of bias frame; 'none' to ignore.

        flat    : string
           Name of flat field frame; 'none' to ignore.

        dark    : string
           Name of dark frame; 'none' to ignore.

        linear : string
           light curve plot linear (else magnitudes)

        ccd     : string [hidden]
           label of the (single) CCD used for the position plot

        smooth_fwhm : float [hidden]
           FWHM to use for smoothing during initial search

        fwhm     : float [hidden]
           the default FWHM to use when fitting, unbinned pixels.

        fwhm_min : float [hidden]
           the default FWHM to use when fitting, unbinned pixels.

        rfac     : float [hidden]
           target aperture radius relative to the FWHM for 'variable' aperture
           photometry. Usual values 1.5 to 2.5.

        rmin     : float [hidden]
           minimum target aperture radius [unbinned pixels]

        rmax     : float [hidden]
           maximum target aperture radius [unbinned pixels]

        sinner   : float [hidden]
           inner sky aperture radius [unbinned pixels]

        souter   : float [hidden]
           outer sky aperture radius [unbinned pixels]

    """

#    print(my_version)

    command, args = utils.script_args(args)

    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('apfile', Cline.LOCAL, Cline.PROMPT)
        cl.register('rfile', Cline.GLOBAL, Cline.PROMPT)
        cl.register('comment', Cline.LOCAL, Cline.PROMPT)
        cl.register('bias', Cline.LOCAL, Cline.PROMPT)
        cl.register('flat', Cline.LOCAL, Cline.PROMPT)
        cl.register('dark', Cline.LOCAL, Cline.PROMPT)
        cl.register('linear', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.HIDE)
        cl.register('fwhm', Cline.LOCAL, Cline.HIDE)
        cl.register('smooth_fwhm', Cline.LOCAL, Cline.HIDE)
        cl.register('fwhm_min', Cline.LOCAL, Cline.HIDE)
        cl.register('rfac', Cline.LOCAL, Cline.HIDE)
        cl.register('rmin', Cline.LOCAL, Cline.HIDE)
        cl.register('rmax', Cline.LOCAL, Cline.HIDE)
        cl.register('sinner', Cline.LOCAL, Cline.HIDE)
        cl.register('souter', Cline.LOCAL, Cline.HIDE)

        # get inputs

        # the aperture file
        apfile = cl.get_value(
            'apfile', 'aperture input file', cline.Fname('aper.ape',hcam.APER)
        )
        # Read the aperture file
        aper = hcam.MccdAper.read(apfile)

        # the reduce file
        rfile = cl.get_value(
            'rfile', 'reduce output file',
            cline.Fname('reduce.red',hcam.RED,cline.Fname.NEW)
        )

        # user comment string
        comment = cl.get_value(
            'comment', 'user comment to add [<cr> for newline to get multilines]',''
        )
        if comment == '':
            comment = '# There was no user comment\n'
        else:
            comment_lines = comment.split('<cr>')
            comment = '# User comment:\n#\n# ' + '\n# '.join(comment_lines)

        # ones you might quite often want to change

        # bias frame
        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        bias = '' if bias is None else bias

        # flat field frame
        flat = cl.get_value(
            'flat', "flat field frame ['none' to ignore]",
            cline.Fname('flat', hcam.HCAM), ignore='none'
        )
        flat = '' if flat is None else flat

        # dark frame
        dark = cl.get_value(
            'dark', "dark field frame ['none' to ignore]",
            cline.Fname('dark', hcam.HCAM), ignore='none'
        )
        dark = '' if dark is None else dark

        linear = cl.get_value(
            'linear', 'linear light curve plot?', False
        )
        linear = 'yes' if linear else 'no'

        # hidden parameters
        ccd = cl.get_value('ccd', 'label for the CCD used for the position plot','2')
        if ccd not in aper:
            raise HipercamError(
                'CCD {:s} not found in aperture file {:s}'.format(ccd,apfile)
            )
        smooth_fwhm = cl.get_value('smooth_fwhm','search smoothing FWHM [unbinned pixels]',6.,3.)
        fwhm = cl.get_value('fwhm','starting FWHM, unbinned pixels',5.,1.5)
        fwhm_min = cl.get_value('fwhm_min','minimum FWHM, unbinned pixels',1.5,0.)
        rfac = cl.get_value('rfac','target aperture scale factor',1.8,1.0)
        rmin = cl.get_value('rmin','minimum target aperture radius [unbinned pixels]',6.,1.)
        rmax = cl.get_value('rmax','maximum target aperture radius [unbinned pixels]',30.,rmin)
        sinner = cl.get_value('sinner','inner sky aperture radius [unbinned pixels]',30.,rmax)
        souter = cl.get_value('souter','outer sky aperture radius [unbinned pixels]',50.,sinner+1)

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff

    # Generate the extraction lines
    extraction = ''
    for cnam in aper:
        extraction += (
            '{:s} = variable normal'
            ' {:.2f} {:.1f} {:.1f}'
            ' 2.5 {:.1f} {:.1f}'
            ' 3.0 {:.1f} {:.1f}\n').format(
                cnam, rfac, rmin, rmax,
                sinner, sinner, souter, souter
            )

    # standard colours for CCDs
    CCD_COLS = {
        '1' : 'purple',
        '2' : 'green',
        '3' : 'orange',
        '4' : 'red',
        '5' : 'darkred'
    }

    # Generate the light curve plot lines
    light_plot = ''
    for cnam in aper:
        ccdaper = aper[cnam]
        if '1' in ccdaper and '2' in ccdaper:
            light_plot += (
                'plot = {:s} 1 2 0 1 {:10s} !  '
                ' # ccd, targ, comp, off, fac, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
        elif '1' in ccdaper and '2' not in ccdaper:
            light_plot += (
                'plot = {:s} 1 ! 0 1 {:10s} !  '
                ' # ccd, targ, comp, off, fac, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )

        if '2' in ccdaper and '3' in ccdaper:
            light_plot += (
                'plot = {:s} 3 2 0 1 {:10s} !  '
                ' # ccd, targ, domp, off, fac, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )

    # Generate the position plot lines
    position_plot = ''
    ccdaper = aper[ccd]
    if '2' in ccdaper:
        position_plot += (
            'plot = {:s} 2 {:10s} !  '
            ' # ccd, targ, dcol, ecol\n').format(
                ccd, CCD_COLS[ccd]
            )
    elif '3' in ccdaper:
        position_plot += (
            'plot = {:s} 3 {:10s} !  '
            ' # ccd, targ, dcol, ecol\n').format(
                ccd, CCD_COLS[ccd]
            )
    elif '1' in ccdaper:
        position_plot += (
            'plot = {:s} 1 {:10s} !  '
            ' # ccd, targ, dcol, ecol\n').format(
                ccd, CCD_COLS[ccd]
            )
    else:
        raise hcam.HipercamError(
            'Targets 1, 2 and 3 not found; cannot make position plot'
        )

    # Generate the transmission plot lines
    transmission_plot = ''
    for cnam in aper:
        ccdaper = aper[cnam]
        if '2' in ccdaper:
            transmission_plot += (
                'plot = {:s} 2 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
        elif '3' in ccdaper:
            transmission_plot += (
                'plot = {:s} 3 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
        elif '1' in ccdaper:
            transmission_plot += (
                'plot = {:s} 1 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
        else:
            raise hcam.HipercamError(
                'Targets 1, 2 and 3 not found; cannot make transmission plot'
            )

    # Generate the seeing plot lines
    seeing_plot = ''
    for cnam in aper:
        ccdaper = aper[cnam]
        if '1' in ccdaper and not ccdaper['1'].is_linked():
            seeing_plot += (
                'plot = {:s} 2 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
        elif '2' in ccdaper and not ccdaper['2'].is_linked():
            seeing_plot += (
                'plot = {:s} 3 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
        elif '3' in ccdaper  and not ccdaper['3'].is_linked():
            seeing_plot += (
                'plot = {:s} 1 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
        else:
            raise hcam.HipercamError(
                'Targets 1, 2 and 3 not found (or they are linked); cannot make seeing plot'
            )

    # monitor targets (whole lot by default)
    targs = set()
    for cnam in aper:
        ccdaper = aper[cnam]
        for targ in ccdaper:
            targs.add(targ)
    monitor = ''
    for targ in targs:
        monitor += '{:s} = DATA_SATURATED TARGET_OFF_EDGE NO_SKY SKY_OFF_EDGE\n'.format(
            targ
        )

    # time stamp
    tstamp = strftime("%d %b %Y %H:%M:%S (UTC)", gmtime())

    # finally write out the reduce file.
    with open(rfile, 'w') as fout:
        # write out file
        fout.write(
            TEMPLATE.format(
                version=hcam.REDUCE_FILE_VERSION, apfile=apfile,
                fwhm=fwhm, fwhm_min=fwhm_min, extraction=extraction,
                bias=bias, flat=flat, dark=dark,
                smooth_fwhm=smooth_fwhm, linear=linear,
                light_plot=light_plot, position_plot=position_plot,
                transmission_plot=transmission_plot, seeing_plot=seeing_plot,
                monitor=monitor, comment=comment, tstamp=tstamp,
                hipercam_version=hipercam_version,
            )
        )

    print('Reduce file written to {:s}'.format(rfile))

#################################################################
#
# Below is the template that defines reduce files. It is a single
# string with various string format section {} which get replaced
# by the script. All are done as keys to make the format statement
# easier to follow.
#
#################################################################

TEMPLATE = """#
# This is a HiPERCAM "reduce file" which defines the operation of the reduce
# script. It was written by the HiPERCAM pipeline command 'genred'.  Basically
# it consists of a series of sections each of which contains a number of
# parameters. This file explains the meaning of these parameters. The idea is
# that these are to large extent unchanging and it would be annoying to be
# prompted every time for them.
#
# File written on {tstamp}
#
# HiPERCAM pipeline version: {hipercam_version}
#
{comment}


# Start with some general items that tend not to change much. 'version' is the
# version of the reduce file format. It automatically matches 'reduce' at the
# time of creation, but old versions of the reduce file may become incompatible
# with later versions of reduce. Either they will reauire updating to be used,
# or the software version can be rolled back to give a compatible version of
# reduce using 'git'.

[general]
version = {version}  # must be compatible with the version in reduce

ldevice  = 1/xs  # PGPLOT plot device for light curve plots
lwidth   = 0     # light curve plot width, inches, 0 to let program choose
lheight  = 0     # light curve plot height, inches

idevice  = 2/xs  # PGPLOT plot device for image plots [if implot True]
iwidth   = 0     # image curve plot width, inches, 0 to let program choose
iheight  = 0     # image curve plot height, inches

satval   = 65000 # Level at which to flag saturated data.


# The next section defines how the apertures are re-positioned from frame to
# frame. Apertures are re-positioned through a combination of a search near a
# start location followed by a 2D fit. Several parameters below are associated
# with this process. If there are reference apertures, they are located first
# to give a mean shift. The search on non-reference apertures can then be
# tightened.

[apertures]
aperfile   = {apfile}  # file of software apertures for each CCD

search_half_width_ref  = 15   # for initial search around reference aperture, unbinned pixels
search_half_width_non  = 5    # for initial search around non-reference aperture, unbinned pixels
search_smooth_fwhm     = {smooth_fwhm:.1f}    # smoothing FWHM, binned pixels

fit_method     = moffat    # gaussian or moffat
fit_beta       = 4.        # Moffat exponent
fit_fwhm       = {fwhm:.1f}       # FWHM, unbinned pixels
fit_fwhm_min   = {fwhm_min:.1f}       # Minimum FWHM, unbinned pixels
fit_fwhm_fixed = no        # Slightly faster not to fit the FWHM.
fit_half_width = 15        # for fit, unbinned pixels
fit_thresh     = 4         # rejection threshold for fits
fit_height_min = 50        # minimum height to accept a fit


# The next lines define how the apertures will be re-sized and how the flux
# will be extracted from the aperture. There is one line per CCD which starts
# with "CCD label =" and then is followed by
#
# variable | fixed : how the aperture size will be determined. If variable it
# will be scaled relative to the FWHM, so there needs to be a FWHM.
#
# normal | optimal : how the flux will be extracted. 'normal' means a straight
# sum of sky subtracted flux over the aperture. This is often the best choice.
# 'optimum' tries to weight more towards regions of high flux, this is only
# possible if you have fit the profile earlier.
#
# Then follows a series of numbers in three triplets, each of which is a scale
# factor relative to the FWHM for the aperture radius if the 'variable' option
# was chosen, then a minimum and a maximum aperture radius in unbinned pixels.
# The three triples correspond to the innermost target aperture radius, the 
# inner sky radius and finally the outer sky radius.

[extraction]
{extraction}


# Next lines determine how the sky background level is calculated. Note
# you can only set error = variance if method = clipped. 'median' should
# usually be avoided as it can cause noticable steps in light curves. It's
# here as a comparator.
[sky]
method = clipped    # 'clipped' | 'median'
error  = variance   # 'variance' | 'photon': first uses actual variance of sky
thresh = 3.5        # threshold in terms of RMS for 'clipped'


# Calibration frames and constants
[calibration]
crop = yes    # Crop calibrations to match the data
bias = {bias}    # Bias frame, blank to ignore
flat = {flat}    # Flat field frame, blank to ignore
dark = {dark}    # Dark frame, blank to ignore
readout = 3.  # RMS ADU. Float or string name of a file
gain = 1.     # Gain, electrons/ADU. Float or string name of a file


# The light curve plot (includes transmission & seeing as well)

[lcplot]
xrange  = 0    # maximum range in X to plot (minutes), <= 0 for everything
extend_x = 10  # amount by which to extend xrange, minutes.

# light curve panel (must be present). Mostly obvious, then a series of lines,
# each starting 'plot' which specify one light curve to be plotted giving CCD,
# target, comparison ('!' if you don't want a comparison), an additive offset,
# a multiplicative scaling factor and then a colour for the data and a colour
# for the error bar There will always be a light curve plot, whereas later
# elements are optional, therefore the light curve panel is defined to have
# unit height and all others are scaled relative to this.

[light]
linear  = {linear}  # linear vertical scale (else magnitudes): 'yes' or 'no'
y_fixed = no   # keep a fixed vertical range or not: 'yes' or 'no'
y1 = 0         # initial lower y value
y2 = 0         # initial upper y value. y1=y2 for auto scaling
extend_y = 0.1 # fraction of plot height to extend when rescaling

# line or lines defining the targets to plot
{light_plot}


# Configures the position plot. Can be commented out if you don't want one
# but make sure to comment it out completely, section name and all parameters.
# You can have multiple plot lines

[position]
height  = 0.5    # height relative to light curve plot
x_fixed = no     # keep X-position vertical range fixed
x_min   = -5     # lower limit for X-position
x_max   = +5     # upper limit for X-position
y_fixed = no     # keep Y-position vertical range fixed
y_min   = -5     # lower limit for Y-position
y_max   = +5     # upper limit for Y-position
extend_y = 0.2  # Vertical extension fraction if limits exceeded

# line or lines defining the targets to plot
{position_plot}


# Configures the transmission plot. Can be commented out if you don't want one
# but make sure to comment it out completely, section name and all parameters.
# You can have multiple plot lines

[transmission]
height = 0.5      # height relative to the light curve plot
ymax   = 110      # Maximum transmission to plot (>= 100 to slow replotting)

# line or lines defining the targets to plot
{transmission_plot}


# Configures the seeing plot. Can be commented out if you don't want one but
# make sure to comment it out completely, section name and all parameters you
# can have multiple plot lines. Don't choose linked targets as their FWHMs are
# not measured.

[seeing]
height = 0.5          # height relative to the light curve plot
ymax = 1.999          # Initial maximum seeing
y_fixed = yes         # fix the seeing scale (or not)
scale = 0.3           # Arcsec per unbinned pixel
extend_y = 0.2   # Y extension fraction if out of range and not fixed

# line or lines defining the targets to plot
{seeing_plot}

# Monitor section. This section allows you to monitor particular targets
# for problems. If they occur, then messages will be printed to the terminal
# during reduce. The messages are determined by the bitmask flag set during
# the extraction of each target. Ones worth testing for are:
#
#  NO_SKY          : no sky pixels at all
#  SKY_OFF_EDGE    : sky aperture off edge of window
#  TARGET_OFF_EDGE : target aperture off edge of window
#  DATA_SATURATED  : at least one pixel in target aperture saturated
#
# For a target you want to monitor, type its label, '=', then the bitmask
# patterns you want to be flagged up if they are set. This is designed mainly
# for observing, as there is less you can do once the data have been taken.

[monitor]
{monitor}
"""
