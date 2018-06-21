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
    """``genred apfile rfile comment bias flat dark linear inst [ncpu extendx
    ccd location smoothfwhm method beta betamax fwhm fwhmmin searchwidth thresh
    hminref hminnrf rfac rmin rmax sinner souter scale]``

    Generates a reduce file as needed by |reduce|. You give it the name of an
    aperture file and a few other parameters and it will write out a reduce
    file which you can then refine by hand. A few simplifying assumptions are
    made, e.g. that the target is called '1'; see below for more. This script
    effectively defines the format of reduce files. The script attempts above
    all to generate a self-consistent reduce file.  e.g. if there are
    apertures in CCD 5, it does not attempt to plot any corresponding light
    curves.

    To avoid excessive prompting, |genred| has many hidden parameters. The
    very first time you use it on a run, specify ``prompt`` on the command line
    to see all of these.  They are chosen to be the parameters most likely to
    vary with telescope or conditions; many others are left at default values
    and require editing to change. If you find yourself repeatedly editing a
    parameter, let me know and I will add it to this routine.

    Parameters:

        apfile   : string
           the input aperture file created using |setaper| (default extension
           .ape). This will be read for the targets. The main target will be
           assumed to have been called '1', the main comparison '2'. If there
           is a '3' it will be plotted relative to '2'; all others will be
           ignored for plotting purposes. Target '2' will be used to define
           the position and transmission plots for one CCD only [user
           definable]. Target '1' will be used for the seeing plot unless it
           is linked when target '2' will be used instead.

        rfile    : string
           the output reduce file created using |setaper|. This will be read
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

        linear  : string
           light curve plot linear (else magnitudes)

        inst    : string
           the instrument (needed to set nonlinearity and saturation levels for
           warning purposes. Possible options listed.

        ncpu    : int [hidden]
           some increase in speed can be obtained by running the reduction in
           parallel. This parameter is the number of CPUs to use. The
           parellisation is over the CCDs, and there is no point having ncpu
           greater than the number of CCDs.

        extendx : float [hidden]
           how many minutes to extend light curve plot by at a time

        ccd     : string [hidden]
           label of the (single) CCD used for the position plot

        location : string [hidden]
           whether to reposition apertures or leave them fixed.

        smoothfwhm : float [hidden]
           FWHM to use for smoothing during initial search

        method   : string
           profile fitting method. 'g' for gaussian, 'm' for moffat

        beta     : float [hidden]
           default Moffat exponent to use to start fitting

        betamax  : float [hidden]
           maximum Moffat exponent to pass on to subsequent fits. Prevents
           it wandering to silly values which can happen.

        fwhm     : float [hidden]
           the default FWHM to use when fitting, unbinned pixels.

        fwhmmin  : float [hidden]
           the default FWHM to use when fitting, unbinned pixels.

        searchwidth : int [hidden]
           half width in (binned) pixels for the target searches

        fitwidth : int [hidden]
           half width in (binned) pixels for the profile fits

        maxshift : float [hidden]
           maximum shift of non-reference targets relative to the initial
           positions derived from the reference targets. The reference targets
           should give good initial positions, thus this can be set to quite a
           small value to improve robustness against duff positions, caused
           e.g. by cosmic rays or field crowding in periods of bad seeing. Use
           linked apertures for more difficult cases still.

        thresh : float [hidden]
           RMS rejection threshold for profile fits.

        hminref : float [hidden]
           minimum peak height for a fit to a reference aperture to be accepted

        hminnrf : float [hidden]
           minimum peak height for a fit to a non-reference aperture to be accepted

        rfac : float [hidden]
           target aperture radius relative to the FWHM for 'variable' aperture
           photometry. Usual values 1.5 to 2.5.

        rmin : float [hidden]
           minimum target aperture radius [unbinned pixels]

        rmax     : float [hidden]
           maximum target aperture radius [unbinned pixels]

        sinner   : float [hidden]
           inner sky aperture radius [unbinned pixels]

        souter   : float [hidden]
           outer sky aperture radius [unbinned pixels]

        readout : float | string [hidden]
           readout noise, RMS ADU. Can either be a single value or an hcm file.

        gain : float [hidden]
           gain, electrons per ADU. Can either be a single value or an hcm file.

        scale    : float [hidden]
           image scale in arcsec/pixel

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
        cl.register('inst', Cline.LOCAL, Cline.HIDE)
        cl.register('ncpu', Cline.LOCAL, Cline.HIDE)
        cl.register('extendx', Cline.LOCAL, Cline.HIDE)
        cl.register('ccd', Cline.LOCAL, Cline.HIDE)
        cl.register('location', Cline.LOCAL, Cline.HIDE)
        cl.register('beta', Cline.LOCAL, Cline.HIDE)
        cl.register('betamax', Cline.LOCAL, Cline.HIDE)
        cl.register('fwhm', Cline.LOCAL, Cline.HIDE)
        cl.register('smoothfwhm', Cline.LOCAL, Cline.HIDE)
        cl.register('method', Cline.LOCAL, Cline.HIDE)
        cl.register('fwhmmin', Cline.LOCAL, Cline.HIDE)
        cl.register('searchwidth', Cline.LOCAL, Cline.HIDE)
        cl.register('fitwidth', Cline.LOCAL, Cline.HIDE)
        cl.register('maxshift', Cline.LOCAL, Cline.HIDE)
        cl.register('thresh', Cline.LOCAL, Cline.HIDE)
        cl.register('hminref', Cline.LOCAL, Cline.HIDE)
        cl.register('hminnrf', Cline.LOCAL, Cline.HIDE)
        cl.register('rfac', Cline.LOCAL, Cline.HIDE)
        cl.register('rmin', Cline.LOCAL, Cline.HIDE)
        cl.register('rmax', Cline.LOCAL, Cline.HIDE)
        cl.register('sinner', Cline.LOCAL, Cline.HIDE)
        cl.register('souter', Cline.LOCAL, Cline.HIDE)
        cl.register('readout', Cline.LOCAL, Cline.HIDE)
        cl.register('gain', Cline.LOCAL, Cline.HIDE)
        cl.register('scale', Cline.LOCAL, Cline.HIDE)

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
            'comment', 'user comment to add [<cr>'
            ' for newline to get multilines]',''
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

        inst = cl.get_value(
            'inst', 'instrument (hipercam, ultracam, ultraspec)',
            'hipercam', lvals=['hipercam', 'ultracam', 'ultraspec','ignore']
        )

        if inst == 'hipercam':
            warn_levels = """# Warning levels for instrument = HiPERCAM
warn = 1 50000 64000
warn = 2 50000 64000
warn = 3 50000 64000
warn = 4 50000 64000
warn = 5 50000 64000
"""
            maxcpu = 5
        elif inst == 'ultracam':
            warn_levels = """# Warning levels for instrument = ULTRACAM
warn = 1 28000 64000
warn = 2 28000 64000
warn = 3 50000 64000
"""
            maxcpu = 3
        elif inst == 'ultraspec':
            warn_levels = """# Warning levels for instrument = ULTRASPEC
warn = 1 60000 64000
"""
            maxcpu = 1

        else:
            warn_levels = """# No warning levels have been set!!"""
            maxcpu = 20

        if maxcpu > 1:
            ncpu = cl.get_value(
                'ncpu', 'number of CPUs to use (<= number of CCDs)',
                1, 1, maxcpu
            )
        else:
            ncpu = 1

        linear = cl.get_value(
            'linear', 'linear light curve plot?', False
        )
        linear = 'yes' if linear else 'no'

        # hidden parameters

        extendx = cl.get_value(
            'extendx', 'how much to extend light curve plot [mins]',
            10.,0.01
        )

        ccd = cl.get_value(
            'ccd', 'label for the CCD used for the position plot','2'
        )
        if ccd not in aper:
            raise hcam.HipercamError(
                'CCD {:s} not found in aperture file {:s}'.format(ccd,apfile)
            )

        # hidden parameters
        location = cl.get_value(
            'location', 'aperture location, f(ixed) or v(ariable)',
            'v', lvals=['f','v']
        )
        location = 'variable' if location == 'v' else 'fixed'
        comm_seeing = '#' if location == 'fixed' else ''
        comm_position = '#' if location == 'fixed' else ''

        smooth_fwhm = cl.get_value(
            'smoothfwhm','search smoothing FWHM [binned pixels]',6.,3.
        )

        profile_type = cl.get_value(
            'method', 'profile fit method, g(aussian) or m(offat)',
            'g', lvals=['g','m']
        )
        profile_type = 'gaussian' if profile_type == 'g' else 'moffat'

        beta = cl.get_value(
            'beta','starting value of beta', 4., 3.
        )

        beta_max = cl.get_value(
            'betamax','maximum value of beta to start consecutive fits',
            20., beta
        )

        fwhm = cl.get_value(
            'fwhm','starting FWHM, unbinned pixels',5.,1.5
        )

        fwhm_min = cl.get_value(
            'fwhmmin','minimum FWHM, unbinned pixels',1.5,0.
        )

        search_half_width = cl.get_value(
            'searchwidth', 'half width for initial searches, unbinned pixels', 11, 3
        )

        fit_half_width = cl.get_value(
            'fitwidth', 'half width for profile fits, unbinned pixels', 21, 5
        )

        fit_max_shift = cl.get_value(
            'maxshift', 'maximum non-reference shift, unbinned pixels', 15, 0
        )

        thresh = cl.get_value(
            'thresh', 'RMS rejection threshold for fits (sigma)', 5., 2.
        )

        height_min_ref = cl.get_value(
            'hminref',
            'minimum peak height for a fit to reference aperture [counts]', 50, 1
        )

        height_min_nrf = cl.get_value(
            'hminnrf',
            'minimum peak height for a fit to non-reference aperture [counts]', 10, 1
        )

        rfac = cl.get_value(
            'rfac','target aperture scale factor',1.8,1.0
        )

        rmin = cl.get_value(
            'rmin','minimum target aperture radius [unbinned pixels]',6.,1.
        )

        rmax = cl.get_value(
            'rmax','maximum target aperture radius [unbinned pixels]',30.,rmin
        )

        sinner = cl.get_value(
            'sinner','inner sky aperture radius [unbinned pixels]',30.,rmax
        )

        souter = cl.get_value(
            'souter','outer sky aperture radius [unbinned pixels]',50.,sinner+1
        )


        readout = cl.get_value(
            'readout', 'readout noise, RMS ADU (float or file name)', '4.5'
        )

        gain = cl.get_value(
            'gain', 'gain, electrons/ADU, (float or file name)', '1.1'
        )

        scale = cl.get_value(
            'scale','image scale [arcsec/unbinned pixel]',0.3,0.001
        )

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff

    # Generate the extraction lines. Note that the aperture location
    # parameter maps into the same names as the aperture re-size
    # parameter
    extraction = ''
    for cnam in aper:
        extraction += (
            '{:s} = {:s} normal'
            ' {:.2f} {:.1f} {:.1f}'
            ' 2.5 {:.1f} {:.1f}'
            ' 3.0 {:.1f} {:.1f}\n').format(
                cnam, location, rfac, rmin, rmax,
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
    no_light = True
    for cnam in aper:
        ccdaper = aper[cnam]
        if '1' in ccdaper and '2' in ccdaper:
            light_plot += (
                'plot = {:s} 1 2 0 1 {:10s} !  '
                ' # ccd, targ, comp, off, fac, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
            no_light = False

        elif '1' in ccdaper and '2' not in ccdaper:
            light_plot += (
                'plot = {:s} 1 ! 0 1 {:10s} !  '
                ' # ccd, targ, comp, off, fac, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
            no_light = False

        if '2' in ccdaper and '3' in ccdaper:
            light_plot += (
                'plot = {:s} 3 2 0 1 {:10s} !  '
                ' # ccd, targ, domp, off, fac, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
            no_light = False

    if no_light:
        raise hcam.HipercamError(
            'Found no targets for light curve plots in any CCD; cannot make light curve plot'
        )

    # Generate the position plot lines
    position_plot = ''
    ccdaper = aper[ccd]
    no_position = True
    if '2' in ccdaper:
        position_plot += (
            '{:s}plot = {:s} 2 {:10s} !  '
            ' # ccd, targ, dcol, ecol\n').format(
                comm_position, ccd, CCD_COLS[ccd]
            )
        no_position = False

    elif '3' in ccdaper:
        position_plot += (
            '{:s}plot = {:s} 3 {:10s} !  '
            ' # ccd, targ, dcol, ecol\n').format(
                comm_position, ccd, CCD_COLS[ccd]
            )
        no_position = False

    elif '1' in ccdaper:
        position_plot += (
            '{:s}plot = {:s} 1 {:10s} !  '
            ' # ccd, targ, dcol, ecol\n').format(
                comm_position, ccd, CCD_COLS[ccd]
            )
        no_position = False

    if no_position:
        raise hcam.HipercamError(
            'Targets 1, 2 and 3 not found in '
            'CCD = {:s}; cannot make position plot'.format(ccd)
        )

    # Generate the transmission plot lines
    transmission_plot = ''
    no_transmission = True
    for cnam in aper:
        ccdaper = aper[cnam]
        if '2' in ccdaper:
            transmission_plot += (
                'plot = {:s} 2 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
            no_transmission = False

        elif '3' in ccdaper:
            transmission_plot += (
                'plot = {:s} 3 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
            no_transmission = False

        elif '1' in ccdaper:
            transmission_plot += (
                'plot = {:s} 1 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    cnam, CCD_COLS[cnam]
                )
            no_transmission = False

    if no_transmission:
        raise hcam.HipercamError(
            'Targets 1, 2 and 3 not found in any CCDs;'
            ' cannot make transmission plot'
        )

    # Generate the seeing plot lines
    seeing_plot = ''
    no_seeing = True
    for cnam in aper:
        ccdaper = aper[cnam]
        if '1' in ccdaper and not ccdaper['1'].is_linked():
            seeing_plot += (
                '{:s}plot = {:s} 1 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    comm_seeing, cnam, CCD_COLS[cnam]
                )
            no_seeing = False

        elif '2' in ccdaper and not ccdaper['2'].is_linked():
            seeing_plot += (
                '{:s}plot = {:s} 2 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    comm_seeing, cnam, CCD_COLS[cnam]
                )
            no_seeing = False

        elif '3' in ccdaper  and not ccdaper['3'].is_linked():
            seeing_plot += (
                '{:s}plot = {:s} 3 {:10s} !  '
                ' # ccd, targ, dcol, ecol\n').format(
                    comm_seeing, cnam, CCD_COLS[cnam]
                )
            no_seeing = False

    if no_seeing:
        raise hcam.HipercamError(
            'Targets 1, 2 and 3 not found in any CCD'
            ' (or they are linked); cannot make seeing plot'
        )

    # monitor targets (whole lot by default)
    targs = set()
    for cnam in aper:
        ccdaper = aper[cnam]
        for targ in ccdaper:
            targs.add(targ)
    monitor = ''
    for targ in sorted(targs):
        monitor += ('{:s} = NO_EXTRACTION TARGET_SATURATED TARGET_AT_EDGE'
                    ' TARGET_NONLINEAR NO_SKY SKY_AT_EDGE\n').format(targ)

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
                hipercam_version=hipercam_version, location=location,
                comm_seeing=comm_seeing, extendx=extendx,
                comm_position=comm_position, scale=scale,
                warn_levels=warn_levels, ncpu=ncpu,
                search_half_width=search_half_width,
                fit_half_width=fit_half_width, profile_type=profile_type,
                height_min_ref=height_min_ref, height_min_nrf=height_min_nrf,
                beta=beta, beta_max=beta_max, thresh=thresh, readout=readout,
                gain=gain, fit_max_shift=fit_max_shift
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
# This is a HiPERCAM "reduce file" which defines the operation
# of the reduce script. It was written by the HiPERCAM pipeline command
# 'genred'.  It consists of a series of sections each of which contains a
# number of parameters. The file is self-documenting on the meaning of these
# parameters. The idea is that these are to large extent unchanging and it
# would be annoying to be prompted every time for them, but it also acts as a
# record of how reduction was carried out and is fed into the log file produce
# by 'reduce'.
#
# File written on {tstamp}
#
# HiPERCAM pipeline version: {hipercam_version}
#
{comment}

# Start with some general items that tend not to change much. 'version' is the
# version of the reduce file format which changes more slowly than the
# software does. It must match the same parameter in 'reduce' for reduction to
# proceed. This is automatically the case at the time of creation, but old
# versions of the reduce file may become incompatible with later versions of
# reduce. Either they will require updating to be used, or the software
# version can be rolled back to give a compatible version of reduce using
# 'git'.

[general]
version = {version} # must be compatible with the version in reduce

ldevice = 1/xs # PGPLOT plot device for light curve plots
lwidth = 0 # light curve plot width, inches, 0 to let program choose
lheight = 0 # light curve plot height, inches

idevice = 2/xs # PGPLOT plot device for image plots [if implot True]
iwidth = 0 # image curve plot width, inches, 0 to let program choose
iheight = 0 # image curve plot height, inches

# series of count levels at which warnings will be triggered for (a) non
# linearity and (b) saturation. Each line starts 'warn =', and is then
# followed by the CCD label, the non-linearity level and the saturation level

{warn_levels}

# The aperture reposition and extraction stages can be run in separate CPUs in
# parallel for each CCD. The next parameter is the number of CPUs to use. The
# maximum useful number is the number of CCDs in the instrument, e.g. 5 for
# HiPERCAM. You probably also want to leave at least one CPU to do other stuff,
# but if you have more than 2 CPUs, this parameter may help speed things.
ncpu = {ncpu}

# The next section defines how the apertures are re-positioned from frame to
# frame. Apertures are re-positioned through a combination of a search near a
# start location followed by a 2D profile fit. Several parameters below are
# associated with this process and setting these right can be the key to a
# successful reduction. If there are reference apertures, they are located
# first to give a mean shift. This is used to bypass the initial search on any
# non-reference apertures. The search is carried out by first extracting a
# square sub-window centred on the last good position of a target. This is
# then smoothed by a gaussian, and the peak is taken as the initial position
# for later profile fits. The gaussian should have a width comparable to the
# targets FWHM and is to make the process more robust against cosmic rays. The
# width of the search box depends on how good the telescope guiding is. In
# particular if at some point there is a sudden jump in position, the box may
# have to be large enough to cope. Well-chosen reference targets, which above
# all should be isolated, can help this process a great deal. The boxes for
# the fits need to be large enough to include the target and a bit of sky to
# ensure that the FWHM is accurately measured. Remember that seeing can flare
# of course. If your target was defocussed, a gaussian or Moffat function will
# be a poor fit and you may be better keeping the FWHM fixed at a large value
# comparable to the widths of your defeoccused images. If the apertures are
# chosen to be fixed, there will be no search or fit carried out in which case
# you must choose 'fixed' as well when it comes the extraction since otherwise
# it needs a FWHM. The most obscure parameter below is probably 'fit_ndiv'. If
# this is made > 0, the fit routine attempts to allow for pixellation by
# evaluating the profile at multiple point in each pixel of the fit. First it
# will evaluate the profile for every unbinned pixel with a binned pixel if
# the pixels are binned; second, it will evaluate the profile over an ndiv by
# ndiv square grid within each unbinned pixel. Obviously this will slow
# things, but it could help if your images are under-sampled. Finally if you
# use reference targets, you should get good initial positions for the
# non-reference targets. You can then guard against problems using the
# parameter 'fit_max_shift' to reject positions that shift too far from the
# initial guess.

[apertures]
aperfile = {apfile} # file of software apertures for each CCD
location = {location} # aperture locations: 'fixed' or 'variable'

search_half_width = {search_half_width:d} # for initial search for objects around previous position, unbinned pixels
search_smooth_fwhm = {smooth_fwhm:.1f}    # smoothing FWHM, binned pixels

fit_method = {profile_type} # gaussian or moffat
fit_beta = {beta:.1f} # Moffat exponent
fit_beta_max = {beta_max:.1f} # max Moffat expt for later fits
fit_fwhm = {fwhm:.1f} # FWHM, unbinned pixels
fit_fwhm_min = {fwhm_min:.1f} # Minimum FWHM, unbinned pixels
fit_ndiv = 0 # sub-pixellation factor
fit_fwhm_fixed = no # Might want to set = 'yes' for defocussed images
fit_half_width = {fit_half_width:d} # for fit, unbinned pixels
fit_thresh = {thresh:.2f} # RMS rejection threshold for fits
fit_height_min_ref = {height_min_ref:.1f} # minimum height to accept a fit, reference aperture
fit_height_min_nrf = {height_min_nrf:.1f} # minimum height to accept a fit, non-reference aperture
fit_max_shift = {fit_max_shift:.1f} # max. non-ref. shift, unbinned pixels.

# The next lines define how the apertures will be re-sized and how the flux
# will be extracted from the aperture. There is one line per CCD with format:
#
# <CCD label> = <resize> <extract method> [scale min max] [scale min max]
#               [scale min max]
#
# where: <CCD label> is the CCD label; <resize> is either 'variable' or
# 'fixed' and sets how the aperture size will be determined -- if variable it
# will be scaled relative to the FWHM, so profile fitting will be attempted;
# <extract method> is either 'normal' or 'optimal' to say how the flux will be
# extracted -- 'normal' means a straight sum of sky subtracted flux over the
# aperture, 'optimal' use Tim Naylor's profile weighting, and requires profile
# fits to work. Finally there follow a series of numbers in three triplets,
# each of which is a scale factor relative to the FWHM for the aperture radius
# if the 'variable' option was chosen, then a minimum and a maximum aperture
# radius in unbinned pixels.  The three triples correspond to the target
# aperture radius, the inner sky radius and finally the outer sky radius. The
# mininum and maximum also apply if you choose 'fixed' apertures and can be
# used to override whatever value comes from the aperture file. A common
# approach is set them equal to each other to give a fixed value, especially
# for the sky where one does not necessarily want the radii to vary.

[extraction]
{extraction}

# Next lines determine how the sky background level is calculated. Note
# you can only set error = variance if method = 'clipped'. 'median' should
# usually be avoided as it can cause noticable steps in light curves. It's
# here as a comparator.

[sky]
method = clipped # 'clipped' | 'median'
error  = variance # 'variance' | 'photon': first uses actual variance of sky
thresh = 3. # threshold in terms of RMS for 'clipped'

# Calibration frames and constants
[calibration]
crop = yes # Crop calibrations to match the data
bias = {bias} # Bias frame, blank to ignore
flat = {flat} # Flat field frame, blank to ignore
dark = {dark} # Dark frame, blank to ignore
readout = {readout} # RMS ADU. Float or string name of a file
gain = {gain} # Gain, electrons/ADU. Float or string name of a file

# The light curve plot which consists of light curves, X & Y poistions,
# the transmission and seeing. All but the light curves can be switched
# off by commenting them out (in full). First a couple of general
# parameters.

[lcplot]
xrange  = 0 # maximum range in X to plot (minutes), <= 0 for everything
extend_x = {extendx:.2f} # amount by which to extend xrange, minutes.

# The light curve panel (must be present). Mostly obvious, then a series of
# lines, each starting 'plot' which specify one light curve to be plotted
# giving CCD, target, comparison ('!' if you don't want a comparison), an
# additive offset, a multiplicative scaling factor and then a colour for the
# data and a colour for the error bar There will always be a light curve plot,
# whereas later elements are optional, therefore the light curve panel is
# defined to have unit height and all others are scaled relative to this.

[light]
linear  = {linear} # linear vertical scale (else magnitudes): 'yes' or 'no'
y_fixed = no # keep a fixed vertical range or not: 'yes' or 'no'
y1 = 0 # initial lower y value
y2 = 0 # initial upper y value. y1=y2 for auto scaling
extend_y = 0.1 # fraction of plot height to extend when rescaling

# line or lines defining the targets to plot
{light_plot}


# The X,Y position panel. Can be commented out if you don't want it but make
# sure to comment it out completely, section name and all parameters.  You can
# have multiple plot lines.

{comm_position}[position]
{comm_position}height = 0.5 # height relative to light curve plot
{comm_position}x_fixed = no # keep X-position vertical range fixed
{comm_position}x_min = -5 # lower limit for X-position
{comm_position}x_max = +5 # upper limit for X-position
{comm_position}y_fixed = no # keep Y-position vertical range fixed
{comm_position}y_min = -5 # lower limit for Y-position
{comm_position}y_max = +5 # upper limit for Y-position
{comm_position}extend_y = 0.2 # Vertical extension fraction if limits exceeded

# line or lines defining the targets to plot
{position_plot}

# The transmission panel. Can be commented out if you don't want one but make
# sure to comment it out completely, section name and all parameters.  You can
# have multiple plot lines. This simply plots the flux in whatever apertures
# are chosen, scaling them by their maximum (hence one can sometimes find that
# what you thought was 100% transmission was actually only 50% revealed as the
# cloud clears).

[transmission]
height = 0.5 # height relative to the light curve plot
ymax = 110 # Maximum transmission to plot (>= 100 to slow replotting)

# line or lines defining the targets to plot
{transmission_plot}

# The seeing plot. Can be commented out if you don't want one but make sure to
# comment it out completely, section name and all parameters. You can have
# multiple plot lines. Don't choose linked targets as their FWHMs are not
# measured.

{comm_seeing}[seeing]
{comm_seeing}height = 0.5 # height relative to the light curve plot
{comm_seeing}ymax = 1.999 # Initial maximum seeing
{comm_seeing}y_fixed = no # fix the seeing scale (or not)
{comm_seeing}scale = {scale:.2f} # Arcsec per unbinned pixel
{comm_seeing}extend_y = 0.2 # Y extension fraction if out of range and not fixed

# line or lines defining the targets to plot
{seeing_plot}

# Monitor section. This section allows you to monitor particular targets
# for problems. If they occur, then messages will be printed to the terminal
# during reduce. The messages are determined by the bitmask flag set during
# the extraction of each target. Ones worth testing for are:
#
#  NO_SKY            : no sky pixels at all
#  SKY_AT_EDGE       : sky aperture off edge of window
#  TARGET_AT_EDGE    : target aperture off edge of window
#  TARGET_SATURATED  : at least one pixel in target above saturation level
#  TARGET_NONLINEAR  : at least one pixel in target above nonlinear level
#
# For a target you want to monitor, type its label, '=', then the bitmask
# patterns you want to be flagged up if they are set. This is designed mainly
# for observing, as there is less you can do once the data have been taken, but
# it still may prove useful.

[monitor]
{monitor}
"""
