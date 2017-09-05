import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta

from trm.pgplot import *
import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ = ['rtplot',]

######################################
#
# rtplot -- display of multiple images
#
######################################

def rtplot(args=None):
    """Plots a sequence of images as a movie in near 'real time', hence
    'rt'. Designed to be used to look at images coming in while at the
    telescope, 'rtplot' comes with many options, a large number of which are
    hidden by default. If you want to see them all, invoke as 'rtplot PROMPT'.

    rtplot can source data from both the ULTRACAM and HiPERCAM servers, from
    local raw ULTRACAM and HiPERCAM files and from lists of HiPERCAM '.hcm'
    files.

    Arguments::

        device  : (string) [hidden]
          Plot device. PGPLOT is used so this should be a PGPLOT-style name,
          e.g. '/xs', '1/xs' etc. At the moment only ones ending /xs are
          supported.

        width   : (float) [hidden]
           plot width (inches). Set = 0 to let the program choose.

        height  : (float) [hidden]
           plot height (inches). Set = 0 to let the program choose. BOTH width
           AND height must be non-zero to have any effect

        source  : (string) [hidden]
           's' = server, 'l' = local, 'f' = file list (ucm for ULTRACAM /
           ULTRASPEC, FITS files for HiPERCAM)

        inst    : (string) [hidden]
           If 's' = server, 'l' = local, name of instrument. Choices: 'u' for
           ULTRACAM / ULTRASPEC, 'h' for HiPERCAM. This is needed because of the
           different formats.

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
           CCD(s) to plot, '0' for all, '1 3' to plot '1' and '3' only, etc.

        nx      : (int)
           number of panels across to display.

        bias    : (string)
           Name of bias frame to subtract, 'none' to ignore.

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

        profit  : (bool) [if plotting a single CCD]
           carry out profile fits or not. If you say yes, then on the first
           plot, you will have the option to pick objects with a cursor. The
           program will then attempt to track these from frame to frame, and
           fit their profile. You may need to adjust 'first' to see anything.
           The parameters used for profile fits are hidden and you may want to
           invoke the command with 'PROMPT' the first time you try profile
           fitting.

        fdevice : (string) [if profit; hidden]
           plot device for profile fits, PGPLOT-style name.
           e.g. '/xs', '2/xs' etc. 

        fwidth  : (float) [hidden]
           fit plot width (inches). Set = 0 to let the program choose.

        fheight : (float) [hidden]
           fit plot height (inches). Set = 0 to let the program choose.
           BOTH fwidth AND fheight must be non-zero to have any effect

        method : (string) [if profit; hidden]
           this defines the profile fitting method, either a gaussian or a
           moffat profile. The latter is usually best.

        beta   : (float) [if profit and method == 'm'; hidden]
           default Moffat exponent

        fwhm   : (float) [if profit; hidden]
           default FWHM, unbinned pixels.

        swidth : (float) [if profit; hidden]
           half width of box for searching for a star, binned pixels.

        smooth : (float) [if profit; hidden]
           FWHM for gaussian smoothing, binned pixels. The initial position
           for fitting is determined by finding the maximum flux in a smoothed
           version of the image in a box of width swidth around the starter
           position.

        fwidth : (float) [if profit; hidden]
           half width of box for profile fit, binned pixels. The fit box is centred
           on the position located by the initial search. It should normally be
           smaller than 'swidth'.


    """

    if args is None:
        args = sys.argv[1:]

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'rtplot', args) as cl:

        # register parameters
        cl.register('device', Cline.LOCAL, Cline.HIDE)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)
        cl.register('source', Cline.LOCAL, Cline.HIDE)
        cl.register('inst', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('twait', Cline.LOCAL, Cline.HIDE)
        cl.register('tmax', Cline.LOCAL, Cline.HIDE)
        cl.register('flist', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('pause', Cline.LOCAL, Cline.HIDE)
        cl.register('nx', Cline.LOCAL, Cline.PROMPT)
        cl.register('bias', Cline.GLOBAL, Cline.PROMPT)
        cl.register('iset', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ilo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ihi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('plo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('phi', Cline.LOCAL, Cline.PROMPT)
        cl.register('profit', Cline.LOCAL, Cline.PROMPT)
        cl.register('fdevice', Cline.LOCAL, Cline.PROMPT)
        cl.register('fwidth', Cline.LOCAL, Cline.HIDE)
        cl.register('fheight', Cline.LOCAL, Cline.HIDE)
        cl.register('method', Cline.LOCAL, Cline.HIDE)
        cl.register('beta', Cline.LOCAL, Cline.HIDE)
        cl.register('fwhm', Cline.LOCAL, Cline.HIDE)
        cl.register('swidth', Cline.LOCAL, Cline.HIDE)
        cl.register('smooth', Cline.LOCAL, Cline.HIDE)
        cl.register('fwidth', Cline.LOCAL, Cline.HIDE)
        cl.register('xlo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('xhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ylo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('yhi', Cline.GLOBAL, Cline.PROMPT)

        # get inputs

        # image plot
        device = cl.get_value('device', 'plot device', '1/xs')
        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

        source = cl.get_value('source', 'data source [s(erver), l(ocal), f(ile list)]',
                              'l', lvals=('s','l','f'))

        if source == 's' or source == 'l':
            inst = cl.get_value('inst', 'instrument [h(ipercam), u(ltracam/spec)]',
                                'h', lvals=('h','u'))
            run = cl.get_value('run', 'run name', 'run005')
            first = cl.get_value('first', 'first frame to plot', 1, 1)

            if source == 's':
                twait = cl.get_value('twait', 'time to wait for a new frame [secs]', 1., 0.)
                tmax = cl.get_value('tmax', 'maximum time to wait for a new frame [secs]', 10., 0.)

        else:
            # set inst = 'h' as only lists of HiPERCAM files are supported
            inst = 'h'
            run = cl.get_value('flist', 'file list', cline.Fname('files.lis',hcam.LIST))
            first = 1

        flist = source == 'f'
        server = source == 's'
        if inst == 'u':
            instrument = 'ULTRA'
        elif inst == 'h':
            instrument = 'HIPER'

        # define the panel grid. first get the labels and maximum dimensions
        ccdinf = hcam.get_ccd_pars(instrument, run, flist)

        try:
            nxdef = cl.get_default('nx')
        except:
            nxdef = 3

        if len(ccdinf) > 1:
            ccd = cl.get_value('ccd', 'CCD(s) to plot [0 for all]', '0')
            if ccd == '0':
                ccds = list(ccdinf.keys())
            else:
                ccds = ccd.split()

            if len(ccds) > 1:
                nxdef = min(len(ccds), nxdef)
                cl.set_default('nx', nxdef)
                nx = cl.get_value('nx', 'number of panels in X', 3, 1)
            else:
                nx = 1
        elif len(ccdinf) == 1:
            nx = 1
            ccds = list(ccdinf.keys())

        cl.set_default('pause', 0.)
        pause = cl.get_value('pause', 'time delay to add between frame plots [secs]', 0., 0.)

        # bias frame (if any)
        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        if bias is not None:
            # read the bias frame
            bframe = hcam.MCCD.rfits(bias)

        # define the display intensities
        iset = cl.get_value(
            'iset', 'set intensity a(utomatically), d(irectly) or with p(ercentiles)?',
            'a', lvals=['a','A','d','D','p','P'])
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == 'd':
            ilo = cl.get_value('ilo', 'lower intensity limit', 0.)
            ihi = cl.get_value('ihi', 'upper intensity limit', 1000.)
        elif iset == 'p':
            plo = cl.get_value('plo', 'lower intensity limit percentile', 5., 0., 100.)
            phi = cl.get_value('phi', 'upper intensity limit percentile', 95., 0., 100.)

        if len(ccds) == 1:
            # many parameters for profile fits, although most are not plotted
            # by default
            profit = cl.get_value('profit', 'do you want profile fits?', False)

            if profit:
                fdevice = cl.get_value('fdevice', 'plot device for fits', '2/xs')
                fwidth = cl.get_value('fwidth', 'fit plot width (inches)', 0.)
                fheight = cl.get_value('fheight', 'fit plot height (inches)', 0.)
                method = cl.get_value('method', 'fit method g(aussian) or m(offat)', 'm', lvals=['g','m'])
                if method == 'm':
                    cl.get_value('beta', 'initial exponent for Moffat fits', 5., 0.5)
                fwhm = cl.get_value('fwhm', 'initial FWHM [unbinned pixels] for profile fits', 6., 2.)
                swidth = cl.get_value(
                    'swidth', 'half width of box for initial location of target [unbinned pixels]', 51., 3.)
                smooth = cl.get_value(
                    'smooth', 'FWHM for smoothing for initial object detection [binned pixels]', 6.)
                fwidth = cl.get_value(
                    'fwidth', 'half width of box for profile fit [unbinned pixels]', 21., 3.)

        else:
            profit = False

        nxmax, nymax = 0, 0
        for cnam in ccds:
            nxtot, nytot = ccdinf[cnam]
            nxmax = max(nxmax, nxtot)
            nymax = max(nymax, nytot)

        xlo = cl.get_value('xlo', 'left-hand X value', 0., 0., nxmax+1)
        xhi = cl.get_value('xhi', 'right-hand X value', float(nxmax), 0., nxmax+1)
        ylo = cl.get_value('ylo', 'lower Y value', 0., 0., nymax+1)
        yhi = cl.get_value('yhi', 'upper Y value', float(nymax), 0., nymax+1)

    # open image plot device
    imdev = hcam.pgp.Device(device)
    if width > 0 and height > 0:
        pgpap(width,height/width)

    # set up panels and axes
    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    # slice up viewport
    pgsubp(nx,ny)

    # plot axes, labels, titles
    for cnam in ccds:
        pgsci(hcam.pgp.Params['axis.ci'])
        pgsch(hcam.pgp.Params['axis.number.ch'])
        pgenv(xlo, xhi, ylo, yhi, 1, 0)
        pglab('X','Y','CCD {:s}'.format(cnam))

    # plot images
    with hcam.data_source(instrument, run, flist, server, first) as spool:

        # 'spool' is an iterable source of MCCDs 
        for n, frame in enumerate(spool):
            print('Exposure {:d}'.format(n+1))

            # display the CCDs chosen
            for nc, cnam in enumerate(ccds):
                ccd = frame[cnam]

                # subtract bias
                if bias is not None:
                    ccd -= bframe[cnam]

                # set to the correct panel and then plot CCD
                ix = (nc % nx) + 1
                iy = nc // nx + 1
                pgpanl(ix,iy)
                hcam.pgp.pccd(ccd,iset,plo,phi,ilo,ihi)

                if n == 0 and profit:

                    class Fpars:
                        """Container class for fit parameters"""
                        def __init__(self, x, y, wind):
                            self.x, self.y, self.wind = x, y, wind

                    # cursor selection of targets after first plot, if profit
                    # accumulate list of starter positions
                    fpos = []
                    print('Please select targets for profile fitting. You can select as many as you like.')
                    x, y, reply = (xlo+xhi)/2, (ylo+yhi)/2, ''
                    while reply != 'Q':
                        print("Place cursor on fit target. Any key to register, 'q' to quit")
                        x, y, reply = pgcurs(x, y)
                        if reply == 'q':
                            break
                        else:
                            # check that the position is inside a window
                            wnam = ccd.inside(x, y, 2)
                            if wnam is not None:
                                # store the position and the window name
                                fpos.append(Fpars(x,y,ccd[wnam]))

                    if len(fpos):
                        # open fit plot device
                        fdev = hcam.pgp.Device(fdevice)
                        if fwidth > 0 and fheight > 0:
                            pgpap(fwidth,fheight/fwidth)

                if profit and len(fpos):
                    # carry out fits
                    for fpar in fpos:
                        # create search region
                        swind = fpar.wind.shrink(
                            fpar.x-swidth, fpar.x+swidth, fpar.y-swidth, fpar.y+swidth)

                        # carry out initial search
                        x,y = detect(swind.data, smooth)

                        # plot location on image
                        imdev.select()
                        pgsci(3)
#                        pgpt1(swind.llx+swind


#    try{
#                                    float sky, peak;
#                                    double xm=x, ym=y;
#                                    Ultracam::profit_init(data[nccd], dvar[nccd], xm, ym, initial_search, fwhm1d, hwidth1d, hwidth, sky, peak, true);#
#
#                                    // Initial estimate of 'a' from FWHM
#                                    double a = 1./2./Subs::sqr(fwhm/Constants::EFAC);

