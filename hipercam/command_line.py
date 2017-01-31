import sys
import os
import math
from collections import OrderedDict

import numpy as np
from astropy.io import fits

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

def hplot(args=None):
    """
    Plots a multi-CCD image.

    """

    if args is None:
        args = sys.argv[1:]

    # create a Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'hplot', args)

    # register parameters
    cl.register('frame', Cline.LOCAL, Cline.PROMPT)
    cl.register('nccd', Cline.LOCAL, Cline.PROMPT)
    cl.register('nx', Cline.LOCAL, Cline.PROMPT)

    try:
        # get inputs
        frame = cl.get_value('frame', 'frame to plot',
                             cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(frame)
        max_ccd = len(mccd)
        print('max ccd =',max_ccd)
        if max_ccd > 1:
            nccd = cl.get_value('nccd', 'CCD number to plot', 0, 0)
            if nccd == 0:
                nx = cl.get_value('nx', 'number of panels in X', 3, 1)
            else:
                nx = 1
        else:
            nccd = 1

    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    import matplotlib.pyplot as plt

    if nccd == 0:
        ny = int(math.ceil(max_ccd / nx))
        fig = plt.figure()
        i = 1
        for label, ccd in mccd.items():
            if i == 1:
                ax = fig.add_subplot(ny,nx,i)
                asave = ax
            else:
                ax = fig.add_subplot(ny,nx,i,sharex=asave,sharey=asave)
            print(i,label)
            i += 1
            hcam.mpl.pccd(ax, ccd)
            plt.xlim(0.5,ccd.nxtot+0.5)
            plt.ylim(0.5,ccd.nytot+0.5)
            plt.title('CCD ' + str(label))
            plt.xlabel('X pixels')
            plt.ylabel('Y pixels')
            plt.tight_layout(pad=1.01)
    else:
        ccd = mccd[nccd]
        hcam.mpl.pccd(ax, ccd)
        plt.xlim(0.5,ccd.nxtot+0.5)
        plt.ylim(0.5,ccd.nytot+0.5)
        plt.title('CCD ' + str(label))
        plt.xlabel('X pixels')
        plt.ylabel('Y pixels')

    plt.show()

def makefield(args=None):
    """Script to generate an artificial star field which is saved to disk file, a
    first step in generating fake data. A previously generated star field can
    be loaded and added to. The targets are distributed at random, with random
    peak heights based on constant luminosity objects distributed throughout
    3D space. All targets have the same shape thus multiple calls are needed
    to generate a field of objects of multiple shapes. Ellisoidal "Moffat"
    functions [1/(1+r^2)^beta] are used.

    Arguments::

       fname : (string)
          file to add to. Will be created if it does not exist. [input, optional / output]

       ntarg : (int)
          The number of targets to add to the field.

       x1 : (float)
          The left-hand limit of the field [unbinned pixels]

       x2 : (float)
          The right-hand limit of the field [unbinned pixels]

       y1 : (float)
          The lowest limit of the field [unbinned pixels]

       y2 : (float)
          The highest limit of the field [unbinned pixels]

       h1 : (float)
          Minimum peak height [counts per unbinned pixel]

       h2 : (float)
          Maximum peak height [counts per unbinned pixel]

       fwmax : (float)
          Major-axis FWHM [unbinned pixels]

       fwmin : (float)
          Minor-axis FWHM [unbinned pixels]

       angle : (float)
          Angle, anti-clockwise from X-axis [degrees]

       beta : (float)
          Moffat function exponent

    """
    if args is None:
        args = sys.argv[1:]

    # create Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'makefield', args)

    # register parameters
    cl.register('fname', Cline.LOCAL, Cline.PROMPT)
    cl.register('ntarg', Cline.LOCAL, Cline.PROMPT)
    cl.register('x1', Cline.LOCAL, Cline.PROMPT)
    cl.register('x2', Cline.LOCAL, Cline.PROMPT)
    cl.register('y1', Cline.LOCAL, Cline.PROMPT)
    cl.register('y2', Cline.LOCAL, Cline.PROMPT)
    cl.register('h1', Cline.LOCAL, Cline.PROMPT)
    cl.register('h2', Cline.LOCAL, Cline.PROMPT)
    cl.register('fwmax', Cline.LOCAL, Cline.PROMPT)
    cl.register('fwmin', Cline.LOCAL, Cline.PROMPT)
    cl.register('angle', Cline.LOCAL, Cline.PROMPT)
    cl.register('beta', Cline.LOCAL, Cline.PROMPT)

    try:
        # get inputs
        fname = cl.get_value('fname', 'file to save field to',
                                cline.Fname('field', hcam.FIELD, exist=False))
        if os.path.exists(fname):
            # Initialise the field from a file
            field = hcam.Field.rjson(fname)
            print('>> Loaded a field of',len(field),'objects from',fname)
        else:
            # Create an empty field
            field = hcam.Field()
            print('>> Created an empty field.')

        ntarg = cl.get_value('ntarg', 'number of targets', 100, 1)
        x1 = cl.get_value('x1', 'left-hand limit of field', -10.)
        x2 = cl.get_value('x2', 'right-hand limit of field', 2000., x1)
        y1 = cl.get_value('y1', 'lower limit of field', -10.)
        y2 = cl.get_value('y2', 'upper limit of field', 1000., y1)
        h1 = cl.get_value('h1', 'lower peak height limit', 0.1, 1.e-6)
        h2 = cl.get_value('h2', 'upper peak height limit', 1000., h1)
        fwmax = cl.get_value('fwmax', 'FWHM along major axis', 4., 1.e-6)
        fwmin = cl.get_value('fwmin', 'FWHM along minor axis', 4., 1.e-6)
        angle = cl.get_value('angle', 'angle of major axis', 0., -360., 360.)
        beta = cl.get_value('beta', 'Moffat exponent', 4., 1.0)

    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    # add targets
    field.add_random(ntarg, x1, x2, y1, y2, h1, h2, fwmax, fwmin, angle, beta)

    # save result
    field.wjson(fname)

    print('>> Saved a field of',len(field),'objects to',fname)

def makehcam(args=None):
    """Script to generate a fake hipercam frame. Use this to create artificial bias and flat field frames.
    This will split each CCD up into 2x2 windows, so the total dimensions must be multiples of 2.

    Arguments::

       frame : (string)
          name of frame [output]

       nccd : (int)
          number of CCDs

       nxtot : (int)
          total X dimension, a multiple of 2 [unbinned pixels]

       nytot : (int)
          total Y dimension, a multiple of 2 [unbinned pixels]

       mean : (float)
          mean level of frame [counts]

       sigma : (float)
          RMS gaussian noise level [counts]

       gradx : (float)
          change in level in X-direction from left-to-right [counts]

       grady : (float)
          change in level in Y-direction from bottom-to-top [counts]

    """
    if args is None:
        args = sys.argv[1:]

    # create Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'makehcam', args)

    # register parameters
    cl.register('frame', Cline.LOCAL, Cline.PROMPT)
    cl.register('nccd', Cline.LOCAL, Cline.PROMPT)
    cl.register('nxtot', Cline.LOCAL, Cline.PROMPT)
    cl.register('nytot', Cline.LOCAL, Cline.PROMPT)
    cl.register('xbin', Cline.LOCAL, Cline.PROMPT)
    cl.register('ybin', Cline.LOCAL, Cline.PROMPT)
    cl.register('mean', Cline.LOCAL, Cline.PROMPT)
    cl.register('sigma', Cline.LOCAL, Cline.PROMPT)
    cl.register('gradx', Cline.LOCAL, Cline.PROMPT)
    cl.register('grady', Cline.LOCAL, Cline.PROMPT)

    try:
        # get cls
        frame = cl.get_value('frame', 'name of output frame',
                                cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW))
        nccd = cl.get_value('nccd', 'number of CCDs/frame', 5, 1)
        nxtot = cl.get_value('nxtot', 'number of CCDs/frame', 2048, 2, multipleof=2)
        nytot = cl.get_value('nytot', 'number of CCDs/frame', 1024, 2, multipleof=2)
        xbin = cl.get_value('xbin', 'X-binning factor', 1, 1)
        ybin = cl.get_value('ybin', 'Y-binning factor', 1, 1)
        mean = cl.get_value('mean', 'mean counts/pixel', 0.)
        sigma = cl.get_value('sigma', 'RMS counts', 0., 0.)
        gradx = cl.get_value('gradx', 'left-to-right change, counts', 0.)
        grady = cl.get_value('grady', 'bottom-to-top change, counts', 0.)

    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    nx = nxtot // 2 // xbin
    ny = nytot // 2 // ybin

    print('Creating 4 Windows at corners of CCD(s)')
    win1  = hcam.Window(1,1,nx,ny,xbin,ybin)
    win2  = hcam.Window(nxtot-xbin*nx+1,1,nx,ny,xbin,ybin)
    win3  = hcam.Window(1,nytot-ybin*ny+1,nx,ny,xbin,ybin)
    win4  = hcam.Window(nxtot-xbin*nx+1,nytot-ybin*ny+1,nx,ny,xbin,ybin)

    print('Creating equivalent Windats')
    winds = (hcam.Windat(win1), hcam.Windat(win2), hcam.Windat(win3), hcam.Windat(win4))

    def func(x,y):
        xcen = (1+nxtot)/2.
        ycen = (1+nytot)/2.
        return gradx*(x-xcen)/(nxtot-1) + grady*(y-ycen)/(nytot-1)

    print('Adding mean, slope and noise')
    for wind in winds:
        wind.add_fxy(func)
        wind.data += np.random.normal(mean,sigma,(ny,nx))

    # Top-level header
    head = fits.Header()
    head['OBJECT'] = ('Fake star field', 'Target name')
    head['INSTRUME'] = ('HiPERCAM', 'Instrument name')
    head['RUN'] = (3,'Run number')
    head['MEAN'] = (mean, 'mean level used by makehcam')
    head['SIGMA'] = (sigma, 'RMS used by makehcam')
    head['GRADX'] = (gradx, 'Change in X used by makehcam')
    head['GRADY'] = (grady, 'Change in Y used by makehcam')
    head.add_history('Generated by makehcam')

    # Prepare Windats to be made into CCDs
    winds = OrderedDict(zip(range(1,nccd+1),winds))

    ihead = fits.Header()

    print('Creating CCDs')
    ccds = OrderedDict()

    ihead['FILTER'] = ('u','Filter name')
    ccds[1] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('g','Filter name')
    ccds[2] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('r','Filter name')
    ccds[3] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('i','Filter name')
    ccds[4] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('z','Filter name')
    ccds[5] = hcam.CCD(winds,nxtot,nytot,ihead)

    print('Creating the MCCD')
    mccd = hcam.MCCD(ccds, head)

    print('Writing the MCCD')
    mccd.wfits(frame,True)

def makedata(args=None):
    """Script to generate fake data given an artificial star field (generated
    e.g. with `makefield`), a bias and a matching flat field frame. This
    allows a sequence of images to be generated allowing for drifting with
    time, jitter and variable transparency.

    Arguments::

       field : (string)
          star field

       bias : (string)
          bias frame

       flat : (string)
          flat field frame

    """
    if args is None:
        args = sys.argv[1:]

    # create Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'makefield', args)

    # register parameters
    cl.register('fname', Cline.LOCAL, Cline.PROMPT)
    cl.register('ntarg', Cline.LOCAL, Cline.PROMPT)
    cl.register('x1', Cline.LOCAL, Cline.PROMPT)
    cl.register('x2', Cline.LOCAL, Cline.PROMPT)
    cl.register('y1', Cline.LOCAL, Cline.PROMPT)
    cl.register('y2', Cline.LOCAL, Cline.PROMPT)
    cl.register('h1', Cline.LOCAL, Cline.PROMPT)
    cl.register('h2', Cline.LOCAL, Cline.PROMPT)
    cl.register('fwmax', Cline.LOCAL, Cline.PROMPT)
    cl.register('fwmin', Cline.LOCAL, Cline.PROMPT)
    cl.register('angle', Cline.LOCAL, Cline.PROMPT)
    cl.register('beta', Cline.LOCAL, Cline.PROMPT)

    try:
        # get cls
        fname = cl.get_value('fname', 'file to save field to',
                                cline.Fname('field', '.fld', exist=False))
        if os.path.exists(fname):
            # Initialise the field from a file
            field = hcam.Field.rjson(fname)
            print('>> Loaded a field of',len(field),'objects from',fname)
        else:
            # Create an empty field
            field = hcam.Field()
            print('>> Created an empty field.')

        ntarg = cl.get_value('ntarg', 'number of targets', 100, 1)
        x1 = cl.get_value('x1', 'left-hand limit of field', -10.)
        x2 = cl.get_value('x2', 'right-hand limit of field', 2000., x1)
        y1 = cl.get_value('y1', 'lower limit of field', -10.)
        y2 = cl.get_value('y2', 'upper limit of field', 1000., y1)
        h1 = cl.get_value('h1', 'lower peak height limit', 0.1, 1.e-6)
        h2 = cl.get_value('h2', 'upper peak height limit', 1000., h1)
        fwmax = cl.get_value('fwmax', 'FWHM along major axis', 4., 1.e-6)
        fwmin = cl.get_value('fwmin', 'FWHM along minor axis', 4., 1.e-6)
        angle = cl.get_value('angle', 'angle of major axis', 0., -360., 360.)
        beta = cl.get_value('beta', 'Moffat exponent', 4., 1.0)

    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    # add targets
    field.add_random(ntarg, x1, x2, y1, y2, h1, h2, fwmax, fwmin, angle, beta)

    # save result
    field.wjson(fname)

    print('>> Saved a field of',len(field),'objects to',fname)

def wmccd(args=None):
    """
    Creates and writes a fake multi-CCD image with multiple windows to disk
    """

    # Dimensions
    nxtot, nytot = 2048, 1024

    # Create a target field
    print('Creating a target field')
    targs = hcam.Field()
    targs.add_random(200, -5, nxtot+5, -5, nytot+5, 1, 1000, 5., 4., 70, 4)

    print('Creating 5 Windows')
    win1  = hcam.Window(11,20,100,150,1,1)
    win2  = hcam.Window(601,20,100,150,1,1)
    win3  = hcam.Window(21,201,100,150,1,1)
    win4  = hcam.Window(801,20,200,200,1,1)
    win5  = hcam.Window(301,700,200,200,2,2)

    print('Creating 5 Windats')
    winds = (hcam.Windat(win1), hcam.Windat(win2), hcam.Windat(win3),
             hcam.Windat(win4), hcam.Windat(win5))

    print('Adding stars and noise')
    for wind in winds:
        wind.add_fxy(targs)
        wind.add_noise(2.5,1.1)

    # Top-level header
    head = fits.Header()
    head['OBJECT'] = ('Fake star field', 'Target name')
    head['INSTRUME'] = ('HiPERCAM', 'Instrument name')
    head['TELESCOP'] = ('GTC', 'Telescope name')
    head['RUN'] = (3,'Run number')
    head['NEXP'] = (332,'Exposure number')

    # Prepare Windats to be made into CCDs
    winds = OrderedDict(zip((1,2,3,4,5),winds))

    ihead = fits.Header()

    print('Creating CCDs')
    ccds = OrderedDict()

    ihead['FILTER'] = ('u','Filter name')
    ccds[1] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('g','Filter name')
    ccds[2] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('r','Filter name')
    ccds[3] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('i','Filter name')
    ccds[4] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('z','Filter name')
    ccds[5] = hcam.CCD(winds,nxtot,nytot,ihead)

    print('Creating the MCCD')
    mccd = hcam.MCCD(ccds, head)

    print('Writing the MCCD')
    mccd.wfits('mfake.fits',True)

def ptarg(args=None):
    """
    Creates and writes a fake CCD image with
    multiple windows to disk
    """

    # Create a Window
    win = hcam.Window(11,20,100,150,1,1)

    # Create a Windat
    wind = hcam.Windat(win)

    # Create a target field
    targs = hcam.Field()
    targs.add_random(100, 3, 110, 15, 175, 1, 100, 5., 4., 70, 4)
    targs.wjson('targets.json')

    ctargs = hcam.Field.rjson('targets.json')
    print(ctargs)

#    print(type(targs))

    # Add the field to the Windat
#    wind.add_fxy(targs)

#    hcam.mpl.pwind(plt, wind)
#    plt.show()

#    hcam.makefield(20,1,100,1,100,1,100,5,5,0,4)
