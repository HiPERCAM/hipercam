import sys
import os
from collections import OrderedDict

import numpy as np
#import matplotlib.pyplot as plt
from astropy.io import fits
import hipercam as hcam
import hipercam.input as inp
from hipercam.input import Input

def makefield(args=None):
    """Entry point script to generate an artificial star field which is saved to
    disk file, a first step in generating fake data. A previously generated
    star field can be loaded and added to. The targets are distributed at
    random, with random peak heights based on constant luminosity objects
    distributed throughout 3D space. All targets have the same shape thus
    multiple calls are needed to generate a field of objects of multiple
    shapes. Ellisoidal "Moffat" functions [1/(1+r^2)^beta] are used. Arguments
    can be sent through as a list of strings. When used as an entry point,
    they will be taken from the command line.  Any arguments not supplied will
    be prompted for.

    Arguments::

       fname : (string)
          file to add to. Will be created if it does not exist.

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

    # create Input object
    input = Input('HIPERCAM_ENV', '.hipercam', 'makefield', args)

    # register parameters
    input.register('fname', Input.LOCAL, Input.PROMPT)
    input.register('ntarg', Input.LOCAL, Input.PROMPT)
    input.register('x1', Input.LOCAL, Input.PROMPT)
    input.register('x2', Input.LOCAL, Input.PROMPT)
    input.register('y1', Input.LOCAL, Input.PROMPT)
    input.register('y2', Input.LOCAL, Input.PROMPT)
    input.register('h1', Input.LOCAL, Input.PROMPT)
    input.register('h2', Input.LOCAL, Input.PROMPT)
    input.register('fwmax', Input.LOCAL, Input.PROMPT)
    input.register('fwmin', Input.LOCAL, Input.PROMPT)
    input.register('angle', Input.LOCAL, Input.PROMPT)
    input.register('beta', Input.LOCAL, Input.PROMPT)

    try:
        # get inputs
        fname = input.get_value('fname', 'file to save field to', 
                                inp.Fname('field', '.fld', exist=False))
        if os.path.exists(fname):
            # Initialise the field from a file
            field = hcam.Field.rjson(fname)
            print('>> Loaded a field of',len(field),'objects from',fname)
        else:
            # Create an empty field
            field = hcam.Field()
            print('>> Created an empty field.')

        ntarg = input.get_value('ntarg', 'number of targets', 100, 1)
        x1 = input.get_value('x1', 'left-hand limit of field', -10.)
        x2 = input.get_value('x2', 'right-hand limit of field', 2000., x1)
        y1 = input.get_value('y1', 'lower limit of field', -10.)
        y2 = input.get_value('y2', 'upper limit of field', 1000., y1)
        h1 = input.get_value('h1', 'lower peak height limit', 0.1, 1.e-6)
        h2 = input.get_value('h2', 'upper peak height limit', 1000., h1)
        fwmax = input.get_value('fwmax', 'FWHM along major axis', 4., 1.e-6)
        fwmin = input.get_value('fwmin', 'FWHM along minor axis', 4., 1.e-6)
        angle = input.get_value('angle', 'angle of major axis', 0., -360., 360.)
        beta = input.get_value('beta', 'Moffat exponent', 4., 1.0)

    except inp.InputError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    # add targets
    field.add_random(ntarg, x1, x2, y1, y2, h1, h2, fwmax, fwmin, angle, beta)

    # save result
    field.wjson(fname)

    print('>> Saved a field of',len(field),'objects to',fname)

def hplot(args=None):
    """
    Creates and plots a fake CCD image with
    multiple windows
    """

    # Dimensions
    nxtot, nytot = 2048, 1024

    # Generate targets
    NTARG = 100
    targs = []
    for nt in range(NTARG):
        x = np.random.uniform(0.,nxtot)
        y = np.random.uniform(0.,nytot)
        h = np.random.uniform(10.,1000.)
        targs.append((x,y,h))

    print('Creating 5 Windows')
    win1  = hcam.Window(11,20,100,150,1,1)
    win2  = hcam.Window(601,20,100,150,1,1)
    win3  = hcam.Window(21,201,100,150,1,1)
    win4  = hcam.Window(801,20,200,200,1,1)
    win5  = hcam.Window(301,400,200,200,2,2)

    print('Creating 5 Windats')
    wind1 = hcam.Windat(win1)
    wind2 = hcam.Windat(win2)
    wind3 = hcam.Windat(win3)
    wind4 = hcam.Windat(win4)
    wind5 = hcam.Windat(win5)

    print('Adding stars')
    wind1.add_stars(targs, 2.5, 3., 0.05)
    wind1.add_noise(2.5,1.1)
    wind2.add_stars(targs, 2.5, 3., 0.05)
    wind2.add_noise(2.5,1.1)
    wind3.add_stars(targs, 2.5, 3., 0.05)
    wind3.add_noise(2.5,1.1)
    wind4.add_stars(targs, 2.5, 3., 0.05)
    wind4.add_noise(2.5,1.1)
    wind5.add_stars(targs, 2.5, 3., 0.05)
    wind5.add_noise(2.5,1.1)

    print('Creating a CCD')
    ccd = hcam.CCD({1: wind1, 2: wind2, 3: wind3, 4: wind4, 5: wind5},
                       nxtot, nytot)

    print('Plotting the CCD')
    hcam.mpl.pccd(plt, ccd)
    plt.xlim(0.5,nxtot+0.5)
    plt.ylim(0.5,nytot+0.5)
    plt.show()

def wccd(args=None):
    """
    Creates and writes a fake CCD image with
    multiple windows to disk
    """

    # Dimensions
    nxtot, nytot = 2048, 1024

    # Generate targets
    NTARG = 100
    targs = []
    for nt in range(NTARG):
        x = np.random.uniform(0.,nxtot)
        y = np.random.uniform(0.,nytot)
        h = np.random.uniform(10.,1000.)
        targs.append((x,y,h))

    print('Creating 5 Windows')
    win1  = hcam.Window(11,20,100,150,1,1)
    win2  = hcam.Window(601,20,100,150,1,1)
    win3  = hcam.Window(21,201,100,150,1,1)
    win4  = hcam.Window(801,20,200,200,1,1)
    win5  = hcam.Window(301,400,200,200,2,2)

    print('Creating 5 Windats')
    wind1 = hcam.Windat(win1)
    wind2 = hcam.Windat(win2)
    wind3 = hcam.Windat(win3)
    wind4 = hcam.Windat(win4)
    wind5 = hcam.Windat(win5)

    print('Adding stars')
    wind1.add_stars(targs, 2.5, 3., 0.05)
    wind1.add_noise(2.5,1.1)
    wind2.add_stars(targs, 2.5, 3., 0.05)
    wind2.add_noise(2.5,1.1)
    wind3.add_stars(targs, 2.5, 3., 0.05)
    wind3.add_noise(2.5,1.1)
    wind4.add_stars(targs, 2.5, 3., 0.05)
    wind4.add_noise(2.5,1.1)
    wind5.add_stars(targs, 2.5, 3., 0.05)
    wind5.add_noise(2.5,1.1)

    head = fits.Header()
    head['OBJECT'] = 'Fake star field'

    print('Creating a CCD')
    ccd = hcam.CCD(
        OrderedDict(
            {1: wind1, 2: wind2, 3: wind3, 4: wind4, 5: wind5}),
        nxtot, nytot, head)

    ccd.wfits('fake.fits', overwrite=True)


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
