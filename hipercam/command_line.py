# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import sys
import argparse
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from astropy.io import fits
import hipercam as hcam
import hipercam.input as inp

def makefield(args=None):
    """Generate an artificial star field which is saved to disk file, a first step
    in generating fake data. A previously generated star field can be loaded
    and added to.

    """

    # create Input object
    input = inp.Input('HIPERCAM_ENV', '.hipercam', args)

    # register parameters
    input.register('ntarg', inp.LOCAL, inp.PROMPT)
    input.register('x1', inp.LOCAL, inp.PROMPT)
    input.register('x2', inp.LOCAL, inp.PROMPT)
    input.register('y1', inp.LOCAL, inp.PROMPT)
    input.register('y2', inp.LOCAL, inp.PROMPT)
    input.register('h1', inp.LOCAL, inp.PROMPT)
    input.register('h2', inp.LOCAL, inp.PROMPT)
    input.register('fwmax', inp.LOCAL, inp.PROMPT)
    input.register('fwmin', inp.LOCAL, inp.PROMPT)
    input.register('angle', inp.LOCAL, inp.PROMPT)
    input.register('beta', inp.LOCAL, inp.PROMPT)
    input.register('fname', inp.LOCAL, inp.PROMPT)

    # Create a target field
    targs = hcam.Field()
    targs.add_random(ntarg, x1, x2, y1, y2, h1, h2, 
                     fwmax, fwmin, angle, beta)

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
