# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import hipercam as hcam

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

    print('Creating 5 Windatas')
    wind1 = hcam.Windata(win1)
    wind2 = hcam.Windata(win2)
    wind3 = hcam.Windata(win3)
    wind4 = hcam.Windata(win4)
    wind5 = hcam.Windata(win5)

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
