import sys
import os

import numpy as np

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ = ['stats',]

#################################################################
#
# stats -- lists basic stats of each window of a multi-CCD image.
#
#################################################################

def stats(args=None):
    """Lists basic stats of a multi-CCD image, i.e. the minimum, maximum,
    mean, median and standard deviation of each window of each CCD. The output
    format can be altered to suit preference.

    Arguments::

      input  : (string)
         name of the MCCD file

      format : (string) [hidden, defaults to 9.3f]
         C-style format code as used in Python format statements for output of
         the numerical values. e.g. '300.00' is '6.2f' (6 characters toal, 2 after
         the decimal point), '1.22e24' is '.2e' (as many characters as needed, 2
         after the decimal point)

    """

    if args is None:
        args = sys.argv[1:]

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', 'stats', args) as cl:

        # register parameters
        cl.register('input', Cline.LOCAL, Cline.PROMPT)
        cl.register('format', Cline.LOCAL, Cline.HIDE)

        # get inputs
        frame = cl.get_value('input', 'frame to lists stats of',
                             cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(frame)

        cl.set_default('format','9.3f')
        format = cl.get_value('format', 'output format for numbers', '9.3f')

    for cnam, ccd in mccd.items():
        for wnam, wind in ccd.items():
            print(
                'CCD {0:s}, window {1:s}: min = {3:{2:s}}, max = {4:{2:s}}, mean = {5:{2:s}}, median = {6:{2:s}}, std = {7:{2:s}}'.format(
                    cnam, wnam, format, wind.min(), wind.max(), wind.mean(), wind.median(), wind.std()
                    )
                )

############################################################################
#
# From this point on come helper methods and classes that are not externally
# visible
#
############################################################################

class Dust:
    """This to generate gaussian dust specks on a flat field using
    the add_fxy method of Windats"""

    def __init__(self, xcen, ycen, rms, depth):
        self.xcen = xcen
        self.ycen = ycen
        self.rms = rms
        self.depth = depth

    def __call__(self, x, y, f):
        ok = (x > self.xcen-5.*self.rms) & (x < self.xcen+5.*self.rms) & \
             (y > self.ycen-5.*self.rms) & (y < self.ycen+5.*self.rms)

        if ok.any():
            rsq = (x[ok]-self.xcen)**2 + (y[ok]-self.ycen)**2
            f[ok] -= self.depth*np.exp(-rsq/(2.*self.rms**2))

class Transform:
    """Field transformation class. This is a callable that can be sent to the
    :class:`Field` method `modify`.
    """

    def __init__(self, nxtot, nytot, xcen, ycen, angle, scale, xoff, yoff):
        self.nxtot = nxtot
        self.nytot = nytot
        self.xcen = xcen
        self.ycen = ycen
        self.angle = angle
        self.scale = scale
        self.xoff = xoff
        self.yoff = yoff

    def __call__(self, x, y):
        # Rotator centre
        xc = (self.nxtot+1)/2. + self.xcen
        yc = (self.nytot+1)/2. + self.ycen

        # Rotate around CCD centre
        cosa = np.cos(np.radians(self.angle))
        sina = np.sin(np.radians(self.angle))
        xcen =  cosa*(x-xc) + sina*(y-yc)
        ycen = -sina*(x-xc) + cosa*(y-yc)

        # Change image scale
        xcen *= self.scale
        ycen *= self.scale

        # Apply offset
        xcen += xc + self.xoff
        ycen += yc + self.yoff

        # Return change in position
        return (xcen-x,ycen-y)

# globals and a function used in the multiprocessing used in 'makedata'. The globals
# are to reduce the need for pickling / unpickling large bits of data
_gframe = None
_gfield = None

def worker(arg):
    global _gframe, _gfield
    cnam, ndiv = arg
    ccd = _gframe[cnam]
    field = _gfield[cnam]
    for wind in ccd.values():
        field.add(wind, ndiv)
    return ccd
