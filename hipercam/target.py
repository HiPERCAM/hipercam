# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to represent astronomical object fields for
simulating data.
"""

# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import math
import json
from .core import *

class Target(object):
    """
    Callable object to represent an astronomical target in terms of its appearance
    in a CCD image. The representation chosen is 2D "Moffat" function, which
    takes the form h/(1+(r/r0)**2)**beta, where 'h' is the central height,
    'r' is the distance from the centre, 'r0' is a scale, and 'beta' is an
    exponent which for large values makes the profile Gaussian.

    To allow for trailing, poor focus and to roughly simulate galaxies,
    :class:`Target` objects can be elliptical in shape.
    """

    def __init__(self, xcen, ycen, height, fwmax, fwmin, angle, beta):
        """
        Constructor. Arguments::

        xcen : (float)
            X position of target (unbinned pixels)

        ycen : (float)
            X position of target (unbinned pixels)

        height : (float)
            height of central peak in counts

        fwmax : (float)
            FWHM, unbinned pixel, along minor axis of ellipse

        fwmin : (float)
            FWHM, unbinned pixel, along minor axis of ellipse

        angle : (float)
            angle, degrees, clockwise from X-axis.

        beta : (float)
            exponent that controls how quickly the wings fall. beta=1
            is a Lorentzian, beta=infinity is a Gausssin.
        """

        self.xcen = xcen
        self.ycen = ycen
        self.height = height
        self.fwmin = fwmin
        self.fwmax = fwmax
        self.angle = angle
        self.beta = beta

        # compute parameters a / b / c used
        # to describe the ellipsoid

        # alpha is the value the quadratic form in x,y needs
        # to reach to drop the height by a factor of 2.
        alpha = 2**(1./beta)-1.

        # cos(theta), sin(theta)
        ct = math.cos(math.radians(angle))
        st = math.sin(math.radians(angle))

        # a, b, c appear in a*x**2 + b*y**2 + c*x*y which
        # is used instead of the (r/r0)**2 of the Moffat
        # function description.
        self._a = 4*alpha*((ct/fwmax)**2+(st/fwmin)**2)
        self._b = 4*alpha*((st/fwmax)**2+(ct/fwmin)**2)
        self._c = 8*alpha*ct*st*(1/fwmax**2-1/fwmin**2)

    def __call__(self, x, y):
        """Computes data representing the :class:`Target` at
        position x, y (which can be arrays of positions)

        Example code::

          >> x = np.linspace(1,100,100)
          >> y = np.linspace(1,100,100)
          >> X,Y = np.meshgrid(x,y)
          >> targ = Target(100., 150., 100., 5., 3.5, 30., 4.)
          >> img = targ(x,y)

        This creates a target, generates 2D arrays of X and Y positions
        then calculates the value of the target at every pixel of these
        arrays.

        """
        xd = x-self.xcen
        yd = y-self.ycen
        rsq = self._a*xd**2 + self._b*yd**2 + self._c*xd*yd
        return self.height/(1+rsq)**self.beta

    def __repr__(self):
        return 'Target(xcen=' + repr(self.xcen) + ', ycen=' + repr(self.ycen) + \
               ', height=' + repr(self.height) + ', fwmin=' + repr(self.fwmin) + \
               ', fwmax=' + repr(self.fwmax) + ', angle=' + repr(self.angle) + \
               ', beta=' + repr(self.beta) + ')'


