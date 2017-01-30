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
import numpy as np
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
        self._fwmin = fwmin
        self._fwmax = fwmax
        self._angle = angle
        self._beta = beta
        self._comp_abc()

    def _comp_abc(self):
        # compute hidden parameters _a, _b and _c used
        # to describe the ellipsoidal target shape. This
        # must be invoked any time beta, fwmax, fwmin, or
        # angle are changed.

        # First compute, alpha the value the quadratic form in x,y needs
        # to reach to drop the height by a factor of 2.
        alpha = 2**(1./self.beta)-1.

        # Then the cosine and sine of the angle
        ct = math.cos(math.radians(self.angle))
        st = math.sin(math.radians(self.angle))

        # _a, _b, _c appear in a*x**2 + b*y**2 + c*x*y which
        # is used instead of the (r/r0)**2 of the Moffat
        # function description.
        self._a = 4*alpha*((ct/self.fwmax)**2+(st/self.fwmin)**2)
        self._b = 4*alpha*((st/self.fwmax)**2+(ct/self.fwmin)**2)
        self._c = 8*alpha*ct*st*(1/self.fwmax**2-1/self.fwmin**2)

    @property
    def fwmax(self):
        return self._fwmax

    @fwmax.setter
    def fwmax(self, fwmax):
        self._fwmax = fwmax
        self._comp_abc()

    @property
    def fwmin(self):
        return self._fwmin

    @fwmin.setter
    def fwmin(self, fwmin):
        self._fwmin = fwmin
        self._comp_abc()

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, angle):
        self._angle = angle
        self._comp_abc()

    @property
    def beta(self):
        return self._beta

    @beta.setter
    def beta(self, beta):
        self._beta = beta
        self._comp_abc()

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

class Field(list):
    """Object to represent a star field as a list of :class:`Target`s. This is to help
    with simulation of star fields.
    """

    def __init__(self):
        super(Field,self).__init__()

    def add_random(self, ntarg, x1, x2, y1, y2, h1, h2, fwmax, fwmin, angle, beta):
        """Adds a random field of ntarg :class:`Target` objects scattered over the
        range x1, x2, y1, y2, with peak heights from h1 to h2 and fixed width
        parameters. The peak heights are chosen between the limits to have a
        distribution appropriate to a 3D Euclidean distribution of constant
        luminosity objects.

        Arguments::

           ntarg : (int)
              number of targets to add

           x1 : (float)
              left limit of region to add targets to

           x2 : (float)
              right limit of region to add targets to

           y1 : (float)
              lower limit of region to add targets to

           y2 : (float)
              upper limit of region to add targets to

           h1 : (float)
              lowest peak height

           h2 : (float)
              highest peak height

           fwmax : (float)
              FWHM along major axis of objects

           fwmin : (float)
              FWHM along minor axis of objects

           angle : (float)
              angle degrees clockwise from the X axis

           beta : (float)
              exponent that appears in Moffat functions

        """

        for nt in range(ntarg):
            xcen = np.random.uniform(x1, x2)
            ycen = np.random.uniform(y1, y2)
            height = 1./(1./math.sqrt(h1)-(1./math.sqrt(h1)-1./math.sqrt(h2))*
                         np.random.uniform())**2
            self.append(Target(xcen, ycen, height, fwmax, fwmin, angle, beta))

    def wjson(self, fname):
        """Writes a :class:`Field` to a file in json format. This is provided as a
        straightforward way to store all the info of a :class:`Field`.

        Arguments::

           fname : (string)
              file to write to

        """
        with open(fname, 'w') as fp:
            json.dump(self, fp, cls=TargetEncoder)

    @classmethod
    def rjson(cls, fname):
        """Creates a :class:`Field` from a json format file as saved by `wjson`.

        Arguments::

           fname : (string)
              the file to read to create the :class:`Field`
        """
        with open(fname) as fp:
            data = json.load(fp)
        field = cls()
        for dat in data:
            field.append(
                Target(dat['xcen'], dat['ycen'], dat['height'],
                       dat['fwmax'], dat['fwmin'], dat['angle'], dat['beta'])
            )
        return field

class TargetEncoder(json.JSONEncoder):
    """Provides a translation to json for :class:`Target` as required
    for writing :class:`Field`s to json
    """

    def default(self, obj):

        if isinstance(obj, Target):
            # Catch Targets
            return {
                '_comment': 'hipercam.Target',
                'xcen' : obj.xcen, 'ycen' : obj.ycen, 'height' : obj.height,
                'fwmax' : obj.fwmax, 'fwmin' : obj.fwmin, 'angle' : obj.angle,
                'beta' : obj.beta}

        # Default
        return json.JSONEncoder.default(self, obj)
