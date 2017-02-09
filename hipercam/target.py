# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to represent astronomical object fields for
simulating data.
"""

import math
import json

import numpy as np

from .core import *

__all__ = ('Target','Field')

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

    def __init__(self, xcen, ycen, height, fwhm1, fwhm2, angle, beta, fmin):
        """
        Constructor. Arguments::

        xcen : (float)
            X position of target (unbinned pixels)

        ycen : (float)
            X position of target (unbinned pixels)

        height : (float)
            height of central peak (counts)

        fwhm1 : (float)
            FWHM along one axis of the ellipse (unbinned pixels)

        fwhm2 : (float)
            FWHM along other axis of the ellipse (unbinned pixels)

        angle : (float)
            angle of axis 1 of ellipse, clockwise from X-axis (degrees)

        beta : (float)
            exponent that controls how quickly the wings fall. beta=1
            is a Lorentzian, beta=infinity is a Gausssin.

        fmin : (float)
            in order to define the region affected by a target, this is the
            minimum level to which the profile will be computed. It must be
            less than height in terms of absolute value.
        """

        if (height > 0. and height < fmin) or \
           (height < 0. and height > fmin):
            raise ValueError(
                'hipercam.Target.__init__: fmin vs height conflict: {0:f} vs {1:f}'.format(fmin,height))

        self.xcen = xcen
        self.ycen = ycen
        self.height = height
        self.fmin = fmin
        self._fwhm1 = fwhm1
        self._fwhm2 = fwhm2
        self._angle = angle
        self._beta = beta
        self._comp_abc()

    def copy(self, memo=None):
        """Returns a copy of the :class:`Target`"""
        return Target(self.xcen,self.ycen,self.height, self.fwhm1,
                      self.fwhm2, self.angle, self.beta, self.fmin)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy(memo)

    def _comp_abc(self):
        # compute hidden parameters _a, _b and _c used
        # to describe the ellipsoidal target shape. This
        # must be invoked any time beta, fwhm1, fwhm2, or
        # angle are changed.

        # First compute, alpha the value the quadratic form in x,y needs
        # to reach to drop the height by a factor of 2.
        alpha = 2**(1./self.beta)-1.

        # Then the cosine and sine of the angle
        ct = math.cos(math.radians(self.angle))
        st = math.sin(math.radians(self.angle))

        # _a, _b, _c appear in a*x**2 + b*y**2 + 2*c*x*y which
        # is used instead of the (r/r0)**2 of the Moffat
        # function description.
        self._a = 4*alpha*((ct/self.fwhm1)**2+(st/self.fwhm2)**2)
        self._b = 4*alpha*((st/self.fwhm1)**2+(ct/self.fwhm2)**2)
        self._c = 4*alpha*ct*st*(1/self.fwhm1**2-1/self.fwhm2**2)

        # Maximum value of the rsq parameter (see __call__) to
        # calculate to. This defines the region to which the profile
        # needs to be calculated which will be stored for fast lookup
        rsqmax = (self.height/self.fmin)**(1./self.beta)-1
        det = self._a*self._b-self._c**2
        xmax = np.sqrt(rsqmax*self._b/det)
        ymax = np.sqrt(rsqmax*self._a/det)
        self._x1 = self.xcen - xmax
        self._x2 = self.xcen + xmax
        self._y1 = self.ycen - ymax
        self._y2 = self.ycen + ymax

    @property
    def fwhm1(self):
        return self._fwhm1

    @fwhm1.setter
    def fwhm1(self, fwhm1):
        self._fwhm1 = fwhm1
        self._comp_abc()

    @property
    def fwhm2(self):
        return self._fwhm2

    @fwhm2.setter
    def fwhm2(self, fwhm2):
        self._fwhm2 = fwhm2
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

    def offset(self, dx, dy):
        """This generates a new :class:`Target` by offsetting the position of the
:class:`Target`

        Arguments::

           dx : (float)
             X offset to add to the positions of all targets [unbinned pixels]

           dy : (float)
             Y offset to add to the positions of all targets [unbinned pixels]

        Returns:: a shifted version of the :class:`Target`

        """
        targ = self.copy()
        targ.xcen += dx
        targ.ycen += dy
        return targ

    def __call__(self, x, y, f=None, scale=1.):
        """Computes data representing the :class:`Target` at
        position x, y (which can be arrays of positions). If f=None,
        the result is simply returned. Otherwise if f is an array
        matching x and y, it will be modified on exit as well as
        being returned. This can save some overheads.

        Example code::

          >> x = np.linspace(1,100,100)
          >> y = np.linspace(1,100,100)
          >> X,Y = np.meshgrid(x,y)
          >> I = np.zeros_like(X)
          >> targ = Target(100., 150., 100., 5., 3.5, 30., 4.)
          >> targ(X,Y,I)

        This creates a target, generates 2D arrays of X and Y positions
        then calculates the value of the target at every pixel of these
        arrays, adding it to I.

        Arguments::

          x : (array / float)
            x position or positions to calculate the profile for

          y : (array / float)
            y position or positions to calculate the profile for

          f : (array / None)
            if x and y are arrays, and f a matching array, then the target will
            be added to it. If f=None, the result is returned. f is ignored
            if x and y are floats.

          scale : (float)
            factor to scale how much is added in.
        """

        if isinstance(x,np.ndarray) and isinstance(y,np.ndarray):

            if f is None:
                f = np.zeros_like(x)

            ok = (x > self._x1) & (x < self._x2) & \
                 (y > self._y1) & (y < self._y2)

            if ok.any():
                xd = x[ok]-self.xcen
                yd = y[ok]-self.ycen
                rsq = self._a*xd**2 + self._b*yd**2 + (2*self._c)*xd*yd
                f[ok] += scale*self.height/(1+rsq)**self.beta

            return f

        elif ok:
            # x, y both floats.
            rsq = self._a*xd**2 + self._b*yd**2 + (2*self._c)*xd*yd
            return scale*self.height/(1+rsq)**self.beta
        else:
            return 0.

    def __repr__(self):
        return 'Target(xcen=' + repr(self.xcen) + ', ycen=' + repr(self.ycen) + \
               ', height=' + repr(self.height) + ', fwhm2=' + repr(self.fwhm2) + \
               ', fwhm1=' + repr(self.fwhm1) + ', angle=' + repr(self.angle) + \
               ', beta=' + repr(self.beta) + ')'

class Field(list):
    """Object to represent a star field as a list of :class:`Target`s. This is to help
    with simulation of star fields.
    """

    def __init__(self):
        super(Field,self).__init__()

    def add_random(self, ntarg, x1, x2, y1, y2, h1, h2, angle1, angle2,
                   fwhm1, fwhm2, beta, fmin):
        """Adds a random field of ntarg :class:`Target` objects scattered over the
        range x1, x2, y1, y2, with peak heights from h1 to h2, random angles
        from angle1 to angle2, but with fixed width parameters. The peak
        heights are chosen between the limits to have a distribution
        appropriate to a 3D Euclidean distribution of constant luminosity
        objects. The distributions over x, y and angle are uniform.

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

           angle1 : (float)
              lowest angle degrees clockwise from the X axis

           angle2 : (float)
              highest angle degrees clockwise from the X axis

           fwhm1 : (float)
              FWHM along major axis of objects

           fwhm2 : (float)
              FWHM along minor axis of objects

           beta : (float)
              exponent that appears in Moffat functions

           fmin : (float)

              parameter used to restrict region affected by a target: it's the
              minimum count level to which the profile will be computed.

        """

        for nt in range(ntarg):
            xcen = np.random.uniform(x1, x2)
            ycen = np.random.uniform(y1, y2)
            height = 1./(1./math.sqrt(h1)-(1./math.sqrt(h1)-1./math.sqrt(h2))*
                         np.random.uniform())**2
            angle = np.random.uniform(angle1, angle2)
            self.append(Target(xcen, ycen, height, fwhm1, fwhm2, angle,
                               beta, fmin))

    def modify(self, transform, fscale):
        """This generates a new :class:`Field` by applying `transform` to the x,y
        locations of each target to change them. It also scales the targets'
        peak heights by a constant factor fscale.

        Arguments::

           transform : (callable)
               Should have signature transform(x,y) and when given the centre
               x,y of a target, it should return the change in x,y needed to
               get to the new position which will be fed to Target.offset

           fscale : (float)
               Scaling factor  It should also have an attribute `fscale`
               that is a scaling factor that will be applied to the targt peak
               heights.

        See `hipercam.command_line.makedata` for an example.

        """
        field = Field()
        for target in self:
            dx,dy = transform(target.xcen, target.ycen)
            targ = target.offset(dx,dy)
            targ.height *= fscale
            field.append(targ)
        return field

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
                       dat['fwhm1'], dat['fwhm2'], dat['angle'], dat['beta'])
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
                'fwhm1' : obj.fwhm1, 'fwhm2' : obj.fwhm2, 'angle' : obj.angle,
                'beta' : obj.beta}

        # Default
        return json.JSONEncoder.default(self, obj)
