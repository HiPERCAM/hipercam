# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to represent sub-windows of a CCD and associated
functions.
"""

import json
import numpy as np
from astropy.io import fits
from .core import *

class Window:
    """
    Class representing a window of a CCD. This represents an arbitrary
    rectangular region of binned pixels. The lower-left pixel of the CCD
    is assumed to have coordinates (x,y) = (1,1). :class:`Window` dimensions are
    in binned pixels.

        >>> from hipercam import Window
        >>> win = Window(12, 6, 100, 150, 2, 3)
        >>> print(win)
        Window(12, 6, 100, 150, 2, 3)

    """

    def __init__(self, llx, lly, nx, ny, xbin, ybin):
        """
        Constructor. Arguments::

        llx : (int)
            X position of lower-left pixel of window (unbinned pixels)

        lly : (int)
            Y position of lower-left pixel of window (unbinned pixels)

        nx : (int)
            X dimension of window, binned pixels

        ny : (int)
            Y dimension of window, binned pixels

        xbin : (int)
            Binning factor in X

        ybin : (int)
            Binning factor in Y
        """

        self.llx   = llx
        self.lly   = lly
        self.nx    = nx
        self.ny    = ny
        self.xbin  = xbin
        self.ybin  = ybin

    def __repr__(self):
        return 'Window(llx=' + repr(self.llx) + ', lly=' + repr(self.lly) + \
            ', nx=' + repr(self.nx) + ', ny=' + repr(self.ny) + \
            ', xbin=' + repr(self.xbin) + ', ybin=' + repr(self.ybin) + ')'


    @property
    def urx(self):
        """
        Returns (unbinned) X pixel at upper-right of :class:`Window`
        """
        return self.llx-1+self.nx*self.xbin

    @property
    def ury(self):
        """
        Returns (unbinned) Y pixel at upper-right of :class:`Window`
        """
        return self.lly-1+self.ny*self.ybin

    def extent(self):
        """
        Returns (left,right,bottom,top) boundaries of :class:`Window`
        """
        return (self.llx-0.5,self.urx+0.5,self.lly-0.5,self.ury+0.5)

    def outside(self, win):
        """
        Returns True if `self` contains `win` in such a way
        that it could be cut down to it. This implies that
        even if binned, its pixels are "in step", or that its
        bin factors are divisors of those of `win`.

        See also :func:`inside`.
        """
        return win.xbin % self.xbin == 0 and win.ybin % self.ybin == 0 and \
            self.llx <= win.llx and self.urx >= win.urx and \
            self.lly <= win.lly and self.ury >= win.ury and \
            (win.llx-self.llx) % self.xbin == 0 and \
            (win.lly-self.lly) % self.ybin == 0

    def inside(self, win):
        """
        Returns True if `win` contains `self` in such a way
        that it could be cut down to it. This implies that
        even if binned, its pixels are "in step", or that its
        bin factors are divisors of those of `self`.

        See also :func:`outside`.
        """
        return self.xbin % win.xbin == 0 and self.ybin % win.ybin == 0 and \
            win.llx <= self.llx and win.urx >= self.urx and \
            win.lly <= self.lly and win.ury >= self.ury and \
            (self.llx-win.llx) % win.xbin == 0 and \
            (self.lly-win.lly) % win.ybin == 0

    def clash(self, win):
        """Raises a HipercamError exception if two :class: `Window`s are considered to
        'clash'.  In this case this means if they have any pixels in common.
        This method is required for :class: `Window`s to be collected into
        :class:`Group`s

        """
        if self.llx <=  win.urx and self.urx >= win.llx and \
           self.lly <=  win.ury and self.ury >= win.lly:
            raise HipercamError('Window.clash: self = ' + str(self) + ' clashes with win = ' + str(win))

    def xy(self):
        """Returns two 2D arrays containing the x and y values at the centre
        of each pixel defined by the :class:`Window`. See numpy.meshgrid to
        see what this means.
        """
        # generate 1D x and y arrays along the edges
        x1 = self.llx+self.xbin/2-0.5
        x2 = x1 + self.xbin*(self.nx-1)
        x = np.linspace(x1,x2,self.nx)

        y1 = self.lly+self.ybin/2-0.5
        y2 = y1 + self.ybin*(self.ny-1)
        y = np.linspace(y1,y2,self.ny)

        return np.meshgrid(x,y)

    def matches(self, win):
        """Tests that the :class:`Window` matches another. If all OK, returns None,
        otherwise raises a HipercamError reporting the two
        :class:`Window`s. See also `__eq__`

        """
        if self != win:
            raise HipercamError(
                'hipercam.Window.matches: self = ' + str(self) +
                ' clashes with win = ' + str(win))

    def __eq__(self, win):
        """
        Defines equality. Two :class:`Window`s are equal if they match exactly
        (same lower left corner, dimensions and binning factors)
        """
        return self.llx == win.llx and  self.lly == win.lly and \
            self.nx == win.nx and self.ny == win.ny and \
            self.xbin == win.xbin and self.ybin == win.ybin

    def __ne__(self, win):
        """Defines equality. Two :class:`Window`s are equal if they match
        exactly (same lower left corner, dimensions and binning factors)

        """
        return not (self.llx == win)

class WindowEncoder (json.JSONEncoder):
    """
    Provides a default that can be used to write Window
    objects to json files
    """
    def default(self, obj):
        if instance(obj, Window):
            obj = {'_comment': 'This is a hipercam.Window JSON file.',
                   'llx' : self.llx, 'lly' : self.lly,
                   'nx' : self.nx, 'ny' : self.ny,
                   'xbin' : self.xbin, 'ybin' : self.ybin}
        else:
            super(WindowEncoder, self).default(obj)

class Windat(Window):
    """
    Class representing a CCD window with data. Constructed from
    a :class:`Window` and a :class:`numpy.ndarray` which is stored
    in an attribute called `data`.

        >>> import numpy as np
        >>> from hipercam import Window, Windat
        >>> win = Window(12, 6, 100, 150, 2, 3)
        >>> data = np.ones((150,100))
        >>> wind = Windat(win,data)
        >>> wind += 0.5
        >>> wind *= 2

    Note that it is the user's reponsibility to make sure that the data
    and Window remain compatible.
    """

    def __init__(self, win, data=None):
        """
        Constructs a :class:`Windat`

        Arguments::

          win : (Window)
              the Window

          data : (numpy.ndarray)
              the data (2D). The dimension must match those in win unless
              data is None in which case a zero array of the correct size
              will be created.
        """
        if data is None:
            self.data = np.zeros((win.ny,win.nx))
        else:
            # Run a couple of checks
            if data.ndim != 2:
                raise HipercamError(
                    'Windat.__init__: data must be 2D. Found {0:d}'.format(data.ndim))
            ny, nx = data.shape
            if nx != win.nx or ny != win.ny:
                raise HipercamError(
                    'Windat.__init__: win vs data dimension conflict. NX: {0:d} vs {1:d}, NY: {2:d} vs {3:d}'.format(win.nx,nx,win.ny,ny))

            self.data = data

        super(Windat, self).__init__(win.llx, win.lly,
                                      win.nx, win.ny, win.xbin, win.ybin)

    @classmethod
    def rhdu(cls, hdu):
        """Constructs a :class:`Windat` from an ImageHdu. Requires
        header parameters 'LLX', 'LLY', 'XBIN' and 'YBIN' to be defined.
        """
        head = hdu.header
        data = hdu.data

        # construct Window
        llx = head['LLX']
        lly = head['LLY']
        xbin = head['XBIN']
        ybin = head['YBIN']
        ny, nx = data.shape
        win = Window(llx, lly, nx, ny, xbin, ybin)

        return cls(win, data)

    @property
    def size(self):
        """
        Number of pixels
        """
        return self.data.size

    @property
    def win(self):
        """A copy of the :class:`Window` underlying the :class:`Windat`"""
        return Window(self.llx, self.lly, self.nx, self.ny, self.xbin, self.ybin)

    def add_noise(self, readout, gain):
        """Adds noise to a :class:`Windat` according to a variance
        calculated from V = readout**2 + counts/gain.
        Arguments::

          readout : (float)
              RMS readout noise in counts

          gain : (float)
              Gain in electrons per count.
        """
        sig = np.sqrt(readout**2+self.data/gain)
        self.data += np.random.normal(scale=sig)

    def min(self):
        """
        Returns the minimum value of the :class:`Windat`.
        """
        return self.data.min()

    def max(self):
        """
        Returns the maximum value of the :class:`Windat`.
        """
        return self.data.max()

    def mean(self):
        """
        Returns the mean value of the :class:`Windat`.
        """
        return self.data.mean()

    def percentile(self, q):
        """
        Computes percentile(s) of a :class:`Windat`.

        Arguments::

        q : (float or sequence of floats)
          Percentile(s) to use, in range [0,100]
        """
        return np.percentile(self.data, q)

    def whdu(self, fhead=None):
        """
        Writes the :class:`Windat` to an :class:`astropy.io.fits.ImageHDU`

        Arguments::

          fhead : (astropy.io.fits.Header)
             An initial FITS header object. The location of the lower-left pixel
             and the binning factors will be added to it, so it will be modified
             on exit. 

        Returns::

          hdu : (astropy.io.fits.ImageHDU)
             The HDU
        """

        if fhead is None:
            fhead = fits.Header()

        fhead['LLX'] = (self.llx,  'X-coordinate of lower-left pixel')
        fhead['LLY'] = (self.lly,  'Y-coordinate of lower-left pixel')
        fhead['XBIN'] = (self.xbin, 'Binning factor in X')
        fhead['YBIN'] = (self.ybin, 'Binning factor in Y')

        return fits.ImageHDU(self.data, fhead)

    def add_fxy(self, funcs):
        """Routine to add in the results of evaluating a function
        or a list of functions of x & y to the :class:`Windat`.
        Each function must take 2D arrays of x and y values
        for each pixel of the :class:`Windat` and return an
        array of values for each point. If you have lots of things
        to add, this routine saves some overhead by only generating
        the x,y pixel arrays once at the start.

        Arguments::

          funcs : (a callable or a list of callables)
              the callables must have signature arr = func(x,y) where x
              and y are 2D arrays containing the x and y values for every
              pixel in the :class:`Windat`
        """
        # generate X,Y arrays
        X,Y = self.xy()
        try:
            for func in funcs:
                self.data += func(X,Y)
        except TypeError:
            self.data += funcs(X,Y)

    def __iadd__(self, other):
        """+= in-place addition operator. 'other' can be another :class:`Windat`,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windat):
            self.data += other.data
        else:
            self.data += other
        return self

    def __isub__(self, other):
        """-= in-place subtraction operator. 'other' can be another Windat,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windat):
            self.data -= other.data
        else:
            self.data -= other
        return self

    def __imul__(self, other):
        """*= in-place multiplication operator. 'other' can be another Windat,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windat):
            self.data *= other.data
        else:
            self.data *= other
        return self

    def __itruediv__(self, other):
        """/= in-place division operator. 'other' can be another Windat,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windat):
            self.data /= other.data
        else:
            self.data /= other
        return self

    def __add__(self, other):
        """Adds `other` to a :class:`Windat` as `wind + other`.  Here `other` can be a
        compatible :class:`Windat` (identical window) or any object that can
        be added to a :class:`numpy.ndarray`, e.g. a float, another matching
        array, etc.

        """

        if isinstance(other, Windat):
            # test compatibility between the windows (raises an exception)
            self.matches(other)
            return Windat(self.win, self.data + other.data)
        else:
           return Windat(self.win, self.data + other)

    def __radd__(self, other):
        """Adds `other` to a :class:`Windat` as `other + wind`.  Here `other` isany
        object that can be added to a :class:`numpy.ndarray`, e.g. a float,
        another matching array, etc.

        """
        return Windat(self.win, self.data + other)

    def __sub__(self, other):
        """Subtracts `other` from a :class:`Windat` as `wind - other`.  Here `other`
        can be a compatible :class:`Windat` (identical window) or any object
        that can be subtracted from a :class:`numpy.ndarray`, e.g. a float, another matching
        array, etc.

        """

        if isinstance(other, Windat):
            # test compatibility between the windows (raises an exception)
            self.matches(other)
            return Windat(self.win, self.data - other.data)
        else:
            return Windat(self.win, self.data - other)

    def __rsub__(self, other):
        """Subtracts a :class:`Windat` from `other` as `other - wind`.  Here `other`
        is any object that can have a :class:`numpy.ndarray` subtracted from it, e.g. a
        float, another matching array, etc.

        """
        return Windat(self.win, other - self.data)

    def __mul__(self, other):
        """Multiplies a :class:`Windat` by `other` as `wind * other`.  Here `other`
        can be a compatible :class:`Windat` (identical window) or any object
        that can multiply a :class:`numpy.ndarray`, e.g. a float, another
        matching array, etc.

        """

        if isinstance(other, Windat):
            # test compatibility between the windows (raises an exception)
            self.matches(other)
            return Windat(self.win, self.data * other.data)
        else:
            return Windat(self.win, self.data * other)

    def __rmul__(self, other):
        """Multiplies a :class:`Windat` by `other` as `other * wind`.  Here `other` is
        any object that can multiply a :class:`numpy.ndarray`, e.g. a float,
        another matching array, etc.

        """
        return Windat(self.win, other * self.data)

    def __div__(self, other):
        """Divides a :class:`Windat` by `other` as `wind / other`.  Here `other`
        can be a compatible :class:`Windat` (identical window) or any object
        that can divide into a :class:`numpy.ndarray`, e.g. a float, another
        matching array, etc.

        """

        if isinstance(other, Windat):
            # test compatibility between the windows (raises an exception)
            self.matches(other)
            return Windat(self.win, self.data / other.data)
        else:
            return Windat(self.win, self.data / other)

    def __rdiv__(self, other):
        """Divides `other` by a :class:`Windat` as `other / wind`.  Here `other` is
        any object that can be divided by a :class:`numpy.ndarray`, e.g. a
        float, another matching array, etc.

        """
        return Windat(self.win, other / self.data)

    def __repr__(self):
        return 'Windat(' + super(Windat,self).__repr__() + \
            ', data=' + repr(self.data) + ')'




