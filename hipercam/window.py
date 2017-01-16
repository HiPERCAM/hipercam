# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to represent sub-windows of a CCD and associated
functions.
"""

# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *


import json
import numpy as np
from astropy.io import fits
from .core import *

class Window(object):
    """
    Class representing a window of a CCD. This represents an arbitrary
    rectangular region of binned pixels. The lower-left pixel of the CCD
    is assumed to have coordinates (x,y) = (1,1). :class:`Window` dimensions are
    in binned pixels.

        >>> from hipercam import Window
        >>> win = Window(12, 6, 100, 150, 2, 3)
        >>> print(win)
        Window(llx=12, lly=6, nx=100, ny=150, xbin=2, ybin=3)

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
        """Returns True if two :class: `Window`s are considered to 'clash'.
        In this case this means if they have any pixels in common.
        This method allows :class: `Window`s to be collected in
        :class:`Group`s
        """
        return self.llx <=  win.urx and self.urx >= win.llx and \
            self.lly <=  win.ury and self.ury >= win.lly


    def wjson(self, file):
        """
        Writes a Window to a json file which gives a fairly easily
        read and editable form of file. "file" is either a file name
        of a file object.
        """
        print('wjson not implemented yet')

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

class Windata(Window):
    """
    Class representing a CCD window with data. Constructed from
    a :class:`Window` and a :class:`numpy.ndarray` which is stored
    in an attribute called `data`.

        >>> import numpy as np
        >>> from hipercam import Window, Windata
        >>> win = Window(12, 6, 100, 150, 2, 3)
        >>> data = np.ones((150,100))
        >>> wind = Windata(win,data)
        >>> wind += 0.5
        >>> wind *= 2

    Note that it is the user's reponsibility to make sure that the data
    and Window remain compatible.
    """

    def __init__(self, win, data=None):
        """
        Constructs a :class:`Windata`

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
                    'Windata.__init__: data must be 2D. Found {0:d}'.format(data.ndim))
            ny, nx = data.shape
            if nx != win.nx or ny != win.ny:
                raise HipercamError(
                    'Windata.__init__: win & data dimensions conflict. NX: {0:d} vs {1:d}, NY: {2:d} vs {3:d}'.format(win.nx,nx,win.ny,ny))

            self.data = data

        super(Windata, self).__init__(win.llx, win.lly,
                                      win.nx, win.ny, win.xbin, win.ybin)

    @property
    def size(self):
        """
        Number of pixels
        """
        return self.data.size

    def add_stars(self, targs, sxx, syy, rxy):
        """Routine to add gaussians to a :class:`Windata`.
        Arguments::

          targs : (array of 3 element arrays)
              array of (xc,yc,h) values marking the central position
              and height of each gaussian to add

          sxx : (float)
              sigma in X direction

          syy : (float)
              sigma in Y direction

          rxy : (float)
              X-Y correlation coefficient
        """
        x = np.linspace(self.llx,self.urx,self.nx)
        y = np.linspace(self.lly,self.ury,self.ny)
        X,Y = np.meshgrid(x,y)
        V = np.matrix([[sxx*sxx,sxx*syy*rxy],[sxx*syy*rxy,syy*syy]])
        A = np.linalg.inv(V)
        for xc,yc,h in targs:
            xd = X-xc
            yd = Y-yc
            dsq = A[0,0]*xd**2+(A[1,0]+A[0,1])*xd*yd+A[1,1]*yd**2
            self.data += h*np.exp(-dsq/2.)

    def add_noise(self, readout, gain):
        """Adds noise to a :class:`Windata` according to a variance
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
        Returns the minimum value of the :class:`Windata`.
        """
        return self.data.min()

    def max(self):
        """
        Returns the maximum value of the :class:`Windata`.
        """
        return self.data.max()

    def mean(self):
        """
        Returns the mean value of the :class:`Windata`.
        """
        return self.data.mean()

    def percentile(self, q):
        """
        Computes percentile(s) of a :class:`Windata`.

        Arguments::

        q : (float or sequence of floats)
          Percentile(s) to use, in range [0,100]
        """
        return np.percentile(self.data, q)

    def whdu(self, fhead=None):
        """
        Writes the Windata to an astropy.io.fits.ImageHDU

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

    def __iadd__(self, other):
        """+= in-place addition operator. 'other' can be another Windata,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windata):
            self.data += other.data
        else:
            self.data += other
        return self

    def __isub__(self, other):
        """-= in-place subtraction operator. 'other' can be another Windata,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windata):
            self.data -= other.data
        else:
            self.data -= other
        return self

    def __imul__(self, other):
        """*= in-place multiplication operator. 'other' can be another Windata,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windata):
            self.data *= other.data
        else:
            self.data *= other
        return self

    def __itruediv__(self, other):
        """/= in-place division operator. 'other' can be another Windata,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windata):
            self.data /= other.data
        else:
            self.data /= other
        return self

    def __repr__(self):
        return 'Windata(win=' + super(Windata,self).__repr__() + \
                              ', data=' + repr(self.data) + ')'
