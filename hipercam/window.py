# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to represent sub-windows of a CCD and associated
functions.
"""

# Standard pre-amble from astropy
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.extern import six

import numpy as np
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

    def wfhead(self, nccd, nwin, fhead, nxccd, nyccd, NX, NY):
        """
        Writes window data into a FITS header. Assumes HiPERCAM type
        data in the sense that the amplifier is computed according to
        the window number. Non-HiPERCAM data will still work but the
        amplifier assignment will be wrong.

        Arguments::

          nccd : (int)
             CCD number, starting from 0

          nwin : (int)
             Window number, starting from 0

          fhead : (astropy.io.fits.Header)
             FITS header object. Will be modified on exit.

          nxccd : (int)
             X dimension of CCD

          nyccd : (int)
             Y dimension of CCD

          NX : (int)
             X dimension for iraf mosaic format

          NY : (int)
             Y dimension for iraf mosaic format

        Returns::

          fhead : (astropy.io.fits.Header)
             The modified FITS header object.
        """

        # The next line will likely need modifying
        namp = nwin % 4 + 1

        fhead['NCCD']    = (nccd+1,'CCD number')
        fhead['NWIN']    = (nwin+1,'Window number')
        fhead['LLX']     = (self.llx,  'X-coordinate of lower-left pixel')
        fhead['LLY']     = (self.lly,  'Y-coordinate of lower-left pixel')
        fhead['XBIN']    = (self.xbin, 'Binning factor in X')
        fhead['YBIN']    = (self.ybin, 'Binning factor in Y')

        fhead['CCDNAME'] = (str(nccd+1),'CCD identifier')
        fhead['AMPNAME'] = (str(namp),'Amplifier identifier')
        fhead['CCDSIZE'] = ('[1:' + str(nxccd) + ',1:' + str(nyccd) + ']', 'CCD size')

        # Flags to define the amplifier quadrant
        left   = nwin % 4 == 0 or nwin % 4 == 2
        bottom = nwin % 4 == 0 or nwin % 4 == 1

        # location within the NX by NY array of CCDs
        nx = nccd % NX
        ny = nccd // NY

        # Now we start on the iraf mosaic stuff, see
        # http://iraf.noao.edu/projects/ccdmosaic/imagedef/imagedef.html
        # for more detail.
        # Basically, 4 coordinate systems: CCD, Amplifier, Data, Detector
        #
        # CCD  -- basic unbinned 1,1 in lower-left corner
        # Amp  -- same, but measured with respect to amplifier
        # Data -- 

        # head
        ccdsec = ampsec = datasec = detsec = '['

        # x range
        ccdsec += str(self.llx) + ':' + str(self.llx+self.nx*self.xbin-1)
        if left:
            ampsec += str(self.llx) + ':' + str(self.llx+self.nx*self.xbin-1)
        else:
            ampsec += str(nxccd-self.llx+1) + ':' + str(nxccd-self.llx-self.nx*self.xbin+2)
        datasec += '1:' + str(self.nx)
        detsec  += str(self.llx+nx*nxccd) + ':' + str(self.llx+nx*nxccd+self.nx*self.xbin-1)

        # middle
        ccdsec += ','
        ampsec += ','
        datasec += ','
        detsec += ','

        # y range
        ccdsec += str(self.lly) + ':' + str(self.lly+self.nx*self.ybin-1)
        if bottom:
            ampsec += str(self.lly) + ':' + str(self.lly+self.ny*self.ybin-1)
        else:
            ampsec += str(nyccd-self.lly+1) + ':' + str(nyccd-self.lly-self.ny*self.ybin+2)
        datasec += '1:' + str(self.ny)
        detsec  += str(self.lly+ny*nyccd) + ':' + str(self.lly+ny*nyccd+self.ny*self.ybin-1)

        # tail
        ccdsec += ']'
        ampsec += ']'
        datasec += ']'
        detsec += ']'


        fhead['CCDSEC']   = ccdsec
        fhead['AMPSEC']   = ccdsec
        fhead['DATASEC']  = datasec
        fhead['DETSEC']   = detsec
        fhead['CCDSUM']   = str(self.xbin) + ' ' + str(self.ybin)

#        fhead['CTYPE1']   = 'Pixel'
#        fhead['CTYPE2']   = 'Pixel'
#        fhead['CRPIX1']   = 1
#        fhead['CRPIX2']   = 1
#        fhead['CRVAL1']   = xref
#        fhead['CRVAL2']   = yref
#        fhead['CDELT1']   = self.xbin
#        fhead['CDELT2']   = self.ybin

        if left:
            fhead['ATM1_1']   = 1.
            fhead['ATV1']     = 0.
        else:
            fhead['ATM1_1']   = -1.
            fhead['ATV1']     = nxccd+1

        if bottom:
            fhead['ATM2_2']   = 1.
            fhead['ATV2']     = 0.
        else:
            fhead['ATM2_2']   = -1.
            fhead['ATV2']     = nyccd+1

        fhead['LTM1_1']   = 1./self.xbin
        fhead['LTM2_2']   = 1./self.ybin
        fhead['LTV1']     = (1-self.llx)/self.xbin
        fhead['LTV2']     = (1-self.lly)/self.ybin
        fhead['DTM1_1']   = 1.
        fhead['DTM2_2']   = 1.
        fhead['DTV1']     = nx*nxccd
        fhead['DTV2']     = ny*nyccd

        return fhead

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

def overlap(win1, win2):
    """Returns True if :class:`Window's win1 and win2 overlap, False if they
    don't. Overlap means that they at least one unbinned pixel in common.
    """
    return win1.llx <=  win2.urx and win1.urx >= win2.llx and \
        win1.lly <=  win2.ury and win1.ury >= win2.lly


class Windata(Window):
    """
    Class representing a CCD window with data. Constructed from
    a :class:`Window` and a :class:`numpy.ndarray`

        >>> import numpy as np
        >>> from hipercam import Window, Windata
        >>> win = Window(12, 6, 100, 150, 2, 3)
        >>> data = np.ones((150,100))
        >>> wind = Windata(win,data)
    """
    def __init__(self, win, data=None):
        """
        Constructs a :class:`Windata`

        Arguments::

          win : (Window)
              the Window

          data : (numpy.ndarray)
              the data (2D). The dimension must match those in win unless
              data is None in which case a zero array of the correct size will
              be created.
        """
        if data is None:
            self.data = np.zeros((win.ny,win.nx))
        else:
            # Run a couple of checks
            if data.ndim != 2:
                raise HipercamError('Windata.__init__: data must be 2D. Found {0:d}'.format(data.ndim))
            ny, nx = data.shape
            if nx != win.nx or ny != win.ny:
                raise HipercamError('Windata.__init__: win & data dimensions conflict. NX: {0:d} vs {1:d}, NY: {2:d} vs {3:d}'.format(win.nx,nx,win.ny,ny))

            self.data = data

        super(Windata, self).__init__(win.llx, win.lly,
                                      win.nx, win.ny, win.xbin, win.ybin)

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

    def __iadd__(self, other):
        """+- in-place addition operator. 'other' can be another Windata,
        an numpy.ndarray or a constant. Window parameters unchanged.
        """
        if isinstance(other, Windata):
            self.data += other.data
        else:
            self.data += other
        return self

    def __repr__(self):
        return 'Windata(win=' + super(Windata,self).__repr__() + \
                              ', data=' + repr(self.data) + ')'
