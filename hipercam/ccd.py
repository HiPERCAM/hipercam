# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Class to represent a CCD and a multi-CCD
"""

# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import numpy as np
from astropy.io import fits
from .core import *
from .group import *
from .window import *

class CCD(Group):
    """
    Class representing a CCD as a Group of Windata objects
    plus a FITS header
    """
    def __init__(self, winds, nxtot, nytot, head=None):
        """
        Constructs a :class:`CCD`

        Arguments::

          winds : (dict)
              dictionary of Windata objects. The dictionary keys should be integers.

          nxtot : (int)
              Unbinned X-dimension of CCD

          nytot : (int)
              Unbinned Y-dimension of CCD

          head : (astropy.io.fits.Header)
              a header which will be written along with the first of the Windata
              sub-images if writing to a file.
        """
        super(CCD,self).__init__(winds)
        self.nxtot = nxtot
        self.nytot = nytot
        self.head = head

    def min(self):
        """
        Returns the minimum value of the :class:`CCD`.
        """
        winds = iter(self.values())
        vmin = next(winds).min()
        for wind in winds:
            vmin = min(vmin, wind.min())
        return vmin

    def max(self):
        """
        Returns the maximum value of the :class:`CCD`.
        """
        winds = iter(self.values())
        vmax = next(winds).max()
        for wind in winds:
            vmax = min(vmax, wind.max())
        return vmax

    def mean(self):
        """
        Returns the mean value of the :class:`Windata`.
        """
        npix, total = 0, 0.
        for wind in self.values():
            npix += wind.size
            total += wind.sum()
        return total / float(npix)

    def percentile(self, q):
        """
        Computes percentile(s) of the :class:`CCD`.

        Arguments::

        q : (float or sequence of floats)
          Percentile(s) to use, in range [0,100]
        """

        # Flatten into a single 1D array
        arrs = []
        for wind in self.values():
            arrs.append(wind.data.flatten())
        arr = np.concatenate(arrs)

        # Then compute percentiles
        return np.percentile(arr, q)

    def add_hdulist(self, hdul, start=None):
        """Add the CCD as a series of ImageHDUs to an hdulist.
        This is to allow a series of CCDs to be written. The CCD
        header is written into the first HDU added to the list.

        Arguments::

            hdul : (astropy.io.fits.HDUList)
                 An HDUList to which the CCD will be added as a series
                 of ImageHDUs, the first of which will contain the header.

            start : (astropy.io.fits.Header)
                 Standard header to add at the start of every the header of
                 each ImageHDU. e.g. This allows each one to be flagged with
                 a CCD number.

        Returns:: the modified HDUList.

        """
        first = True
        for key, wind in self.items():
            if start is None:
                fhead = fits.Header()
            else:
                fhead = start.copy()
            fhead['NWIN'] = (int(key), 'Window number')
            if first:
                fhead += self.head.copy()
                first = False
            hdul.append(wind.whdu(fhead))

        return hdul

    def wfits(self, fname, overwrite=False):
        """Writes out the CCD to a FITS file.

        Arguments::

            fname : (string)
                 Name of file to write to.

            overwrite : (bool)
                 True to overwrite pre-existing files
        """

        phdu = fits.PrimaryHDU(header=self.head)
        hdus = [phdu,]
        for key, wind in self.items():
            fhead = fits.Header()
            fhead['NWIN'] = (int(key), 'Window number')
            hdus.append(wind.whdu(fhead))
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(fname,overwrite=overwrite)

    def clash(self, ccd):
        """Returns "false" indicating two CCDs never clash
        """
        return false

    def __repr__(self):
        return 'CCD(winds=' + super(CCD,self).__repr__() + \
                            ', nxtot=' + repr(self.nxtot) + \
                            ', nytot=' + repr(self.nytot) + \
                            ', head=' + repr(self.head) + ')'


class MCCD(Group):
    """
    Class representing a multi-CCD as a Group of CCD objects
    plus a FITS header.
    """
    def __init__(self, ccds, head=None):
        """
        Constructs a :class:`MCCD`

        Arguments::

          ccds : (dict)
              dictionary of CCD objects. Keys should be integers.

          head : (astropy.io.fits.Header)
              a header which will be written as the primary header.
        """
        super(MCCD,self).__init__(ccds)
        for key, ccd in self.items():
            ccd.head['NCCD'] = (int(key), 'CCD number')
        self.head = head

    def wfits(self, fname, overwrite=False):
        """Writes out the MCCD to a FITS file.

        Arguments::

            fname : (string)
                 Name of file to write to.

            overwrite : (bool)
                 True to overwrite pre-existing files
        """

        phdu = fits.PrimaryHDU(header=self.head)
        hdus = [phdu,]
        for key, ccd in self.items():
            start = fits.Header()
            start['NCCD'] = int(key)
            ccd.add_hdulist(hdus, start)
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(fname,overwrite=overwrite)

    def __repr__(self):
        return 'MCCD(ccds=' + super(MCCD,self).__repr__() + \
                            ', head=' + repr(self.head) + ')'
