# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Class to represent a CCD and a multi-CCD
"""

# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import numpy as np
from collections import OrderedDict

from astropy.io import fits
from .core import *
from .group import *
from .window import *

class CCD(Group):
    """
    Class representing a CCD as a :class:`Group` of :class:`Windat`
    objects plus a FITS header.
    """
    def __init__(self, winds, nxtot, nytot, head=None):
        """
        Constructs a :class:`CCD`.

        Arguments::

          winds : (OrderedDict)
              :class:`Windat` objects keyed by integers and/or strings.

          nxtot : (int)
              Unbinned X-dimension of CCD

          nytot : (int)
              Unbinned Y-dimension of CCD

          head : (astropy.io.fits.Header)
              a header which will be written along with the first of the Windat
              sub-images if writing to a file.

        The latter three are stored as identically-named attributes of the CCD.
        """
        super(CCD,self).__init__(winds)
        self.nxtot = nxtot
        self.nytot = nytot
        self.head = head

    def min(self):
        """
        Returns the minimum value of the :class:`CCD`, i.e. the
        minimum of all the minimum values of all :class:`Windat`s
        """
        return min(v.min() for v in self.values())

    def max(self):

        """
        Returns the maximum value of the :class:`CCD`.
        """
        return max(v.max() for v in self.values())

    def mean(self):
        """
        Returns the mean value of the :class:`Windat`.
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
        arrs = [wind.data.flatten() for wind in self.values()]
        arr = np.concatenate(arrs)

        # Then compute percentiles
        return np.percentile(arr, q)

    def add_hdulist(self, hdul, label=None):
        """Add the CCD as a series of ImageHDUs to an hdulist.  The CCD header is
        written into the first HDU added to the list.

        Arguments::

            hdul : (astropy.io.fits.HDUList)
                 An HDUList to which the CCD will be added as a series of ImageHDUs,
                 the first of which will contain the header.

            label : (int / string / None)
                 A label to attach to every HDU associated with the CCD. It appears
                 as a header item 'CCD'

        Returns:: the modified HDUList.

        """
        first = True
        for key, wind in self.items():
            if first:
                fhead = self.head.copy()
                fhead['NXTOT'] = (self.nxtot, 'Total unbinned X dimension')
                fhead['NYTOT'] = (self.nytot, 'Total unbinned Y dimension')
                first = False
            else:
                fhead = fits.Header()

            if label is not None:
                fhead['CCD'] = (label, 'CCD label')
            fhead['WINDOW'] = (key, 'Window label')

            hdul.append(wind.whdu(fhead))

        return hdul

    def wfits(self, fname, label=None, overwrite=False):
        """Writes out the CCD to a FITS file.

        Arguments::

            fname : (string)
                 Name of file to write to.

            label : (int / string / None)
                 Label for the CCD to write into the headers

            overwrite : (bool)
                 True to overwrite pre-existing files
        """

        hdul = fits.HDUList()
        hdul = self.add_hdulist(hdul, label)
        phead = hdul[0].header
        phead.add_blank('-----------------------------')
        phead.add_comment('Data representing a single CCD frame written by hipercam.CCD.wfits.')
        phead.add_comment('Multiple sub-windows are written in a series of HDUs. LLX and LLY give')
        phead.add_comment('the pixel location in unbinned pixels of their lower-left corners. XBIN')
        phead.add_comment('and YBIN are the binning factors. NXTOT and NYTOT are the total unbinned')
        phead.add_comment('dimensions of the CCD. The first HDU contains any extra header items')
        phead.add_comment('associated with the CCD.')

        hdul.writeto(fname,overwrite=overwrite)

    @classmethod
    def from_fits(cls, fname):
        """Builds a :class:`CCD` from a FITS file. Expects a primary HDU,
        containing no data, followed by a series of HDUs each containing data
        for a series of non-overlapping windows. A header is extracted from
        the first HDU with data.
        """
        hdul = fits.open(fname)
        ccd = cls.from_hdul(hdul)
        hdul.close()
        return ccd

    @classmethod
    def from_hdul(cls, hdul, multi=False):
        """Builds a :class:`CCD` or several :class:`CCD`s from an :class:`HDUList`.
        This will usually be passed the :class:`HDUList` from a file. The
        header from the first HDU will be used to create the header for the
        :class:`CCD`. The data from all HDUs will be read into the
        :class:`Windat`s that make up the CCD. Each HDU will be searched for a
        header parameter WINDOW to label the :class:`Windat`s, but will
        attempt to generate a sequential label if WINDOW is not found. If the
        auto-generated label conflicts with one already found, then a
        HipercamError will be raised.

        Arguments::

          hdul : :class:`HDUList`
               each ImageHDU will be read as sub-window of the :class:`CCD`

          multi: (bool)
               if True, the routine will work as a generator, returning a :class:`CCD`
               each time a set of HDUs with identical header parameters 'CCD' have been
               found.

        Returns::

          A :class:`CCD` if multi=False, or a series of (label, :class:`CCD`) 2-element 
          tuples if multi=True

        """
        winds = Group()
        nwin = 0
        first = True

        for hdu in hdul:
            if multi and not first:
                # Check that CCD header item matches, if not
                # then we return the stuff bagged to date and reset
                # to be ready for a new CCD
                if ccd_label != head['CCD']:
                    ccd = cls(winds, nxtot, nytot, first_head)
                    winds = Group()
                    first = True
                    nwin = 0
                    yield (ccd_label,ccd)

            nwin += 1
            head = hdu.header
            if 'WINDOW' in head:
                label = head['WINDOW']
            elif nwin in odict:
                raise HipercamError('CCD.from_hdul: window label conflict')
            else:
                label = nwin

            winds[label] = Windat.from_hdu(hdu)

            if first:
                # Extract header from first HDU with data
                # First remove standard items
                nxtot = head['NXTOT']
                nytot = head['NYTOT']
                if multi:
                    ccd_label = head['CCD']

                del head['NXTOT']
                del head['NYTOT']
                del head['LLX']
                del head['LLY']
                del head['XBIN']
                del head['YBIN']
                if 'WINDOW' in head: del head['WINDOW']
                if 'CCD' in head: del head['CCD']
                first_head = head
                first = False

        if multi:
            return (ccd_label, cls(winds, nxtot, nytot, first_head))
        else:
            return cls(winds, nxtot, nytot, first_head)

    def clash(self, ccd):
        """Simply returns "False" indicating two :class:`CCD`s never clash. Needed
        by the container class :class:`Group`.
        """
        return False

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
              dictionary of CCD objects. Keys should be integers or strings

          head : (astropy.io.fits.Header)
              a header which will be written as the primary header.
        """
        super(MCCD,self).__init__(ccds)
        self.head = head

    def wfits(self, fname, overwrite=False):
        """Writes out the MCCD to a FITS file.

        Arguments::

            fname : (string)
                 Name of file to write to.

            overwrite : (bool)
                 True to overwrite pre-existing files
        """

        phead = self.head.copy()
        phead.add_blank('-----------------------------')
        phead.add_comment('Data representing multiple CCDs written by hipercam.MCCD.wfits.')
        phead.add_comment('Each sub-windows of each CCD is written in an HDU following the primary')
        phead.add_comment('HDU which contains a header only. The pixel location in unbinned pixels')
        phead.add_comment('of their lower-left corners is stored as LLX, LLY. XBIN and YBIN are the')
        phead.add_comment('binning factors. NXTOT and NYTOT are the total unbinned dimensions of')
        phead.add_comment('the CCD. The first HDU of each CCD contains any extra header items')
        phead.add_comment('associated with it. Each HDU associated with a CCD is labelled with a')
        phead.add_comment('header parameter CCD.')

        phdu = fits.PrimaryHDU(header=phead)
        hdul = [phdu,]
        for key, ccd in self.items():
            hdul = ccd.add_hdulist(hdul, key)
        hdulist = fits.HDUList(hdul)
        hdulist.writeto(fname,overwrite=overwrite)

    @classmethod
    def from_fits(cls, fname):
        """Builds a :class:`CCD` from a FITS file. Expects a primary HDU,
        containing no data, followed by a series of HDUs each containing data
        for a series of non-overlapping windows. A header is extracted from
        the first HDU with data.
        """
        hdul = fits.open(fname)
        ccd = cls.from_hdul(hdul)
        hdul.close()
        return ccd

    @classmethod
    def from_hdul(cls, hdul):
        """Builds an :class:`MCCD` from an :class:`HDUList`. This will usually be passed
        the :class:`HDUList` from a file. The header from the first (primary) HDU will
        be used to create the header for the :class:`MCCD`. It is then assumed that the
        data for each :class:`CCD` is contained in the succeeding HDUs, i.e. that there is
        no data in the primary HDU. The data from all
        HDUs will be read into the :class:`Windat`s that make up the CCD. Each
        HDU will be searched for a header parameter WINDOW to label the
        :class:`Windat`s, but will attempt to generate a sequential label if
        WINDOW is not found. If the auto-generated label conflicts with one
        already found, then a HipercamError will be raised.

        Arguments::

          hdul : :class:`HDUList`
               each ImageHDU will be read as sub-window of the :class:`CCD`

        """
        # Get the header from the first hdu, but otherwise ignore it
        head = hdul[0].header

        # Attempt to read rest of HDUs into a series of CCDs
        ccds = Group()
        for label, ccd in CCD.from_hdu(hdul[1:]):
            ccd[label] = ccd

        return cls(ccds, head)

    def __repr__(self):
        return 'MCCD(ccds=' + super(MCCD,self).__repr__() + \
                            ', head=' + repr(self.head) + ')'
