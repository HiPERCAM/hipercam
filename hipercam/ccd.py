# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Class to represent a CCD and a multi-CCD
"""

import warnings
import numpy as np

from astropy.io import fits
from .core import *
from .group import *
from .window import *

__all__ = ('CCD', 'MCCD', 'Hcam')

class CCD(Agroup):
    """
    Class representing a CCD as a :class:`Group` of :class:`Windat`
    objects plus a FITS header.
    """
    def __init__(self, winds, nxtot, nytot, head=None, copy=False):
        """Constructs a :class:`CCD`.

        Arguments::

          winds : (Group)
              Group of :class:`Windat` objects

          nxtot : (int)
              Unbinned X-dimension of CCD

          nytot : (int)
              Unbinned Y-dimension of CCD

          head : (astropy.io.fits.Header)
              a header which will be written along with the first of the
              Windat sub-images if writing to a file. If head=None on input,
              an empty header will be created.

          copy : (bool)
              if True, copy all the data over, otherwise only references are
              held. Holding references is fine if the calling program keeps
              re-generating the data but could cause problems in some
              circumstances.

        nxtot, nytot and head are stored as identically-named attributes of the CCD.

        """
        super().__init__(winds)
        self.nxtot = nxtot
        self.nytot = nytot
        if head is None:
            self.head = fits.Header()
        else:
            self.head = head

        if copy:
            self = self.copy()

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

        # Add comments if not already present.
        comm1 = 'Data representing a single CCD frame written by hipercam.CCD.wfits.'
        if 'COMMENT' not in phead or phead['COMMENT'].find(comm1) == -1:
            phead.add_comment(comm1)
            phead.add_comment('Multiple sub-windows are written in a series of HDUs. LLX and LLY give')
            phead.add_comment('the pixel location in unbinned pixels of their lower-left corners. XBIN')
            phead.add_comment('and YBIN are the binning factors. NXTOT and NYTOT are the total unbinned')
            phead.add_comment('dimensions of the CCD. The first HDU contains any extra header items')
            phead.add_comment('associated with the CCD.')

        hdul.writeto(fname,overwrite=overwrite)

    @classmethod
    def rfits(cls, fname):
        """Builds a :class:`CCD` from a FITS file. Expects a primary HDU,
        containing no data, followed by a series of HDUs each containing data
        for a series of non-overlapping windows. A header is extracted from
        the first HDU with data.
        """
        hdul = fits.open(fname)
        ccd = cls.rhdul(hdul)
        hdul.close()
        return ccd

    @classmethod
    def rhdul(cls, hdul, multi=False):
        """Builds a :class:`CCD` or several :class:`CCD`s from an :class:`HDUList`.
        Given an :class:`HDUList`, the header from the first HDU will be used
        to create the header for the :class:`CCD`. The data from all HDUs will
        be read into the :class:`Windat`s that make up the CCD. Each HDU will
        be searched for a header parameter WINDOW to label the
        :class:`Windat`s, but the routine will attempt to generate a
        sequential label if WINDOW is not found. If the auto-generated label
        conflicts with one already found, then a KeyError will be raised.

        The method can also run in a mode where the :class:`HDUList` is
        assumed to contain several :class:`CCD`s. In this case each
        :class:`CCD` comes in a series of continguous HDUs. Each HDU of a
        :class:`CCD` must be labelled with the keyword 'CCD'.

        Arguments::

          hdul : :class:`HDUList`
               each ImageHDU will be read as sub-window of the :class:`CCD`

          multi: (bool)
               if True, the routine will work as a generator, returning a
               :class:`CCD` each time a set of HDUs with identical header
               parameter 'CCD' has been found and processed.

        Returns::

          A :class:`CCD` if multi=False, or a series of (label, :class:`CCD`)
          2-element tuples if multi=True. In the latter case use as follows::

            for label, ccd in MCCD.rhdul(hdul, True):
               ... do something

        """
        winds = Group()
        nwin = 0
        first = True

        for hdu in hdul:
            head = hdu.header
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
            if 'WINDOW' in head:
                label = head['WINDOW']
            elif nwin in odict:
                raise KeyError('CCD.rhdul: window label conflict')
            else:
                label = nwin

            winds[label] = Windat.rhdu(hdu)

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
            yield (ccd_label, cls(winds, nxtot, nytot, first_head))
        else:
            return cls(winds, nxtot, nytot, first_head)

    def clash(self, ccd):
        """Dummy routine to allow :class:`CCD`s to be added into :class:`Group`
        objects.

        """
        pass

    def matches(self, ccd):
        """Check that the :class:`CCD` matches another, which in this means checking
        that each window of the same label matches the equivalent in the other `CCD`.

        There will be no match if there any windows present in one `CCD` but
        not the other and the windows must have identical location, size and
        binning.

        Raises KeyError or ValueError exception depending on the problem.
        """
        for key in self:
            if key not in ccd:
                raise KeyError(
                    'hipercam.CCD.matches: window {0:d} not found in ccd'.format(key))

        for key in ccd:
            if key not in self:
                raise KeyError(
                    'hipercam.CCD.matches: window {0:d} not found in self'.format(key))

        if self.nxtot != ccd.nxtot or self.nytot != ccd.nytot:
            raise ValueError('hipercam.CCD.matches: self / ccd have conflicting total (nxtot,nytot): ({0:d},{1:d}) vs ({2:d},{3:d})'.format(self.nxtot,self.nytot,ccd.nxtot,ccd.nytot))

        for key, wind in self.items():
            wind.matches(ccd[key])

    def copy(self, memo=None):
        """Make a copy of the :class:`CCD`

        copy.copy and copy.deepcopy of a `CCD` use this method
        """
        return CCD(
            super().copy(memo), self.nxtot, self.nytot, self.head.copy())

    def __repr__(self):
        return 'CCD(winds=' + super().__repr__() + \
                            ', nxtot=' + repr(self.nxtot) + \
                                       ', nytot=' + repr(self.nytot) + \
                                                  ', head=' + repr(self.head) + ')'

class MCCD(Agroup):
    """
    Class representing a multi-CCD as a Group of CCD objects
    plus a FITS header.
    """
    def __init__(self, ccds, head=None, copy=False):
        """
        Constructs a :class:`MCCD`

        Arguments::

          ccds : (Group)
              Group of CCD objects.

          head : (astropy.io.fits.Header)
              a header which will be written as the primary header. If head=None
              on input, and empty header will be created.

          copy : (bool)
              if True, copy all the data over, otherwise only references are
              held. Holding references is fine if the calling program keeps
              re-generating the data but could cause problems in some
              circumstances.
        """
        super().__init__(ccds)
        if head is None:
            self.head = fits.Header()
        else:
            self.head = head
        if copy:
            self = self.copy()

    def wfits(self, fname, overwrite=False):
        """Writes out the MCCD to a FITS file.

        Arguments::

            fname : (string)
                 Name of file to write to.

            overwrite : (bool)
                 True to overwrite pre-existing files
        """

        phead = self.head.copy()

        # Add comments if not already present.
        comm1 = 'Data representing multiple CCDs written by hipercam.MCCD.wfits.'
        if 'COMMENT' not in phead or str(phead['COMMENT']).find(comm1) == -1:
            phead.add_comment(comm1)
            phead.add_comment('Each window of each CCD is written in an HDU following the primary HDU')
            phead.add_comment('which contains a header only. The pixel location in unbinned pixels of')
            phead.add_comment('their lower-left corners is stored as LLX, LLY. XBIN and YBIN are the')
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
    def rfits(cls, fname):
        """Builds a :class:`CCD` from a FITS file. Expects a primary HDU,
        containing no data, followed by a series of HDUs each containing data
        for a series of non-overlapping windows. A header is extracted from
        the first HDU with data.
        """
        hdul = fits.open(fname)
        ccd = cls.rhdul(hdul)
        hdul.close()
        return ccd

    @classmethod
    def rhdul(cls, hdul):
        """Builds an :class:`MCCD` from an :class:`HDUList`. This will usually be
        passed the :class:`HDUList` from a file. The header from the first
        (primary) HDU will be used to create the header for the
        :class:`MCCD`. It is then assumed that the data for each :class:`CCD`
        is contained in the succeeding HDUs, i.e. that there is no data in the
        primary HDU. The data from all HDUs will be read into the
        :class:`Windat`s that make up the CCD. Each HDU will be searched for a
        header parameter WINDOW to label the :class:`Windat`s, but will
        attempt to generate a sequential label if WINDOW is not found.

        Arguments::

          hdul : :class:`HDUList`
               each ImageHDU will be read as sub-window of the :class:`CCD`

        """
        # Get the header from the first hdu, but otherwise ignore it
        head = hdul[0].header

        # Attempt to read rest of HDUs into a series of CCDs
        ccds = Group()
        for label, ccd in CCD.rhdul(hdul[1:],True):
            ccds[label] = ccd

        return cls(ccds, head)

    def copy(self, memo=None):
        """Make a copy of the :class:`MCCD`

        copy.copy and copy.deepcopy of an `MCCD` use this method
        """
        return MCCD(super().copy(memo), self.head.copy())

    def matches(self, mccd):
        """Check that the :class:`MCCD` matches another, which in this means checking
        that each CCD of the same label matches the equivalent in the other `MCCD`.

        There will be no match if there any CCDs present in one `MCCD` which are
        not in the other.

        Raises KeyError or ValueError exception depending on the problem.
        """
        for key in self:
            if key not in mccd:
                raise KeyError(
                    'hipercam.MCCD.matches: CCD {0:d} not found in mccd'.format(key))

        for key in mccd:
            if key not in self:
                raise KeyError(
                    'hipercam.MCCD.matches: CCD {0:d} not found in self'.format(key))

        for key, ccd in self.items():
            ccd.matches(mccd[key])

    def __repr__(self):
        return self.__class__.__name__ + \
            '(ccds=' + super().__repr__() + \
                     ', head=' + repr(self.head) + ')'


class Hcam(MCCD):
    """Specialisation of an MCCD to account for particular features of
    HiPERCAM data"""

    # Expected dimensions of all CCDs
    NXTOT = 1024
    NYTOT = 2048

    # Number of CCDs
    NCCD = 5

    def __init__(self, ccds, head=None, copy=False):
        """
        Constructs a :class:`Hcam`

        Arguments::

          ccds : (Group)
              Group of CCD objects.

          head : (astropy.io.fits.Header)
              a header which will be written as the primary header. If head=None
              on input, and empty header will be created.

          copy : (bool)
              if True, copy all the data over, otherwise only references are
              held. Holding references is fine if the calling program keeps
              re-generating the data but could cause problems in some
              circumstances.
        """
        super().__init__(ccds, head, copy)

        for nccd, ccd in self.items():
            if nccd < 1 or nccd > Hcam.NCCD:
                warnings.warn(
                    'hipercam.Hcam: nccd = {0:d} outside expected range (1-{1:d})'.format(nccd,Hcam.NCCD))

            if ccd.nxtot != Hcam.NXTOT or ccd.nytot != Hcam.NYTOT:
                warnings.warn(
                    'hipercam.Hcam: nxtot,nytot = {0:d},{1:d} do not match expect values ({2:d},{3:d})'.format(ccd.nxtot,ccd.nytot,Hcam.NXTOT,Hcam.NYTOT))


            for nwin, wind in ccd.items():
                if wind.llx < Hcam.NXTOT/2 and wind.urx > Hcam.NXTOT/2:
                    warnings.warn(
                        'hipercam.Hcam: window = {0:s} straddles output boundary in X'.format(wind.format()))
                if wind.lly < Hcam.NYTOT/2 and wind.ury > Hcam.NYTOT/2:
                    warnings.warn(
                        'hipercam.Hcam: window = {0:s} straddles output boundary in Y'.format(wind.format()))
