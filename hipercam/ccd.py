# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Class to represent a CCD and a multi-CCD
"""

import warnings
import numpy as np
from collections import OrderedDict

from astropy.io import fits
from .core import *
from .group import *
from .window import *

__all__ = ('CCD', 'MCCD', 'get_ccd_info')

class CCD(Agroup):
    """Class representing a single CCD as a :class:`Group` of :class:`Windat`
    objects plus a FITS header. It supports a few operations such as returning
    the mean value, arithematic operations, and file I/O to/from FITS files.

    FITS file structure: CCD objects are saved to FITS as a series of
    HDUs. The first HDU contains header information, but no data, while the
    Windats are stored in the suceeding HDUs. The first HDU's header should
    contain keywords NXTOT and NYTOT defining the total unbinned dimensions of
    the CCD and HIPERCAM, which should be set to 'CCD' to distinguish the file
    from the similar :class:`MCCD` objects. The headers of the following HDUs
    should each contain keywords WINDOW, LLX, LLY, XBIN and YBIN. WINDOW
    should be a unique integer label to attach to the Windat. LLX and LLY are
    the coordinates of the lower-left unbinned pixel in the Windat (starting
    at (1,1) for the bottom left of the entire chip). XBIN and YBIN are
    integer binning factors.

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

        if len(self):
            # swift sanity check, 'cos it's easy to get confused
            first = next(iter(self.values()))
            if not isinstance(first, Windat):
                raise ValueError(
                    'CCD.__init__: values of "winds" must all be Windat objects')

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

    def whdul(self, hdul, label=None):
        """Write the :class:`CCD` as a series of HDUs to an hdulist. The header is
        written in the first HDU, with no accompanying data and with keywords
        for the total dimensions (NXTOT, NYTOT) and the numbers of windows
        (NUMWIN) added in. If it is the first HDU added, it will be made a
        PrimaryHDU, otherwise it will be an ImageHDU of extension type
        CCDH. The :class:`Windat`s then follow as ImageHDUs of type
        'WIND'. This makes the HDUs associated with a give CCD stand out in
        e.g. 'fv'.

        Arguments::

            hdul : (astropy.io.fits.HDUList)
                 An HDUList to which the CCD will be added as a series of
                 ImageHDUs, the first of which will contain the header. If
                 empty at the start, the first HDU will be a PrimaryHDU.

            label : (int / string / None)
                 A label to attach to every HDU associated with the CCD. It
                 appears as a header item 'CCD'

        Returns:: the modified HDUList.

        """
        # First write the header
        head = self.head.copy()
        if label is not None:
            head['CCD'] = (label, 'CCD label')
        head['NXTOT'] = (self.nxtot, 'Total unbinned X dimension')
        head['NYTOT'] = (self.nytot, 'Total unbinned Y dimension')
        head['NUMWIN'] = (len(self), 'Total number of windows')
        if len(hdul):
            hdul.append(fits.ImageHDU(header=head,name='CCDH'))
        else:
            hdul.append(fits.PrimaryHDU(header=head))

        # Now the Windats
        for key, wind in self.items():
            whead = fits.Header()
            if label is not None:
                whead['CCD'] = (label, 'CCD label')
            whead['WINDOW'] = (key, 'Window label')
            hdul.append(wind.whdu(whead))

        return hdul

    @classmethod
    def rhdul(cls, hdul, multi=False):
        """
        Builds a :class:`CCD` or several :class:`CCD`s from an
        :class:`HDUList`.  Given an :class:`HDUList` representing a single CCD
        the header from the first HDU will be used to create the header for
        the :class:`CCD`. The data from the remaining HDUs will be read into
        the :class:`Windat`s that make up the CCD. Each data HDU will be
        searched for a header parameter WINDOW to label the :class:`Windat`s,
        but the routine will attempt to generate a sequential label if WINDOW
        is not found. If the auto-generated label conflicts with one already
        found, then a KeyError will be raised.

        The method can also run in a mode where the :class:`HDUList` is
        assumed to contain several :class:`CCD`s. In this case each
        :class:`CCD` comes in a series of continguous HDUs, starting with a
        header-only one followed by the data windows, and all HDUs of a given
        :class:`CCD` must be labelled with the keyword 'CCD' to allow the CCD
        to be defined.

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
        nwin = 1
        first = True

        for hdu in hdul:
            # The header of the HDU
            head = hdu.header

            if multi and not first:
                # If we expect multiple CCDs, then each one is defined by
                # matching values of the header keyword CCD. As soon as a
                # mismatch is found, we assume that we have started on a new
                # CCD and return the old one.
                if ccd_label != head['CCD']:
                    ccd = cls(winds, nxtot, nytot, main_head)
                    winds = Group()
                    first = True
                    nwin = 1
                    yield (ccd_label,ccd)

            if first:
                # Extract header from first HDU (any data it might contain are
                # ignored)
                nxtot = head['NXTOT']
                nytot = head['NYTOT']
                if multi: ccd_label = head['CCD']

                del head['NXTOT']
                del head['NYTOT']
                if 'NUMWIN' in head: del head['NUMWIN']
                main_head = head
                first = False

            else:

                # Except for the first HDU of a CCD, all should contain data,
                # and for a single CCD, we assume all are part of the CCD.
                if 'WINDOW' in head:
                    label = head['WINDOW']
                elif str(nwin) in winds:
                    raise KeyError('CCD.rhdul: window label conflict')
                else:
                    label = str(nwin)

                winds[label] = Windat.rhdu(hdu)

            # step the window counter
            nwin += 1

        # Finishing up having gone through all CCDs.
        if multi:
            yield (ccd_label, cls(winds, nxtot, nytot, main_head))
        else:
            return cls(winds, nxtot, nytot, first_head)


    def wfits(self, fname, label=None, overwrite=False):
        """Writes out the CCD to a FITS file.

        Arguments::

            fname : (string)
                 Name of file to write to.

            label : (int | None)
                 Label for the CCD to write into the headers

            overwrite : (bool)
                 True to overwrite pre-existing files
        """

        # compile the HDUList
        hdul = fits.HDUList()
        hdul = self.whdul(hdul, label)

        # Get the main header (header of first HDU)
        phead = hdul[0].header
        phead['HIPERCAM'] = ('CCD', 'Type of HiPERCAM data (CCD | MCCD)')

        # Add comments if not already present.
        comm1 = 'Data representing a single CCD frame written by hipercam.CCD.wfits.'
        if 'COMMENT' not in phead or phead['COMMENT'].find(comm1) == -1:
            phead.add_comment(comm1)
            phead.add_comment('Multiple sub-windows are written as a series of HDUs. LLX and LLY give')
            phead.add_comment('the pixel location in unbinned pixels of their lower-left corners. XBIN')
            phead.add_comment('and YBIN are the binning factors. The first HDU contains the main header')
            phead.add_comment('along with NXTOT and NYTOT, the total unbinned and NUMWIN the number of')
            phead.add_comment('windows which should match the number of HDUs following the first.')

        hdul.writeto(fname,overwrite=overwrite)

    @classmethod
    def rfits(cls, fname):
        """Builds a :class:`CCD` from a FITS file. Expects a primary HDU,
        containing no data, followed by a series of HDUs each containing data
        for a series of non-overlapping windows.
        """
        with fits.open(fname) as hdul:
            return cls.rhdul(hdul)

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
            super().copy(memo), self.nxtot, self.nytot, self.head.copy()
        )

    def float32(self):
        """Applies :class:Windat.float32 to all Windats of a CCD"""
        for wind in self.values():
            wind.float32()

    def uint16(self):
        """Applies :class:Windat.uint16 to all Windats of a CCD"""
        for wind in self.values():
            wind.uint16()

    def __repr__(self):
        return '{:s}(winds={:s}, nxtot={!r}, nytot={!r}, head={!r})'.format(self.__class__.__name__,super().__repr__(),self.nxtot,self.nytot,self.head)

class MCCD(Agroup):
    """Class representing a multi-CCD as a Group of CCD objects plus a FITS
    header. It supports arithematic operations and file I/O to/from FITS
    files.

    FITS file structure: :class:`MCCD`s are saved to FITS as a series of
    HDUs. The first HDU contains header information, but no data. The
    succeeding HDUs comes in block representing each CCD. See the docs on
    :class:`CCD` to understand their layout and important
    keywords. :class:`MCCD` sets the keyword HIPERCAM in the primary HDU's
    header to 'MCCD' to indicate the nature of the file. It also adds CCD to
    the header of the first HDU (containing no data) of every CCD it
    contains. CCD should be set to a unique integer label for each CCD.

    Here is a schematic summary of the structure of a 3 CCD MCCD file (which
    typically will contain other header data in addition):

      Primary HDU [HDU 1] -- no data
         HIPERCAM = 'MCCD'
         NUMCCD = 3 [for convenience, number of CCDs]

      HDU 2, extension CCDH (for CCD header) [no data]
         CCD = 2 [any integer is OK]
         NUMWIN = 2 [for convenience, number of windows]
         NXTOT = 1024 [max X dimension]
         NYTOT = 2048 [max Y dimension]

      HDU 3, extension WIND for Windat
         CCD = 2 [the CCD it belongs to]
         WINDOW = 3 [any integer is OK]
         LLX = 11 [X of left-hand unbinned coords]
         LLY = 21 [Y of lowest unbinned coords]
         XBIN = 2 [X binning factor]
         YBIN = 3 [Y binning factor]
         + 2D data array contain the data.

      HDU 4, extension WIND for Windat
         CCD = 2 [the CCD it belongs to]
         WINDOW = 4 [any integer is OK]
         LLX = 51 [X of left-hand unbinned coords]
         LLY = 301 [Y of lowest unbinned coords]
         XBIN = 2 [X binning factor]
         YBIN = 1 [Y binning factor]
         + 2D data array contain the data.

      HDU 4, extension CCDH (for CCD header) [no data]
         CCD = 5 [any integer is OK]
         NUMWIN = 2 [for convenience, number of windows]
         NXTOT = 824 [max X dimension]
         NYTOT = 2048 [max Y dimension]

      then two WIND HDUs etc. Note that the MCCD format is flexible
      when it comes to different binning factors etc, although a given
      instrument may well have its own restrictions on this.
    """

    def __init__(self, ccds, head=None, copy=False):
        """
        Constructs a :class:`MCCD`

        Arguments::

          ccds : (Group)
              Group of CCD objects, or something that can be loaded into
              a Group.

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

        if len(self):
            # swift sanity check, 'cos it's easy to get confused
            first = next(iter(self.values()))
            if not isinstance(first, CCD):
                raise ValueError(
                    'MCCD.__init__: values of "ccds" must all be CCD objects')

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
        phead['NUMCCD'] = (len(self), 'Number of CCDs')
        phead['HIPERCAM'] = ('MCCD', 'Type of HiPERCAM data (CCD | MCCD)')

        # Add comments if not already present.
        comm1 = 'Data representing multiple CCDs written by hipercam.MCCD.wfits.'
        if 'COMMENT' not in phead or str(phead['COMMENT']).find(comm1) == -1:
            phead.add_comment(comm1)
            phead.add_comment('Each window of each CCD is written in a series of HDUs following an') 
            phead.add_comment('HDU containing only the header. These follow an overall header for the')
            phead.add_comment('MCCD containing top-level information. The headers of the data window')
            phead.add_comment('HDUs have keywords LLX, LLY giving the pixel location in unbinned')
            phead.add_comment('pixels and XBIN and YBIN for their binning factors. The total unbinned')
            phead.add_comment('dimensions of each CCD are stored under keywords NXTOT and NYTOT for')
            phead.add_comment('each CCD. Each HDU associated with a given CCD is labelled with a')
            phead.add_comment('header keyword CCD.')

        # make the first HDU
        phdu = fits.PrimaryHDU(header=phead)
        hdul = [phdu,]

        # add in the HDUs of all the CCDs
        for key, ccd in self.items():
            hdul = ccd.whdul(hdul, key)
        hdulist = fits.HDUList(hdul)
        hdulist.writeto(fname,overwrite=overwrite)

    @classmethod
    def rfits(cls, fname):
        """Builds an :class:`MCCD` from a FITS file. Expects a primary HDU,
        containing no data, followed by a series of HDU blocks consisting
        of the header for each CCD followed by data HDUs.
        """
        with fits.open(fname) as hdul:
            return cls.rhdul(hdul)

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
        # Get the main header from the first hdu, but otherwise ignore it
        head = hdul[0].header
        if 'NUMCCD' in head: del head['NUMCCD']

        # Attempt to read the rest of HDUs into a series of CCDs
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
        return '{:s}(ccds={:s}, head={!r})'.format(
            self.__class__.__name__, super().__repr__(), self.head)

    def float32(self):
        """Applies :class:Windat.float32 to all Windats of an MCCD"""
        for ccd in self.values():
            ccd.float32()

    def uint16(self):
        """Applies :class:Windat.uint16 to all Windats of an MCCD"""
        for ccd in self.values():
            ccd.uint16()

def get_ccd_info(fname):
    """Routine to return some useful basic information from an MCCD file
    without reading the whole thing in. It returns an OrderedDict keyed on the
    CCD label which returns the maximum X and Y dimensions of the CCD as a
    2-element tuple. If no CCD label [header parameter 'CCD'] is found, the
    rourins assigns a label of '1' on the assumption that the file cotains a
    single CCD.
    """
    info = OrderedDict()

    with fits.open(fname) as hdul:
        first = True
        for hdu in hdul:
            # The header of the HDU
            head = hdu.header

            if not first and (
                ('CCD' in head and ccd_label != head['CCD']) or
                ('NXTOT' in head and 'NYTOT' in head)):
                # Reset because we think we are on the first HDU of a CCD
                first = True

            if first and 'NXTOT' in head and 'NYTOT' in head:
                # Extract header from first HDU of a CCD
                if 'CCD' in head:
                    ccd_label = head['CCD']
                else:
                    ccd_label = '1'
                nxtot = head['NXTOT']
                nytot = head['NYTOT']
                info[ccd_label] = (nxtot,nytot)
                first = False

    return info

