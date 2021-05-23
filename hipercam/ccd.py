# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Classes & methods to represent CCDs and multi-CCDs

A :class:`CCD` is built as a :class:`Group` of :class:`Window`
objects, while a multi-CCD of :class:`MCCD` is built as a
:class:`Group` of :class:`CCD` objects.

"""

import warnings
import numpy as np
from collections import OrderedDict

from astropy.io import fits
from .core import *
from .group import *
from .window import *
from .header import *

__all__ = ("CCD", "MCCD", "get_ccd_info", "trim_ultracam")

# Special keywords that will be stripped from FITS headers
# on input as they are added on output. This is to make
# to make diagnostic print statements easier to follow.
KEYWORDS = (
    "NXTOT",
    "NYTOT",
    "NUMWIN",
    "CCD",
    "WINDOW",
)


class CCD(Agroup):
    """Represents a single CCD as a :class:`Group` of :class:`Window`.

    :class:`CCD` objects support a few operations such as returning
    the mean value, arithematic operations, and file I/O to/from FITS
    files. Header information is passed via the Windows. The first of
    these is assumed to be the location for any data pertaining to the
    CCD as a whole.

    """

    def __init__(self, winds, nxtot, nytot, nxpad=0, nypad=0, copy=False):
        """Constructs a :class:`CCD`.

        Parameters:

          winds : :class:`Group`
              Group of :class:`Window` objects

          nxtot : int
              Unbinned X-dimension of imaging area of CCD

          nytot : int
              Unbinned Y-dimension of imaging area of CCD

          nxpad : int
              extra padding in the X-direction to define plot limits (accounts
              for possible prescans / overscans that have to be displayed
              outside the imaging area)

          nypad : int
              extra padding in the Y-direction to define plot limits (accounts
              for possible prescans / overscans that have to be displayed
              outside the imaging area)

          copy : bool
              if True, copy all the data over, otherwise only references are
              held. Holding references is fine if the calling program keeps
              re-generating the data but could cause problems in some
              circumstances.


        nxtot, nytot, nxpad, nypad are stored as identically-named attributes of
        the CCD.

        """
        super().__init__(Window, winds)

        self.nxtot = nxtot
        self.nytot = nytot
        self.nxpad = nxpad
        self.nypad = nypad

        if copy:
            self = self.copy()

    def __reduce__(self):
        """This is to overcome a problem with pickling CCDs. The arguments for init
        don't get picked up for some reason so this must return them in the
        second element of the tuple

        """
        return (
            self.__class__,
            (list(self.items()), self.nxtot, self.nytot, self.nxpad, self.nypad),
        )

    def flatten(self):
        """Returns all data of a CCD as a single 1D array, analogous to
        numpy.flatten
        """
        arrs = []
        for wind in self.values():
            arrs.append(wind.flatten())
        return np.concatenate(arrs)

    def min(self):
        """Returns the minimum value of the :class:`CCD`, i.e. the minimum of
        all the minimum values of all the :class:`Window`

        """
        return min(v.min() for v in self.values())

    def max(self):
        """
        Returns the maximum value of the :class:`CCD`.
        """
        return max(v.max() for v in self.values())

    def mean(self):
        """
        Returns the mean value of the :class:`CCD`.
        """
        npix, total = 0, 0.0
        for wind in self.values():
            npix += wind.size
            total += wind.sum()
        return total / float(npix)

    def median(self):
        """
        Returns the median value of the :class:`CCD`.
        """

        # Flatten into a single 1D array
        arrs = [wind.data.flatten() for wind in self.values()]
        arr = np.concatenate(arrs)

        # Then compute median
        return np.median(arr)

    def percentile(self, q, xlo=None, xhi=None, ylo=None, yhi=None):
        """
        Computes percentile(s) of the :class:`CCD`.

        Arguments::

          q : float or sequence of floats
            Percentile(s) to use, in range [0,100]

          xlo : int | None
            To restrict range of pixels

          xhi : int | None
            To restrict range of pixels

          ylo : int | None
            To restrict range of pixels

          yhi : int | None
            To restrict range of pixels
        """

        # Flatten into a single 1D array, accounting
        # for the sub-region defined by xlo etc
        arrs = []
        for wind in self.values():
            try:
                wind_small = wind.window(xlo, xhi, ylo, yhi)
                arrs.append(wind_small.data.flatten())
            except HipercamError as err:
                # could get errors if window not aligned
                # with xlo/xhi/ylo/yhi
                pass

        if len(arrs) == 0:
            raise ValueError("supplied xlo,xhi,ylo,yhi region contains no pixels")

        # concatenate into 1D array
        arr = np.concatenate(arrs)

        # Then compute percentiles
        return np.percentile(arr, q)

    def whdul(self, hdul=None, cnam=None, xoff=0, yoff=0):
        """Write the :class:`CCD` as a series of HDUs, one per
        :class:`Window`, adding to and returning an :class:`HDUList`.

        Arguments::

            hdul : astropy.io.fits.HDUList
               The HDUList to add to.

            cnam : string | None
               CCD name.

            xoff : int [if cname is not None]
               X-offset for mosaicing in ds9. Ignored if cnam is None.

        Returns:: an :class:`HDUList`.

        """
        if hdul is None:
            hdul = fits.HDUList()

        # Now 1 HDU per Window
        for n, (wnam, wind) in enumerate(self.items()):

            if cnam is None:
                extnam = "W:{:s}".format(wnam)
            else:
                extnam = "C:{:s}, W:{:s}".format(cnam, wnam)

            # Generate a few extra items to add to the header
            head = Header()
            if cnam:
                head["CCD"] = (cnam, "CCD label")

            if n == 0:
                # Tack a few general items into the first Window
                # to be written for a given CCD.
                head["NXTOT"] = (self.nxtot, "Total unbinned X dimension")
                head["NYTOT"] = (self.nytot, "Total unbinned Y dimension")
                head["NUMWIN"] = (len(self), "Total number of windows")
                head["NXPAD"] = (self.nxpad, "X-padding for display purposes")
                head["NYPAD"] = (self.nypad, "Y-padding for display purposes")

            head["WINDOW"] = (wnam, "Window label")
            hdul.append(wind.whdu(head, xoff, yoff, extnam))

        return hdul

    @classmethod
    def rhdul(cls, hdul, cnam=None):
        """Builds a single :class:`CCD` from an :class:`HDUList`.

        If cnam is None, it is assumed that each HDU represents a Window.
        Each data HDU will be searched for a header parameter 'WINDOW' to label
        the :class:`Window`, but the routine will attempt to generate a
        sequential label if WINDOW is not found. If the auto-generated label
        conflicts with one already found, then a KeyError will be raised.

        If cnam is not None, cnam is taken to be the label of a CCD to be
        extracted from the :class:`HDUList` which may contain multiple
        CCDs. It is assumed in this case that the CCD of interest has a
        continguous set of HDUs all containing the keyword CCD set equal to
        cnam. This allows individual CCDs to be read from hcm MCCD files.

        Arguments::

          hdul  : :class:`HDUList`
               each HDU will be read as sub-window of the :class:`CCD`

          cnam  : (None | string)
               if cnam is not None, the routine will try to find and return a
               :class:CCD built from a contiguous subset of HDUs that matches
               cnam (the keyword CCD is used)

        Returns a :class:`CCD`. Data are converted to float32 unless they are
        read in as float64.

        """

        # Group of Windows
        winds = Group(Window)
        nwin = 1
        first = True

        for hdu in hdul:

            # The header of the HDU
            head = hdu.header

            if cnam is None or ("CCD" in head and cnam == head["CCD"]):
                if first:
                    # Extract maximum dimensions and padding from first
                    # HDU of a CCD
                    nxtot = head["NXTOT"]
                    nytot = head["NYTOT"]
                    nxpad = head.get("NXPAD", 0)
                    nypad = head.get("NYPAD", 0)
                    first = False

                # attempt auto-generation of window labels as an aid
                # to ingesting foreign-format data
                if "WINDOW" in head:
                    label = head["WINDOW"]
                elif str(nwin) in winds:
                    raise KeyError("window label conflict")
                else:
                    label = str(nwin)

                # create the Window
                wind = Window.rhdu(hdu)

                # clean its header
                for key in KEYWORDS:
                    if key in wind:
                        del wind[key]

                # store it
                winds[label] = wind

                # step the window counter
                nwin += 1

            elif not first:
                # have encountered an HDU with a mis-matching label and
                # since we assume each CCD is a continguous block, we stop
                break

        return cls(winds, nxtot, nytot, nxpad, nypad)

    @classmethod
    def rmhdul(cls, hdul):
        """Builds a :class:`CCD` from an :class:`HDUList` potentially containing
        multiple CCDs, i.e. an MCCD object. This function acts as a generator
        to allow the CCDs to be returned in successive calls. It is assumed
        that hdul is positioned at the start of the CCDs (i.e. with the
        initial data-less primary HDU removed).

        Arguments::

          hdul  : :class:`HDUList`
               each ImageHDU will be read as a sub-window of the :class:`CCD`

        Returns::

          A series of (label, :class:`CCD`) 2-element tuples to be used as
          as follows::

            >> for label, ccd in CCD.rmhdul(hdul, True):
            >>    ... do something

        Data are converted to float32 unless they are read in as float64.

        """

        # Group on Windows
        winds = Group(Window)
        nwin = 1
        first = True

        for hdu in hdul:

            # The header of the HDU
            head = hdu.header

            if not first and ccd_label != head["CCD"]:
                # we have found a new CCD. Create and return the old CCD
                ccd = cls(winds, nxtot, nytot, nxpad, nypad)
                yield (ccd_label, ccd)

                # re-initialise for the next CCD
                winds = Group(Window)
                first = True
                nwin = 1

            if first:
                # we are on the first HDU of a CCD, get some essentials
                nxtot = head["NXTOT"]
                nytot = head["NYTOT"]
                nxpad = head.get("NXPAD", 0)
                nypad = head.get("NYPAD", 0)
                ccd_label = head["CCD"]
                first = False

            # store the HDU as a Window
            if "WINDOW" in head:
                label = head["WINDOW"]
            elif str(nwin) in winds:
                raise KeyError("window label conflict")
            else:
                label = str(nwin)

            # create the Window
            wind = Window.rhdu(hdu)

            # clean its header
            for key in KEYWORDS:
                if key in wind:
                    del wind[key]

            # store it
            winds[label] = wind

            # step the window counter
            nwin += 1

        # when we get here, we should have data still to be returned
        yield (ccd_label, cls(winds, nxtot, nytot, nxpad, nypad))

    def write(self, fname, overwrite=False):
        """Writes out the CCD to a FITS file.

        Arguments::

            fname : (string)
                 Name of file to write to.

            overwrite : (bool)
                 True to overwrite pre-existing files
        """

        # compile the HDUList
        hdul = self.whdul()

        # Get the main header (header of first HDU)
        phead = hdul[0].header
        phead["HIPERCAM"] = ("CCD", "Type of HiPERCAM data (CCD | MCCD)")

        # Add comments if not already present.
        comm1 = "Data representing a single CCD frame written by hipercam.CCD.write."
        if "COMMENT" not in phead or phead["COMMENT"].find(comm1) == -1:
            phead.add_comment(comm1)
            phead.add_comment(
                "Multiple sub-windows are written as a series of HDUs. LLX and LLY give"
            )
            phead.add_comment(
                "the pixel location in unbinned pixels of their lower-left corners. XBIN"
            )
            phead.add_comment(
                "and YBIN are the binning factors. The first HDU contains the main header"
            )
            phead.add_comment(
                "along with NXTOT and NYTOT, the total unbinned dimensions of the imaging"
            )
            phead.add_comment(
                "area of the CCD, NXPAD and NYPAD to account for pre- and over-scans for"
            )
            phead.add_comment(
                "display purposes, and NUMWIN the number of windows which should match"
            )
            phead.add_comment("the number of HDUs following the first.")

        hdul.writeto(fname, overwrite=overwrite)

    @classmethod
    def read(cls, fname, cnam=None):
        """Builds a :class:`CCD` from a FITS file. Expects a primary HDU,
        containing no data, followed by a series of HDUs each containing data
        for a series of non-overlapping windows.

        This can also be applied to extract a particular CCD, labelled 'cnam'
        from an MCCD hcm file, on the assumption that all HDUs for a given CCD
        come in a contiguous block with each one labelled with the keyword
        'CCD'.

        Data are converted to float32 unless they are read in as float64.
        """

        with fits.open(fname) as hdul:
            return cls.rhdul(hdul, cnam)

    def matches(self, ccd):
        """Check that the :class:`CCD` matches another, which in this case means
        checking that each window of the same label matches the equivalent in
        the other `CCD`.

        There will be no match if there any windows present in one `CCD` but
        not the other and the windows must have identical location, size and
        binning.

        Raises KeyError or ValueError exception depending on the problem.

        """
        for key in self:
            if key not in ccd:
                raise KeyError("window {0:d} not found in ccd".format(key))

        for key in ccd:
            if key not in self:
                raise KeyError("window {0:d} not found in self".format(key))

        if self.nxtot != ccd.nxtot or self.nytot != ccd.nytot:
            raise ValueError(
                "self / ccd have conflicting total (nxtot,nytot):"
                " ({0:d},{1:d}) vs ({2:d},{3:d})".format(
                    self.nxtot, self.nytot, ccd.nxtot, ccd.nytot
                )
            )

        for key, wind in self.items():
            wind.matches(ccd[key])

    def copy(self, memo=None):
        """Make a copy of the :class:`CCD`

        copy.copy and copy.deepcopy of a `CCD` use this method
        """
        return CCD(super().copy(memo), self.nxtot, self.nytot, self.nxpad, self.nypad)

    def float32(self):
        """Applies :class:Window.float32 to all Windows of a CCD"""
        for wind in self.values():
            wind.float32()

    def float64(self):
        """Applies :class:Window.float64 to all Windows of a CCD"""
        for wind in self.values():
            wind.float64()

    def uint16(self):
        """Applies :class:Window.uint16 to all Windows of a CCD"""
        for wind in self.values():
            wind.uint16()

    def set_const(self, val):
        """Sets all Windows to a constant value"""
        for wind in self.values():
            wind.set_const(val)

    def inside(self, x, y, dmin):
        """Tests whether a point x,y lies within any Window of a CCD, at
        least dmin from its outer edge. It returns with the Window label
        or None if the point is not inside any Window"""

        for wnam, wind in self.items():
            if wind.distance(x, y) > dmin:
                return wnam
        else:
            return None

    def crop(self, ccd):
        """Given a template :class:`CCD` called `ccd`, this tries to modify
        the format of the :class:`CCD`, returning a new
        :class:`CCD`. This is often needed e.g. prior to applying a
        full frame flat field. In general this operation requires
        re-labelling of the :class:`Window` composing the CCD. It
        raises a :class:`HipercamError` if it fails.

        """
        # build an empty CCD to get the headers. this will be added to
        tccd = CCD(Group(Window), self.nxtot, self.nytot)

        # wind through the Windows of the template CCD
        for twnam, twind in ccd.items():

            # for each one, search for a surrounding Window in
            # the CCD we are chopping down to. if it succeeds,
            # we break out of the loop to avoid the exception
            for wind in self.values():
                try:
                    tccd[twnam] = wind.crop(twind)
                    break
                except ValueError:
                    pass
            else:
                raise HipercamError(
                    "failed to find any enclosing window"
                    " for window label = {:s}".format(twnam)
                )

        return tccd

    @property
    def head(self):
        """Returns the first Window or an empty fits.Header if there is
        no first Window. The header of the first Window is where general
        header items of the CCD are stored."""
        return next(iter(self.values()), fits.Header())

    @head.setter
    def head(self, key, value):
        """Use to set general header items which will be stored in the
        first Window. If there is no first Window, a ValueError will be
        raised."""
        try:
            next(iter(self.values()))[key] = value
        except StopIteration:
            raise ValueError("there is no first Window to store any CCD header")

    def is_data(self):
        """Returns True / False according to whether the frame is thought
        to contain data. Uses DSTATUS keyword if present, else returns True
        """
        return len(self) and (self.head["DSTATUS"] if "DSTATUS" in self.head else True)

    def __repr__(self):
        return "{:s}(winds={:s}, nxtot={!r}, nytot={!r}, nxpad={!r}, nypad={!r})".format(
            self.__class__.__name__,
            super().__repr__(),
            self.nxtot,
            self.nytot,
            self.nxpad,
            self.nypad,
        )


class MCCD(Agroup):
    """Class representing a multi-CCD as a Group of CCD objects plus a FITS
    header. It supports arithematic operations and file I/O to/from FITS
    files.

    """

    def __init__(self, ccds, head=None, copy=False):
        """Constructs a :class:`MCCD`

        Arguments::

          ccds : Group
              Group of CCD objects.

          head : astropy.io.fits.Header
              a header which will be written as the primary header. If
              head=None on input, an empty header will be created.

          copy : bool
              if True, copy all the data over, otherwise only references are
              held. Holding references is fine if the calling program keeps
              re-generating the data but could cause problems in some
              circumstances.

        """
        super().__init__(CCD, ccds)
        if head is None:
            head = fits.Header()
        self.head = head
        if copy:
            self = self.copy()

    def __reduce__(self):
        """This is to overcome a problem with pickling MCCDs. The arguments for
        init don't get picked up for some reason so this must return them in
        the second element of the tuple
        """
        return (self.__class__, (list(self.items()), self.head))

    def write(self, fname, overwrite=False, xgap=200, ygap=200):
        """Writes out the MCCD to a FITS file.

        Arguments::

            fname : string or file-like object
                Name of file to write to. Can also be a file, opened in
                writeable binary mode.

            overwrite : bool
                True to overwrite pre-existing files

            xgap  : int
               X-gap used to space CCDs for ds9 mosaicing (unbinned pixels)

            ygap  : int
               Y-gap used to space CCDs for ds9 mosaicing (unbinned pixels)

        """

        phead = self.head.copy()
        phead["NUMCCD"] = (len(self), "Number of CCDs")
        phead["HIPERCAM"] = ("MCCD", "Type of HiPERCAM data (CCD | MCCD)")

        # Add comments if not already present.
        comm1 = "Data representing multiple CCDs written by hipercam.MCCD.write."
        if "COMMENT" not in phead or str(phead["COMMENT"]).find(comm1) == -1:
            phead.add_comment(comm1)
            phead.add_comment(
                "Each window of each CCD is written in a series of HDUs following an"
            )
            phead.add_comment(
                "HDU containing only the header. These follow an overall header for the"
            )
            phead.add_comment(
                "MCCD containing top-level information. The headers of the data window"
            )
            phead.add_comment(
                "HDUs have keywords LLX, LLY giving the pixel location in unbinned"
            )
            phead.add_comment(
                "pixels and XBIN and YBIN for their binning factors. The total unbinned"
            )
            phead.add_comment(
                "dimensions of each CCD are stored under keywords NXTOT and NYTOT for"
            )
            phead.add_comment(
                "each CCD. Each HDU associated with a given CCD is labelled with a"
            )
            phead.add_comment("header keyword CCD.")

        # make the first HDU
        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU(header=fits.Header(phead.cards)))

        # add in the HDUs of all the CCDs NX = 3 specific to HiPERCAM but
        # should do reasonably generally I think.
        xoff, yoff, noff = 0, 0, 0
        NX = 3
        for cnam, ccd in self.items():
            ccd.whdul(hdul, cnam, xoff, yoff)
            noff += 1
            if noff % NX == 0:
                xoff = 0
                yoff -= (ccd.nytot + 2 * ccd.nypad) + ygap
            else:
                xoff += (ccd.nxtot + 2 * ccd.nxpad) + xgap
        hdul.writeto(fname, overwrite=overwrite)

    @classmethod
    def read(cls, fname):
        """Builds an :class:`MCCD` from a FITS file. Expects a primary HDU, containing
        general headers and the first HDU, followed by a series of HDU blocks
        for each CCD with the first one of each containing some general
        per-CCD header items.

        Data are converted to float32 unless they are read in as float64.

        """
        with fits.open(fname) as hdul:
            return cls.rhdul(hdul)

    @classmethod
    def rhdul(cls, hdul):
        """Builds an :class:`MCCD` from an :class:`HDUList`.

        This will usually be passed the :class:`HDUList` from a
        file. The header from the first (primary) HDU will be used to
        create the header for the :class:`MCCD`. It is then assumed
        that the data for each :class:`CCD` is contained in blocks of
        HDUs after the first. The data from all HDUs will be read into
        the :class:`Window` objects that make up the CCD. Each HDU
        will be searched for a header parameter WINDOW to label the
        :class:`Window`, but the code will attempt to generate a
        sequential label if WINDOW is not found.

        Arguments::

          hdul : :class:`HDUList`
               each ImageHDU will be read as sub-window of the :class:`CCD`

        Data are converted to float32 unless they are read in as float64.

        """
        # Get the main header from the first hdu, but otherwise ignore it
        head = hdul[0].header
        if "NUMCCD" in head:
            del head["NUMCCD"]
        if "HIPERCAM" in head:
            del head["HIPERCAM"]

        # Attempt to read the rest of HDUs into a series of CCDs
        ccds = Group(CCD)
        for label, ccd in CCD.rmhdul(hdul[1:]):
            ccds[label] = ccd

        return cls(ccds, head)

    def copy(self, memo=None):
        """Make a copy of the :class:`MCCD`

        copy.copy and copy.deepcopy of an `MCCD` use this method
        """
        return MCCD(super().copy(memo), self.head.copy())

    def matches(self, mccd):
        """Check that the :class:`MCCD` matches another, which in this means checking
        that each CCD of the same label matches the equivalent in the other
        `MCCD`.

        There will be no match if there any CCDs present in one `MCCD` which are
        not in the other.

        Raises KeyError or ValueError exception depending on the problem.

        """
        for key in self:
            if key not in mccd:
                raise KeyError("CCD {0:d} not found in mccd".format(key))

        for key in mccd:
            if key not in self:
                raise KeyError("CCD {0:d} not found in self".format(key))

        for key, ccd in self.items():
            ccd.matches(mccd[key])

    def crop(self, mccd):
        """Crops the :class:MCCD to the same format as a template,
        `mccd`. Returns the new version. Will raise a HipercamError if
        it fails.

        """

        # build an empty CCD to get the headers. this will be added to
        tmccd = MCCD(Group(CCD), self.head)

        # wind through the CCDs of the template MCCD
        for cnam, ccd in mccd.items():
            if cnam in self:
                tmccd[cnam] = self[cnam].crop(ccd)

        return tmccd

    def __repr__(self):
        return "{:s}(ccds={:s}, head={!r})".format(
            self.__class__.__name__, super().__repr__(), self.head
        )

    def float32(self):
        """Applies :class:Window.float32 to all Windows of an MCCD"""
        for ccd in self.values():
            ccd.float32()

    def float64(self):
        """Applies :class:Window.float64 to all Windows of an MCCD"""
        for ccd in self.values():
            ccd.float64()

    def uint16(self):
        """Applies :class:Window.uint16 to all Windows of an MCCD"""
        for ccd in self.values():
            ccd.uint16()

    def set_const(self, val):
        """Sets all CCDs to a constant value"""
        for ccd in self.values():
            ccd.set_const(val)


def get_ccd_info(fname):
    """Routine to return some useful basic information from an MCCD file without
    reading the whole thing in. It returns an OrderedDict keyed on the CCD
    label which returns the maximum X and Y dimensions of the CCD and X and Y
    padding as a 4-element tuple. If no CCD label [header parameter 'CCD'] is
    found, the routine assigns a label of '1' on the assumption that the file
    cotains a single CCD.

    """
    info = OrderedDict()

    with fits.open(fname) as hdul:
        first = True
        for hdu in hdul:
            # The header of the HDU
            head = hdu.header

            if not first and (
                ("CCD" in head and ccd_label != head["CCD"])
                or ("NXTOT" in head and "NYTOT" in head)
            ):
                # Reset because we think we are on the first HDU of a CCD
                first = True

            if first and "NXTOT" in head and "NYTOT" in head:
                # Extract header from first HDU of a CCD
                if "CCD" in head:
                    ccd_label = head["CCD"]
                else:
                    ccd_label = "1"
                nxtot = head["NXTOT"]
                nytot = head["NYTOT"]
                nxpad = head.get("NXPAD", 0)
                nypad = head.get("NYPAD", 0)
                info[ccd_label] = (nxtot, nytot, nxpad, nypad)
                first = False

    return info


def trim_ultracam(mccd, ncol, nrow):
    """Trims columns and rows from an MCCD. The rows and columns are trimmed
    from each nearest to the readout amplifier as indicated by the value of
    the "outamp" attribute. e.g. if outamp='LL' then the rows and columns
    are stripped from the bottom and the left of the window respectively. This
    is useful because for ULTRACAM and ULTRASPEC these rows and columns can
    suffer from ringing effects.


    Arguments::

       mccd : MCCD
         the MCCD to be trimmed. Modified in place, i.e. it's changed on
         output.

       ncol : int
         number of columns to be trimmed.

       nrow : int
         number of rows to be trimmed

    """
    for cnam, ccd in mccd.items():
        for wlab, win in ccd.items():
            if win.outamp == "LL":
                # lower-left
                win.data = win.data[nrow:, ncol:]
                win.llx += ncol * win.xbin
                win.lly += nrow * win.ybin

            elif win.outamp == "LR":
                # lower-right
                if ncol:
                    win.data = win.data[nrow:, :-ncol]
                else:
                    win.data = win.data[nrow:, :]
                win.lly += nrow * win.ybin

            elif win.outamp == "UR":
                # upper-right
                if ncol and nrow:
                    win.data = win.data[:-nrow, :-ncol]
                elif nrow:
                    win.data = win.data[:-nrow, :]
                elif ncol:
                    win.data = win.data[:, :-ncol]

            elif win.outamp == "UL":
                # upper-left
                if nrow:
                    win.data = win.data[:-nrow, ncol:]
                else:
                    win.data = win.data[:, ncol:]
                win.llx += ncol * win.xbin

            else:
                warnings.warn(
                    "encountered a CCD window with no output"
                    " amplifier location defined"
                )
