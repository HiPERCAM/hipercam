"""hlog is a sub-module for reading in the log files written by reduce. It
defines a class to hold log files and a simple one to represent time series to
allow quick development of scripts to plot results.

For example, suppose a log file 'eg.log' has been written with 2 CCDs, '1' and
'2', and that CCD '2' has apertures labelled 't' and 'c' for target and
comparison. Then the following commands would load it, divide target by comparison
and plot the result with matplotlib:

  >> import matplotlib.pyplot as plt
  >> import hipercam as hcam
  >>
  >> hlog = hcam.hlog.Hlog.read('ts.log')
  >> targ = hlog.tseries('2','t')
  >> comp = hlog.tseries('2','c')
  >>
  >> ratio = targ / comp
  >> ratio.mplot(plt, 'r')
  >> plt.show()


"""

import struct
import warnings
import numpy as np
import copy
from astropy.io import fits
from astropy.stats import sigma_clip

from .core import *
from . import utils

__all__ = ("Hlog", "Tseries")

# maps numpy dtype code into struct flagd
NUMPY_TO_STRUCT = {
    "f4": "f",
    "f8": "d",
    "?": "?",
    "i4": "i",
    "u4": "I",
}

# Map of column names to output format for the Hlog.write method.
# Will need adapting if log files are changed. Designed to replicate
# the output in reduction.py

CNAME_TO_FMT = {
    "nframe": "{:d}",
    "MJD": "{:.14f}",
    "MJDok": "{:b}",
    "Exptim": "{:.8g}",
    "mfwhm": "{:.2f}",
    "mbeta": "{:.2f}",
    "x": "{:.4f}",
    "xe": "{:.4f}",
    "y": "{:.4f}",
    "ye": "{:.4f}",
    "fwhm": "{:.3f}",
    "fwhme": "{:.3f}",
    "beta": "{:.3f}",
    "betae": "{:.3f}",
    "counts": "{:.1f}",
    "countse": "{:.1f}",
    "sky": "{:.2f}",
    "skye": "{:.2f}",
    "nsky": "{:d}",
    "nrej": "{:d}",
    "flag": "{:d}",
}


class Hlog(dict):
    """Class to represent a HiPERCAM log as produced by reduce.  Based on
    dictionaries, Hlog files contain numpy structured arrays for each CCD
    which can be accessed by the CCD label name. Each array contains data that
    can be accessed by the label of the column. e.g.

       >> import hipercam as hcam
       >> hlog = hcam.hlog.Hlog.from_ascii('run011.log')
       >> print(hlog['2']['x_1'])

    would print the X-values of aperture '1' from CCD '2'.

    Hlog objects have four attributes:

     1) 'apnames', a dictionary keyed by CCD label giving a list of the 
        aperture labels used for each CCD. i.e.

        >> print(hlog.apnames['2'])

        might return ['1','2','3'].

    2) 'cnames', a dictionary keyed by CCD label giving a list of the column 
        names in the order they appeared in the original file. This is to
        help write out the data

    3) 'comments', a list of strings storing the header comments to
       allow the Hlog to be written out with full information if read
       from an ASCII log.

    4) 'writable', a flag to say whether the Hlog can be written which
       at the moment is only true if it has been read from an ASCII
       |hipercam| log file. Potentially fixable in the future.

    """

    @classmethod
    def read(cls, fname):
        warnings.warn("Hlog.read() is deprecated; use Hlog.rascii() instead.", DeprecationWarning)
        return Hlog.rascii(fname)

    @classmethod
    def rascii(cls, fname):
        """
        Loads a HiPERCAM ASCII log file written by reduce. Each CCD is loaded
        into a separate structured array, returned in a dictionary keyed by
        the CCD label.

        Argument::

           fname : string | list
              Can be a single log file name or a list of log file names. If a list,
              just the first file's comments will be read, the others ignored, and 
              it will be assumed, but not checked, that the later files have the
              same number and order of columns.

        """

        hlog = cls()
        read_cnames = False
        read_dtypes = False
        read_data = False
        dtype_defs = {}
        dtypes = {}
        struct_types = {}
        hlog.apnames = {}
        hlog.cnames = {}
        hlog.comments = []
        start = True

        if isinstance(fname, (list, tuple, np.ndarray)):
            fnames = fname
        else:
            fnames = [
                fname,
            ]

        first = True

        for fnam in fnames:
            with open(fnam) as fin:
                for line in fin:

                    if line.startswith("#"):

                        if first:
                            if start:
                                hlog.comments.append(line)

                            if line.find("Start of column name definitions") > -1:
                                read_cnames = True

                            elif line.find("Start of data type definitions") > -1:
                                read_dtypes = True

                            elif read_cnames:
                                # reading the column names
                                if line.find("End of column name definitions") > -1:
                                    read_cnames = False

                                elif line.find("=") > -1:
                                    # store the column names
                                    cnam = line[1 : line.find("=")].strip()
                                    hlog.cnames[cnam] = (
                                        line[line.find("=") + 1 :].strip().split()[1:]
                                    )
                                    hlog.apnames[cnam] = list(
                                        set(
                                            [
                                                item[item.find("_") + 1 :]
                                                for item in hlog.cnames[cnam]
                                                if item.find("_") > -1
                                            ]
                                        )
                                    )
                                    hlog.apnames[cnam].sort()

                            elif read_dtypes:
                                # reading the data types. We build
                                # numpy.dtype objects and strings for
                                # packing the data with struct.pack to
                                # save memory.
                                if line.find("End of data type definitions") > -1:
                                    read_dtypes = False

                                    # can now create the record array dtype
                                    for cnam in hlog.cnames:
                                        # numpy dtype
                                        dtypes[cnam] = np.dtype(
                                            list(
                                                zip(hlog.cnames[cnam], dtype_defs[cnam])
                                            )
                                        )
                                        # equivalent struct. Use
                                        # native byte order with no
                                        # alignment (initial '=') to
                                        # match numpy packing.
                                        struct_types[cnam] = "=" + "".join(
                                            [
                                                NUMPY_TO_STRUCT[dt]
                                                for dt in dtype_defs[cnam]
                                            ]
                                        )

                                    read_data = True

                                elif line.find("=") > -1:
                                    # store the data types
                                    cnam = line[1 : line.find("=")].strip()
                                    dtype_defs[cnam] = (
                                        line[line.find("=") + 1 :].strip().split()[1:]
                                    )

                                    # get ready for the data
                                    hlog[cnam] = bytearray()

                    elif read_data:

                        # No more storage of comment lines
                        start = False

                        # read a data line
                        items = line.strip().split()
                        cnam = items.pop(0)
                        dts = dtype_defs[cnam]

                        # convert types
                        for n in range(len(items)):
                            dt = dts[n]
                            if dt.startswith("f"):
                                items[n] = float(items[n])
                            elif dt.startswith("i") or dt.startswith("u"):
                                items[n] = int(items[n])
                            elif dt == "?":
                                items[n] = bool(int(items[n]))

                        # store in a list. Although lists are wasteful, they grow
                        # quite fast and each element here is efficiently packed
                        # so it should cope with quite large log files.
                        hlog[cnam].extend(struct.pack(struct_types[cnam], *items))

            # flag up subsequent runs
            first = False

        # for each CCD convert the list of byte data to a numpy array and then
        # ro a record array of the right dtype.
        for cnam in hlog:
            hlog[cnam] = np.frombuffer(hlog[cnam], dtypes[cnam])

        hlog.writable = True
        return hlog

    @classmethod
    def rfits(cls, fname):
        """
        Loads a HiPERCAM FITS log file written by |reduce| and subsequently
        converted by |hlog2fits|. Each CCD is loaded into a separate
        structured array, returned in a dictionary labelled by the CCD key.
        """

        hlog = cls()
        hlog.apnames = {}

        with fits.open(fname) as hdul:
            for hdu in hdul[1:]:
                cnam = hdu.header["CCDNAME"]
                hlog[cnam] = hdu.data
                hlog.apnames[cnam] = list(
                    set(
                        [
                            item[item.find("_") + 1 :]
                            for item in hdu.data.dtype.names
                            if item.find("_") > -1
                        ]
                    )
                )
                hlog.apnames[cnam].sort()

        hlog.writable = False
        return hlog

    @classmethod
    def rulog(cls, fname):
        """
        Creates an Hlog from an ASCII file produced by the C++ ULTRACAM
        pipeline.
        """

        hlog = cls()
        hlog.apnames = {}

        # CCD labels, number of apertures, numpy dtypes, struct types
        cnames = {}
        naps = {}
        dtypes = {}
        stypes = {}

        n = 0
        with open(fname) as fin:
            for line in fin:
                n += 1
                if not line.startswith("#"):
                    arr = line.split()
                    nframe, mjd, tflag, expose, cnam, fwhm, beta = arr[:7]

                    mjd = float(mjd)
                    tflag = bool(tflag)
                    expose = float(expose)
                    fwhm = float(fwhm)
                    beta = float(beta)

                    values = [mjd, tflag, expose, fwhm, beta]

                    if cnam in hlog:
                        # at least one line for this CCD has been read already
                        if len(arr[7:]) // 14 != naps[cnam]:
                            raise hcam.HipercamError(
                                (
                                    "First line of CCD {:s} had {:d} apertures,"
                                    " whereas line {:d} of file has {:d}"
                                ).format(cnam, naps[cnam], len(arr[7:]) // 14)
                            )

                        for nap in range(naps[cnam]):
                            (
                                naper,
                                x,
                                y,
                                xm,
                                ym,
                                exm,
                                eym,
                                counts,
                                sigma,
                                sky,
                                nsky,
                                nrej,
                                worst,
                                error_flag,
                            ) = arr[7 + 14 * nap : 7 + 14 * (nap + 1)]
                            values += [
                                float(x),
                                float(y),
                                float(xm),
                                float(ym),
                                float(exm),
                                float(eym),
                                float(counts),
                                float(sigma),
                                float(sky),
                                int(nsky),
                                int(nrej),
                                int(worst),
                                int(error_flag),
                            ]

                    else:
                        # first time for this CCD
                        hlog[cnam] = bytearray()
                        names = ["MJD", "MJDok", "Expose", "FWHM", "beta"]
                        dts = ["f8", "?", "f4", "f4", "f4"]
                        naps[cnam] = len(arr[7:]) // 14
                        hlog.apnames[cnam] = [str(i) for i in range(1, naps[cnam] + 1)]
                        for nap in range(naps[cnam]):
                            (
                                naper,
                                x,
                                y,
                                xm,
                                ym,
                                exm,
                                eym,
                                counts,
                                sigma,
                                sky,
                                nsky,
                                nrej,
                                worst,
                                error_flag,
                            ) = arr[7 + 14 * nap : 7 + 14 * (nap + 1)]
                            names += [
                                "x_{:s}".format(naper),
                                "y_{:s}".format(naper),
                                "xm_{:s}".format(naper),
                                "ym_{:s}".format(naper),
                                "exm_{:s}".format(naper),
                                "eym_{:s}".format(naper),
                                "counts_{:s}".format(naper),
                                "countse_{:s}".format(naper),
                                "sky_{:s}".format(naper),
                                "nsky_{:s}".format(naper),
                                "nrej_{:s}".format(naper),
                                "worst_{:s}".format(naper),
                                "flag_{:s}".format(naper),
                            ]
                            dts += [
                                "f4",
                                "f4",
                                "f4",
                                "f4",
                                "f4",
                                "f4",
                                "f4",
                                "f4",
                                "f4",
                                "i4",
                                "i4",
                                "i4",
                                "i4",
                            ]

                            values += [
                                float(x),
                                float(y),
                                float(xm),
                                float(ym),
                                float(exm),
                                float(eym),
                                float(counts),
                                float(sigma),
                                float(sky),
                                int(nsky),
                                int(nrej),
                                int(worst),
                                int(error_flag),
                            ]

                        dtypes[cnam] = np.dtype(list(zip(names, dts)))
                        stypes[cnam] = "=" + "".join(
                            [NUMPY_TO_STRUCT[dt] for dt in dts]
                        )

                    # store in a bytearray
                    hlog[cnam].extend(struct.pack(stypes[cnam], *values))

        # convert lists to numpy arrays
        for cnam in hlog:
            hlog[cnam] = np.frombuffer(hlog[cnam], dtype=dtypes[cnam])

        hlog.writable = False
        return hlog

    def tseries(self, cnam, apnam, name="counts"):
        """
        Returns with a Tseries corresponding to CCD cnam and
        aperture apnam. By default it accesses the 'counts',
        but 'x', 'y', 'fwhm', 'beta' and 'sky' are alternative
        choices.

        Arguments:

           cnam    : (string)
              CCD label. 'str' will be used to make a string of non-string
              entries.

           apnam   : (string)
              Aperture label. 'str' will be used to make a string of non-string
              entries.
        """
        ccd = self[str(cnam)]

        # Work out bad times which will be OR-ed onto the aperture
        # specific mask array
        mjdok = ccd["MJDok"]
        tmask = np.zeros_like(mjdok, dtype=np.uint)
        tmask[~mjdok] = BAD_TIME

        return Tseries(
            ccd["MJD"].copy(),
            ccd["{:s}_{:s}".format(name, str(apnam))],
            ccd["{:s}e_{!s}".format(name, apnam)],
            np.bitwise_or(ccd["flag_{!s}".format(apnam)], tmask),
        )

    def write(self, fname):
        """
        Writes out the Hlog to an ASCII file. This is to allow one to read in
        a log file, modify it and then write it out, useful for example to flag
        cloudy data. At the moment, it will only work for an Hlog read from a
        |hipercam| ASCII log file. NB It won't exactly replicate the log file
        input since it writes out in CCD order.
        """
        if self.writable:
            with open(fname, "w") as fout:
                # write out the comments
                for line in self.comments:
                    fout.write(line)

                # write out the data
                for cnam, data in self.items():
                    # generate format string
                    fmts = [
                        "{:s}",
                    ]
                    for name in data.dtype.names:
                        nam = name[: name.find("_")] if name.find("_") > -1 else name
                        fmts.append(CNAME_TO_FMT[nam])
                    fmt = " ".join(fmts) + "\n#\n"
                    for row in data:
                        fout.write(fmt.format(cnam, *row))
        else:
            raise Hipercam_Error("Hlog not writable")


class Tseries:
    """
    Class representing a basic time series with times, y values,
    y errors and flags. Attributes are::

       t    : (ndarray)
         times

       y    : (ndarray)
         y values

       ye   : (ndarray)
         y errors

       mask : (ndarray)
         bitmask propagated through from reduce log

    Errors are set = -1 and values = 0 if some sort of problem occurs. This
    should be used in addition to the mask to spot bad values. The bitmask
    valued are OR-ed when an operation on two Tseries is performed. Thus the
    child of two Tseries inherits all the problems of each of its parents.

    """

    def __init__(self, t, y, ye, mask):
        self.t = t
        self.y = y
        self.ye = ye
        self.mask = mask

    def __len__(self):
        return len(self.t)

    def set_mask(self, mvalue, mask):
        """
        This updates the internal mask by bitwise_or -ing it with the input mask value
        'mvalue' for all values for which mask is true
        """
        self.mask[mask] = np.bitwise_or(self.mask[mask], mvalue)

    def get_mask(self, mvalue=None, invert=False):
        """Returns logical array of all points with positive errors which
        also match the given bitmask mvalue. NB The array returned is not
        the same as the internal property 'mask'.

        Arguments::

           mvalue : bitmask
             use to select points according to the internal mask
             array. A match is made if any of the bits in the mask
             array value for a given point match the set bits in
             mvalue. e.g. mvalue=(NO_FWHM | NO_SKY) would select
             all points without a FWHM or sky aperture, while
             mvalue=ALL_OK selects perfect points only. Everything BUT
             the matching points can be selected using the 'invert'
             flag. See Tseries.report for how to get a full list of
             possible flags. A value of None makes no selection and
             just picks points with positive errors.

           invert : bool
             set True to return everything that does NOT match the supplied
             mvalue. Points with negative errors still ignored. Thus
             'mvalue=BAD_TIME, invert=True' will only plot points with good
             times.

        """
        if mvalue is not None:
            if mvalue == ALL_OK:
                # special case: no bits set at all.
                if invert:
                    return (self.mask != mvalue) & (self.ye > 0.0)
                else:
                    return (self.mask == mvalue) & (self.ye > 0.0)
            else:
                if invert:
                    return (np.bitwise_and(self.mask, mvalue) == 0) & (self.ye > 0.0)
                else:
                    return (np.bitwise_and(self.mask, mvalue) > 0) & (self.ye > 0.0)
        else:
            # no selection other than positive errors
            return self.ye > 0.0

    def get_data(self, mvalue=None, invert=False):
        """
        Returns the data as (times, values, errors). 'mask' is a bitmask
        to define which to select in addition to data with non-positive errors.

        Arguments::

           mvalue : bitmask
             bitmask to select points according to the internal mask array.

           invert : bool
             set True to return everything that does NOT match the supplied
             mvalue. Points with negative errors still ignored.
        """
        ok = self.get_mask(mvalue, invert)
        return (self.t[ok], self.y[ok], self.ye[ok])

    def mplot(
        self,
        axes,
        colour="b",
        fmt=".",
        mvalue=None,
        invert=False,
        capsize=0,
        erry=True,
        trange=None,
        **kwargs
    ):
        """Plots a Tseries to a matplotlib Axes instance, only plotting points
        that match the bitmask `mask` and have positive errors.

        Arguments::

           axes : Axes
              the axis instance to plot to. Can just be matplotlib.pyplot

           colour : valid matplotlib colour
              the colour to use

           fmt : string
              marker to use for points or style for line plots, e.g. ',',
              '.', 'o' give different points while '-' and '--' give
              different lines.

           mvalue : bitmask
              used to select the points to plot. See 'get_mask'

           invert : bool
             set True to return everything BUT the points selected
             by the mvalue

           capsize : float
              if error bars are plotted with points, this sets
              the length of terminals

           erry : boolean
              True / False for error bars or not

           trange : None | (t1,t2)
              Two element tuple to limit the time range

        """
        # Select points
        ok = self.get_mask(mvalue, invert)

        if trange is not None:
            t1, t2 = trange
            ok &= (self.t > t1) & (self.t < t2)

        if erry:
            axes.errorbar(
                self.t[ok],
                self.y[ok],
                self.ye[ok],
                fmt=fmt,
                color=colour,
                capsize=capsize,
                **kwargs
            )
        else:
            ymasked = np.ma.masked_array(self.y, ~ok)
            axes.plot(self.t, ymasked, fmt, color=colour, **kwargs)

    def __repr__(self):
        return "Tseries(t={!r}, y={!r}, ye={!r}, mask={!r}".format(
            self.t, self.y, self.ye, self.mask
        )

    def __truediv__(self, other):
        """Divides the Tseries by 'other' returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if len(self) != len(other):
                raise ValueError(
                    "input lengths [{:d} vs {:d}] do not match".format(
                        len(self), len(other)
                    )
                )

            if not all(self.t == other.t):
                raise ValueError("input times do not match")

            ok = (self.ye > 0) & (other.ye > 0)

            y = np.zeros_like(self.y)
            y[ok] = self.y[ok] / other.y[ok]

            ye = np.empty_like(y)
            ye[ok] = (
                np.sqrt(
                    self.ye[ok] ** 2 + (other.ye[ok] * self.y[ok] / other.y[ok]) ** 2
                )
                / other.y[ok]
            )
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # division by a constant or an array of constants
            y = self.y / other
            ye = np.empty_like(y)
            ok = self.ye > 0
            ye[ok] = self.ye[ok] / other
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t.copy(), y, ye, mask)

    def __mul__(self, other):
        """Multiplies the Tseries by 'other' returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if len(self) != len(other):
                raise ValueError(
                    "input lengths [{:d} vs {:d}] do not match".format(
                        len(self), len(other)
                    )
                )

            if not all(self.t == other.t):
                raise ValueError("input times do not match")

            ok = (self.ye > 0) & (other.ye > 0)

            y = np.zeros_like(self.y)
            y[ok] = self.y[ok] * other.y[ok]

            ye = np.empty_like(y)
            ye[ok] = np.sqrt(
                (other.y[ok] * self.ye[ok]) ** 2 + (self.ye[ok] * other.y[ok]) ** 2
            )
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # multiplication by a constant or an array of constants
            y = self.y * other
            ye = np.empty_like(y)
            ok = self.ye > 0
            ye[ok] = self.ye[ok] * other
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t.copy(), y, ye, mask)

    def __add__(self, other):
        """Add 'other' to the Tseries returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if len(self) != len(other):
                raise ValueError(
                    "input lengths [{:d} vs {:d}] do not match".format(
                        len(self), len(other)
                    )
                )

            if not all(self.t == other.t):
                raise ValueError("input times do not match")

            ok = (self.ye > 0) & (other.ye > 0)

            y = np.ones_like(self.y)
            y[ok] = self.y[ok] + other.y[ok]

            ye = np.empty_like(y)
            ye[ok] = np.sqrt(self.ye[ok] ** 2 + other.ye[ok] ** 2)
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # addition of a constant or an array of constants
            y = self.y + other
            ye = np.empty_like(y)
            ok = self.ye > 0
            ye[ok] = self.ye[ok]
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t.copy(), y, ye, mask)

    def __sub__(self, other):
        """Subtracts 'other' from the Tseries returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if len(self) != len(other):
                raise ValueError(
                    "input lengths [{:d} vs {:d}] do not match".format(
                        len(self), len(other)
                    )
                )

            if not all(self.t == other.t):
                raise ValueError("input times do not match")

            ok = (self.ye > 0) & (other.ye > 0)

            y = np.ones_like(self.y)
            y[ok] = self.y[ok] - other.y[ok]

            ye = np.empty_like(y)
            ye[ok] = np.sqrt(self.ye[ok] ** 2 + other.ye[ok] ** 2)
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # subtraction of a constant or an array of constants
            y = self.y - other
            ye = np.empty_like(y)
            ok = self.ye > 0
            ye[ok] = self.ye[ok]
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t.copy(), y, ye, mask)

    def __getitem__(self, key):
        copy_self = copy.copy(self)
        copy_self.t = self.t[key]
        copy_self.y = self.y[key]
        copy_self.ye = self.ye[key]
        copy_self.mask = self.mask[key]
        return copy_self

    def bin(self, binsize=10, method="mean"):
        """
        Bins the Timeseries using the function defined by `method` in blocks of binsize.

        Parameters
        -----------
        binsize : int
            Number of observations to incude in every bin
        method : str, one of 'mean' or 'median'
            The summary statistic to return

        Returns
        -------
        TSeries : Tseries object
            Binned Timeseries
        Notes
        -----
        - If the ratio between the Tseries length and the binsize is not
            a whole number, then the remainder of the data points will be
            ignored.
        - The binned TSeries will report the root-mean-square error.
        - The bitwise OR of the quality flags will be returned per bin.
        """
        available_methods = ["mean", "median"]
        if method not in available_methods:
            raise ValueError("method must be one of: {}".format(available_methods))
        methodf = np.__dict__["nan" + method]
        n_bins = len(self.t) // binsize

        t = np.array([methodf(a) for a in np.array_split(self.t, n_bins)])
        y = np.array([methodf(a) for a in np.array_split(self.y, n_bins)])
        ye = (
            np.array(
                [np.sqrt(np.nansum(a ** 2)) for a in np.array_split(self.ye, n_bins)]
            )
            / binsize
        )
        mask = np.array(
            [np.bitwise_or.reduce(a) for a in np.array_split(self.mask, n_bins)]
        )
        return Tseries(t, y, ye, mask)

    def remove_outliers(self, sigma=5.0, return_mask=False, **kwargs):
        """
        Removes outlier flux values using sigma-clipping.

        This method returns a new Tseries object from which flux values
        are removed if they are separated from the mean flux by `sigma` times
        the standard deviation.

        Parameters
        ----------
        sigma : float
            The number of standard deviations to use for clipping outliers.
            Defaults to 5.
        return_mask : bool
            Whether or not to return the mask indicating which data points
            were removed. Entries marked as `True` are considered outliers.
        **kwargs : dict
            Dictionary of arguments to be passed to `astropy.stats.sigma_clip`.

        Returns
        -------
        clean_tseries : Tseries object
            A new ``Tseries`` in which outliers have been removed.
        """
        outlier_mask = sigma_clip(data=self.y, sigma=sigma, **kwargs).mask
        if return_mask:
            return self[~outlier_mask], outlier_mask
        return self[~outlier_mask]

    def append(self, others):
        """
        Append Tseries objects

        Parameters
        ----------
        others : Tseries object or list of Tseries objects
            Time series to be appended to the current one

        Returns
        -------
        new_ts : Tseries object
            Concatenated time series
        """
        if not hasattr(others, "__iter__"):
            others = [others]
        new_lc = copy.copy(self)
        for i in range(len(others)):
            new_lc.t = np.append(new_lc.t, others[i].t)
            new_lc.y = np.append(new_lc.y, others[i].y)
            new_lc.ye = np.append(new_lc.ye, others[i].ye)
            new_lc.mask = np.append(new_lc.mask, others[i].mask)
        return new_lc

    def fold(self, period, t0=0.0):
        """Folds the time series at a specified ``period`` and ``phase``.

        This method returns a new ``Tseries`` object in which the time
        values range between -0.5 to +0.5.  Data points which occur exactly
        at ``t0`` or an integer multiple of `t0 + n*period` have time
        value 0.0.

        Parameters
        ----------
        period : float
            The period upon which to fold.
        t0 : float, optional
            Time reference point.

        Returns
        -------
        folded_tseries: Tseries object
            A new ``Tseries`` in which the data are folded and sorted by
            phase.
        """
        fold_time = ((self.t - t0) / period) % 1
        # fold time domain from -.5 to .5
        fold_time[fold_time > 0.5] -= 1
        sorted_args = np.argsort(fold_time)
        return Tseries(
            fold_time[sorted_args],
            self.y[sorted_args],
            self.ye[sorted_args],
            self.mask[sorted_args],
        )

    def normalise(self):
        """Returns a normalized version of the time series.

        The normalized timeseries is obtained by dividing `y` and `ye`
        by the median y value.

        Returns
        -------
        normalized_tseries : Tseries object
            A new ``Tseries`` in which `y` and `ye` are divided
            by the median y value.
        """
        lc = copy.copy(self)
        median_y = np.nanmedian(lc.y)
        lc.ye = lc.ye / median_y
        lc.y = lc.y / median_y
        return lc

    def downsize(self, other):
        """
        Bins the Timeseries down to match the times from 'other'.
        Useful for dealing with data taken with nskips

        Parameters::

          other : Tseries
            A coarser Tseries object you wish to match. Must have
            an integer ratio fewer points.

        Returns
        -------
          TSeries : Tseries object
             Binned Timeseries
        """

        if len(self) % len(other) != 0:
            raise ValueError(
                "length of other ({:d}) not a divisor of sel ({:d})".format(
                    len(other), len(self)
                )
            )

        nblock = len(other)

        y = self.y.copy()
        ye = self.ye.copy()
        bad = ye <= 0
        y[bad] = np.nan
        ye[bad] = np.nan

        y = np.array([np.nanmean(a) for a in np.array_split(y, nblock)])
        yec = np.array(
            [
                np.sqrt(np.nansum(a ** 2)) / max(1, np.count_nonzero(~np.isnan(a))) ** 2
                for a in np.array_split(ye, nblock)
            ]
        )
        bad = np.array(
            [np.count_nonzero(~np.isnan(a)) == 0 for a in np.array_split(ye, nblock)]
        )

        y[bad] = 0
        yec[bad] = -1

        mask = np.array(
            [np.bitwise_or.reduce(a) for a in np.array_split(self.mask, nblock)]
        )
        return Tseries(other.t, y, yec, mask)

    def ymean(self):
        ok = self.ye > 0.0
        return np.mean(self.y[ok])

    def clip_ends(self, nstart, nend):
        """
        Clips points off either end of Tseries, returns a new Tseries,
        leaves old one unchanged.

        Arguments::

           nstart : int
              number to remove from start

           nstart : int
              number to remove from end
        """
        if nstart < 0 or nend < 0:
            raise ValueError("nstart and nend must be >= 0")
        if nend:
            return Tseries(
                self.t[nstart:-nend],
                self.y[nstart:-nend],
                self.ye[nstart:-nend],
                self.mask[nstart:-nend],
            )
        elif nstart:
            return Tseries(
                self.t[nstart:], self.y[nstart:], self.ye[nstart:], self.mask[nstart:]
            )
        else:
            return copy.copy(self)

    def set_clouds(self, nwin, vmax, tmin):
        """This tries automatically to flag points as being affected by clouds using
        two criteria: (1) is a measure of the fractional RMS measured in a
        moving window of width nwin, (2) is a minimum transmission limit
        defined relative to the maximum value of the Tseries.

        It modifies the mask of the Tseries by bitwise OR-ing it with the flag
        hipercam.CLOUDS. No substitute for doing it by hand mind you as it
        tends to spread the period of clouds out to neighbouring points
        rather. It's main use is probably to clean up quick look plots.

        Arguments::

           nwin : int
              moving average window width. Must be odd

           vlim : float
              maximum relative variability (RMS/mean) above
              which a point will be flagged as cloud.

           tmin : float
              minimum transmission below which a point will
              be flagged as cloud.

        Modifies Tseries in place; returns three element tuple
        (clouds,rvar,trans) containing the logical cloud mask, the running
        relative variability and the running transmission. The "clouds"
        array can be used to transfer the clouds to other Tseries of
        the same length.

        For ULTRACAM / HIPERCAM / ULTRASPEC data it is often advisable
        to clip a point off the start and end to avoid end effects before
        applying this routine.

        """
        if nwin % 2 == 0:
            raise ValueError("nwin must be odd")
        WINDOW = np.ones(nwin) / nwin

        # extend ends with end values to get around end
        # effects
        y = np.concatenate(
            [self.y[0] * np.ones(nwin // 2), self.y, self.y[-1] * np.ones(nwin // 2)]
        )

        # apply running mean
        y_mave = np.convolve(y, WINDOW, "valid")

        # form squared differences
        ydsq = (self.y - y_mave) ** 2

        # extend at ends again
        ydsq = np.concatenate(
            [ydsq[0] * np.ones(nwin // 2), ydsq, ydsq[-1] * np.ones(nwin // 2)]
        )

        # compute running relative variability array
        rvar = np.sqrt(np.convolve(ydsq, WINDOW, "valid")) / y_mave
        trans = self.y / self.y.max()
        clouds = (rvar > vmax) | (trans < tmin)

        # modify Tseries mask
        self.set_mask(CLOUDS, clouds)

        return (clouds, rvar, trans)

    def report(self):
        """Reports numbers and types of bad points"""

        for fname, flag in FLAGS:
            match = np.bitwise_and(self.mask, flag) == flag
            print(
                "Flag = {:s} raised in {:d} points out of {:d}".format(
                    fname, int(round(match.sum())), len(self)
                )
            )

        print(
            "There were {:d} non-positive y-errors".format(len(self.ye[self.ye <= 0]))
        )

    def to_mag(self):
        """
        Convert to magnitudes, i.e. take -2.5*log10(y)
        """
        ok = (self.y > 0) & (self.ye > 0)
        new_ts = copy.copy(self)
        new_ts.y[ok] = -2.5 * np.log10(self.y[ok])
        new_ts.y[~ok] = 0
        new_ts.ye[ok] = 2.5 / np.log(10) * self.ye[ok] / self.y[ok]
        new_ts.ye[~ok] = -1
        return new_ts

    def write(self, fname, fmt="%.14e %.6e %.2e  %d", **kwargs):
        """
        Writes out the Tseries to an ASCII file using numpy.savetxt
        Other arguments to savetxt can be passed via kwargs
        """
        header = (
            kwargs.get("header", "")
            + """

Data written by hipercam.hlog.Tseries.write.
Four columns:
times    y-values  y-errors  integer mask
"""
        )
        kwargs["header"] = header
        np.savetxt(
            fname,
            np.column_stack([self.t, self.y, self.ye, self.mask]),
            fmt=fmt,
            **kwargs
        )

    @classmethod
    def read(cls, fname):
        """
        Creates a Tseries from an ASCII column data file
        """
        t, y, ye, mask = np.loadtxt(fname, unpack=True)
        return Tseries(t, y, ye, mask)
