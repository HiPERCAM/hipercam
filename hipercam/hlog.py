"""hlog is a sub-module for reading in the log files written by reduce. It
defines one class `Hlog` to hold log files and another `Tseries` to represent
time series to allow quick development of scripts to plot results.

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

The `Tseries` object know about bad data and carray a bitmask array reflecting problems
flagged during reduction.
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

NaN = float('NaN')

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

        # For backwards compatibility, try to spot bad data in old hipercam logs where
        # negative errors indicated bad data rather than NaNs. 'fixes' indicate the items
        # to look for and then the items to correct as a result
        fixes = {

            'perccd' : {
                'mfwhm' : ('mfwhm',),
                'mbeta' : ('mbeta',),
            },

            'peraper' : {
                'xe' : ('x', 'xe'),
                'ye' : ('y', 'ye'),
                'fwhme' : ('fwhm', 'fwhme'),
                'betae' : ('beta', 'beta'),
                'countse' : ('counts','countse'),
                'skye' : ('sky','skye'),
            },
        }

        for cnam in hlog.cnames:
            ccd = hlog[cnam]

            # per CCD corrections
            for key, corrs in fixes['perccd'].items():
                bad = ccd[key] <= 0.
                for corr in corrs:
                    ccd[corr][bad] = NaN

            # per aperture corrections
            for key, corrs in fixes['peraper'].items():
                for apnam in hlog.apnames[cnam]:
                    bad = ccd[f"{key}_{apnam}"] <= 0.
                    for corr in corrs:
                        ccd[f"{corr}_{apnam}"][bad] = NaN

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
                                float(xm) if float(exm) > 0 else NaN,
                                float(ym) if float(eym) > 0 else NaN,
                                float(exm) if float(exm) > 0 else NaN,
                                float(eym) if float(eym) > 0 else NaN,
                                float(counts) if float(sigma) > 0 else NaN,
                                float(sigma) if float(sigma) > 0 else NaN,
                                float(sky) if int(ssky) > 0 else NaN,
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
                                float(xm) if float(exm) > 0 else NaN,
                                float(ym) if float(eym) > 0 else NaN,
                                float(exm) if float(exm) > 0 else NaN,
                                float(eym) if float(eym) > 0 else NaN,
                                float(counts) if float(sigma) > 0 else NaN,
                                float(sigma) if float(sigma) > 0 else NaN,
                                float(sky) if int(ssky) > 0 else NaN,
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
    """Class representing a basic time series with times, y values, y
    errors and flags, and allowing bad data. Attributes are::

       t : ndarray
         times

       y : ndarray
         y values

       ye : ndarray
         y errors

       mask : ndarray
         bitmask propagated through from reduce log, that indicates possible
         problems affecting each data point.

    Bad data is indicated by any/all of t,y,ye being set = NaN. Errors
    are set = -1 and values = 0 if some sort of problem occurs. This
    should be used in addition to the mask to spot bad values. The
    bitmask valued are OR-ed when an operation on two Tseries is
    performed. Thus the child of two Tseries inherits all the problems
    of each of its parents.

    """

    def __init__(self, t, y, ye, mask):
        self.t = t
        self.y = y
        self.ye = ye
        self.mask = mask

    def __len__(self):
        return len(self.t)

    def set_mask(self, bitmask, mask):
        """
        This updates the internal mask by 'bitwise_or'-ing it with the input mask value
        'bitmask' for all values for which 'mask' is true
        """
        self.mask[mask] = np.bitwise_or(self.mask[mask], bitmask)

    def get_bad(self):
        """Returns with a numpy boolean array of bad data defined as data
        where any one of t, y, ye has been set to NaN.
        """
        return np.isnan(self.t) | np.isnan(self.y) | np.isnan(self.ye)

    def get_mask(self, bitmask=None, flag_bad=True):
        """Returns logical array of all points which contain NaN values or
        which match the bitmask 'bitmask'. The intention is to return
        a mask suitable for creating numpy.ma.masked_array
        objects. i.e. bitmask=(hcam.NOSKY | hcam.TARGET_SATURATED)
        would pick up all points said to have no sky or saturated
        target pixels. Thus the array returned flags bad data.

        Arguments::

           bitmask : None | numpy.uint32
              32-bit bitmask to select points according to the internal mask
              array.  e.g. mvalue=(NO_FWHM | NO_SKY) would select
              all points without a FWHM or sky aperture.
              See Tseries.report for how to get a full list of
              possible flags. A value of None makes no selection and only
              NaN data will be flagged.

           flag_bad : bool
              The usual operating mode (flag_bad=True) is to flag up any bad
              data (affected by NaNs) as well as data matching 'bitmask'. If
              set False, only data matching 'bitmask' is flagged.

        Returns::

          numpy logical array of same length as the Tseries.

        .. Note::

           The logical array returned is not the same as the internal
           attribute 'mask' which contains the bitmask itself.

        .. Note::

           Because of the way the mask is applied using
           'numpy.bitwise_and', it makes no sense to set
           bitmask=ALL_OK (it would actually be equivalent to setting
           bitmask=None in terms of operation) and so to do so raises
           a ValueErrror. In this instance you may want to set it to
           ANY_FLAG instead.

        """

        if bitmask == ALL_OK:
            raise ValueError('bitmask=ALL_OK is invalid; ANY_FLAG may be what you want')

        if flag_bad:
            # Flag bad data
            bad = self.get_bad()
        else:
            # Do not flag bad data
            bad = np.zeros_like(self.t, dtype=np.bool)

        if bitmask is not None:
            # Flag data matching the bitmask
            bad |= np.bitwise_and(self.mask, bitmask) > 0

        return bad

    def get_data(self, bitmask=None, flagged=False):
        """
        Returns the data as (times, values, errors). 'mask' is a bitmask
        to define which the data to exclude in addition to NaN values.

        Arguments::

           bitmask : int | None
             bitmask to identify points flagged in the internal mask array. Use to
             exclude bad data (in addition to data with NaN). See get_mask for usage.

           flagged : bool
             If set True, rather than using 'bitmask' to exclude data, this returns
             data that does match 'bitmask' (but still avoids bad NaN data)
        """
        if flagged:
            bad = self.get_mask(bitmask,False)
            return (self.t[bad], self.y[bad], self.ye[bad])
        else:
            good = ~self.get_mask(bitmask)
            return (self.t[good], self.y[good], self.ye[good])


    def mplot(
        self,
        axes,
        colour="b",
        fmt=".",
        bitmask=None,
        flagged=False,
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

           bitmask : None | int
              bitmask to remove bad points. See 'get_mask' for usage.

           flagged : bool
              set True to plot points that match the bitmask rather than
              ones that don't.

           capsize : float
              if error bars are plotted with points, this sets
              the length of terminals

           erry : boolean
              True / False for error bars or not

           trange : None | (t1,t2)
              Two element tuple to limit the time range

        """
        # Generate a mask to apply
        if flagged:
            mask = ~self.get_mask(bitmask)
        else:
            mask = self.get_mask(bitmask)

        if trange is not None:
            # add in time range limits
            t1, t2 = trange
            mask &= (self.t < t1) | (self.t > t2)

        # create masked arrays
        t = np.ma.masked_array(self.t,mask)
        y = np.ma.masked_array(self.y,mask)
        ye = np.ma.masked_array(self.ye,mask)

        if erry:
            axes.errorbar(
                t, y, ye,
                fmt=fmt,
                color=colour,
                capsize=capsize,
                **kwargs
            )
        else:
            axes.plot(t, y, fmt, color=colour, **kwargs)

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

            y = self.y / other.y
            ye = np.sqrt(
                self.ye**2 + (other.ye*self.y/other.y)**2
            ) / np.abs(other.y)

            mask = np.bitwise_or(self.mask, other.mask)

        else:

            # division by a constant or an array of constants
            y = self.y / other
            ye = self.ye / np.abs(other)
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

            y = self.y * other.y
            ye = np.sqrt(
                (other.y * self.ye) ** 2 + (self.ye * other.y) ** 2
            )
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # multiplication by a constant or an array of constants
            y = self.y * other
            ye = self.ye * np.abs(other)
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

            y = self.y + other.y

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

            y = self.y - other.y
            ye = np.sqrt(self.ye**2 + other.ye**2)
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # subtraction of a constant or an array of constants
            y = self.y - other
            ye = self.ye.copy()
            mask = self.mask.copy()

        return Tseries(self.t.copy(), y, ye, mask)

    def __getitem__(self, key):
        copy_self = copy.copy(self)
        copy_self.t = self.t[key]
        copy_self.y = self.y[key]
        copy_self.ye = self.ye[key]
        copy_self.mask = self.mask[key]
        return copy_self

    def bin(self, binsize, bitmask=None):
        """
        Bins the Timeseries into blocks of binsize; bitmask mvalue
        can be used to skip points.

        Parameters
        -----------
        binsize : int
            Number of observations to incude in every bin

        bitmask : int | None
            Bitmask that selects elements to ignore before binning, in addition to
            bad points. e.g. bitmask=hcam.TARGET_SATURATED will skip points
            flagged as saturated.

        Returns
        -------
        TSeries : Tseries object
            Binned Timeseries. If any output bin has no allowed input points (all
            bad or flagged by bitmask

        Notes
        -----
        - If the ratio between the Tseries length and the binsize is not
            a whole number, then the remainder of the data points will be
            ignored.
        - The binned TSeries will report the root-mean-square error.
        - The bitwise OR of the quality flags will be returned per bin.
        - Any bins with no points will have their errors set negative to mask them
        """

        n_bins = len(self) // binsize
        n_used = n_bins*binsize

        # Turn into numpy masked arrays, reshape to allow operations over x axis
        mask = self.get_mask(bitmask)[:n_used].reshape(n_bins, binsize)

        # compute number of points per bin
        ones = np.ma.masked_array(np.ones_like(mask), mask)
        nbin = np.ma.getdata(np.sum(ones, 1))

        t = np.ma.masked_array(self.t[:n_used].reshape(n_bins, binsize), mask)
        y = np.ma.masked_array(self.y[:n_used].reshape(n_bins, binsize), mask)
        ye = np.ma.masked_array(self.ye[:n_used].reshape(n_bins, binsize), mask)
        bmask = np.ma.masked_array(self.mask[:n_used].reshape(n_bins, binsize), mask)

        # take mean along x-axis, convert back to ordinary arrays
        t = np.ma.getdata(np.mean(t,1))
        y = np.ma.getdata(np.mean(y, 1))
        ye = np.ma.getdata(np.sqrt(np.sum(ye*ye,1)))
        ye[nbin == 0] = -1
        ye[nbin > 0] /= nbin[nbin > 0]
        bmask = np.ma.getdata(np.bitwise_or.accumulate(bmask, 1))

        # Return the binned Tseries
        return Tseries(t, y, ye, bmask)

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
        new_lc = copy.deepcopy(self)
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

    def normalise(self, bitmask=None, method='median'):
        """Returns a normalized version of the time series.

        The normalized timeseries is obtained by dividing `y` and `ye`
        by the median (or mean) y value. Arguments::

          bitmask : None | int
             see get_mask for meaning. Bad data are always ignored.

          method : str
             method to apply, either 'median' or 'mean'

        Returns
        -------
        normalized_tseries : Tseries object
            A new ``Tseries`` in which `y` and `ye` are divided
            by the median / mean y value.
        """
        lc = copy.deepcopy(self)
        ymask = np.ma.masked_array(lc.y, lc.get_mask(bitmask))
        if method == 'median':
            norm_factor = np.median(ymask)
        else:
            norm_factor = np.mean(ymask)
        lc.y /= norm_factor
        lc.ye /= norm_factor
        return lc

    def downsize(self, other, bitmask=None):
        """
        Bins the Tseries down to match the times from 'other'.
        Useful for dealing with data taken with nskips

        Parameters::

          other : Tseries
            A coarser Tseries you wish to match. Must have
            an integer ratio fewer points. e.g. for HiPERCAM this
            might be the u-band while 'self' is 'g'-band data, assuming
            they are commensurate.

          bitmask : None | int
            Bitmask to feed through to the binning routine used; see 'bin'.
            Use to try to exclude bad data where possible.

        Returns
        -------
          TSeries : Tseries
             Binned version of the Tseries. Every output bin will be present
             but some could be flagged as bad if there was no valid input. This
             allows the binned and coarse Tseries to stay 'in step'.
        """

        if (lself := len(self)) % (lother := len(other)) != 0:
            raise ValueError(
                f"length of other ({lother:d}) not a divisor of self ({lself:d})"
            )

        return self.bin(lself // lother, bitmask)

    def ymean(self, bitmask=None):
        ymask = np.ma.masked_array(lc.y, lc.get_mask(bitmask))
        return np.mean(ymask)

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
            return copy.deepcopy(self)

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
        """Reports numbers and types of bad points and details which are problematical."""

        ntot = len(self)
        bad = self.get_mask()
        nbad = len(bad[bad])
        flagged = np.bitwise_and(self.mask, ANY_FLAG) > 0
        nflag = len(flagged[flagged])

        print(f"There were {nbad} bad data points out of {ntot}")
        print(f"There were {nflag} flagged data points out of {ntot}")

        if nflag:
            print('\nFlags raised:\n')
            for fname, flag in FLAGS:
                if flag != ALL_OK and flag != ANY_FLAG:
                    match = np.bitwise_and(self.mask, flag) > 0
                    if (nfl := len(match[match])):
                        print(f"   Flag = {fname} was raised for {nfl} points out of {ntot}")

        if nbad or nflag:
            print('\nPoint-by-point (bad and/or flagged):\n')

            for i in range(len(self.t)):
                if bad[i] or flagged[i]:

                    if flagged[i]:
                        flags_raised = []
                        for fname, flag in FLAGS:
                            if flag != ANY_FLAG: 
                                if self.mask[i] & flag:
                                    flags_raised.append(fname)
                        flags_raised = ', '.join(flags_raised)
                    else:
                        flags_raised = "none"

                    print(f"   Index {i}: (t,y,ye) = ({self.t[i]},{self.y[i]},{self.ye[i]}), bad = {bad[i]}, flags raised = {flags_raised}")

    def to_mag(self):
        """
        Convert to magnitudes, i.e. take -2.5*log10(y)
        """
        new_ts = copy.deepcopy(self)
        new_ts.y = -2.5 * np.log10(self.y)
        new_ts.ye = 2.5 / np.log(10.) * self.ye / self.y
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
