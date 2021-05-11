"""hlog is a sub-module for reading in the log files written by reduce. It
defines one class `Hlog` to hold log files and another `Tseries` to represent
time series to allow quick development of scripts to plot results.

For example, suppose a log file 'eg.log' has been written with 2 CCDs,
'1' and '2', and that CCD '2' has apertures labelled 't' and 'c' for
target and comparison. Then the following commands would load it,
divide target by comparison and plot the result with matplotlib:

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

The `Tseries` object know about bad data and carray a bitmask array
reflecting problems flagged during reduction.

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
    "cmax": "{:d}",
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
                ok = ~np.isnan(ccd[key])
                bad = ccd[key][ok] <= 0.
                for corr in corrs:
                    ccd[corr][ok][bad] = NaN

            # per aperture corrections
            for key, corrs in fixes['peraper'].items():
                for apnam in hlog.apnames[cnam]:
                    ok = ~np.isnan(ccd[f"{key}_{apnam}"])
                    bad = ccd[f"{key}_{apnam}"][ok] <= 0.
                    for corr in corrs:
                        ccd[f"{corr}_{apnam}"][ok][bad] = NaN

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
                                float(sky) if int(nsky) > 0 else NaN,
                                int(nsky),
                                int(nrej),
                                int(worst),
                                int(error_flag),
                            ]

                    else:
                        # first time for this CCD
                        hlog[cnam] = bytearray()
                        names = ["MJD", "MJDok", "Exptim", "FWHM", "beta"]
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
                                "u4",
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
                                float(sky) if int(nsky) > 0 else NaN,
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

    def tseries(self, cnam, apnam, name="counts", ecol=True):
        """
        Returns with a Tseries corresponding to CCD cnam and
        aperture apnam. By default it accesses the 'counts',
        but 'x', 'y', 'fwhm', 'beta' and 'sky' are alternative
        choices. Looks at the column names in a HiPERCAM log.
        Can also access items not specific to apertures.

        Arguments:

           cnam : str
              CCD label. 'str' will be used to make a string of non-string
              entries.

           apnam : str
              Aperture label. Set = 'None' if the item of interest is not
              specific to apertures, e.g. 'mfwhm'.

           name : str
              Item to return. e.g. 'counts', 'x', 'fwhm', 'mfwhm'. If apnam is
              not None, then f"{name}_{apnam}" will be used as the column name,
              else f"{name}"

           ecol : bool
              If True, an attempt will be made to read errors from a column called
              f"{name}e_{apnam}" or f"{name}e" will be made. Otherwise the errors
              are set = 0.
        """
        ccd = self[str(cnam)]

        # Work out bad times which will be OR-ed onto the
        # aperture-specific bitmask array
        mjdok = ccd["MJDok"]
        tmask = np.zeros_like(mjdok, dtype=np.uint)
        tmask[~mjdok] = BAD_TIME

        times = ccd["MJD"].copy()
        texps = ccd["Exptim"].copy()/86400

        if apnam is None:
            data = ccd[f"{name}"].copy()
            if ecol:
                errors = ccd[f"{name}e"].copy()
            else:
                errors = np.zeros_like(mjdok)
            bmask = tmask
        else:
            data = ccd[f"{name}_{apnam}"].copy()
            if ecol:
                errors = ccd[f"{name}e_{apnam}"].copy()
            else:
                errors = np.zeros_like(mjdok)
            bmask = ccd[f"flag_{apnam}"] | tmask

        return Tseries(times, data, errors, bmask, texps)

    def write(self, fname):
        """Writes out the Hlog to an ASCII file. This is to allow one to read
        in a log file, modify it and then write it out, useful for
        example to flag cloudy data. At the moment, it will only work
        for an Hlog read from a |hipercam| ASCII log file. NB It won't
        exactly replicate the log file input since it writes out in
        CCD order.

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
         mid-times

       y : ndarray
         y values

       ye : ndarray
         y errors

       bmask : ndarray
         bitmask propagated through from reduce log, that indicates possible
         problems affecting each data point. See core.py for the defined flags.

       te : None : ndarray
         Exposure times (not much is done with these, but sometimes it's good
         to have them on output.

       cpy : bool
         If True, copy the arrays by value, not reference.

    Data can be bad, or it can be flagged, or it can be both. Some
    flags imply that the associated data will be bad, but there are
    some you may be able to live with. Bad data are indicated by
    any/all of t,y,ye being set = NaN. The bitmask values correspond
    to flags defined in hipercam.core. The bitmask values are OR-ed
    when an operation on two Tseries is performed. Thus the child of
    two Tseries inherits all the problems of each of its parents.

    """

    def __init__(self, t, y, ye=None, bmask=None, te=None, cpy=False):
        """
        ye=None and bmask=None means their arrays will be set=0
        """
        if isinstance(t, np.ndarray) and (
                len(t) != len(y) or \
                (ye is not None and len(t)!= len(ye)) or \
                (bmask is not None and len(t) != len(bmask)) or \
                (te is not None and len(t) != len(te))):
            raise ValueError("problem with one or more of t,y, ye, bmask, te")

        if cpy:
            self.t = copy.copy(t)
            self.y = copy.copy(y)
            self.ye = copy.copy(ye) if ye is not None else \
                np.zeros_like(t)
            self.bmask = copy.copy(bmask) if bmask is not None else \
                np.zeros_like(t,dtype=np.uint32)
            self.te = copy.copy(te)
        else:
            self.t = t
            self.y = y
            self.ye = ye if ye is not None else \
                np.zeros_like(t)
            self.bmask = bmask if bmask is not None else \
                np.zeros_like(t,dtype=np.uint32)
            self.te = te

    def __len__(self):
        return len(self.t)

    def set_bitmask(self, bitmask, mask):
        """This updates the internal bitmask by 'bitwise_or'-ing it with the
        input bitmask value 'bitmask' for all values for which 'mask'
        is true. 'mask' should match the Tseries length.  Arguments::

           bitmask : int
             a bitmask value which will be OR-ed with elements of the current bitmask array

           mask: np.ndarray
             the bitmask will be applied to every element for which mask is True.

        """
        self.bmask[mask] |= bitmask

    def get_bad(self):
        """Returns with a numpy boolean array of bad data defined as data
        where any one of t, y, ye has been set to NaN or Inf
        """
        return (
            np.isnan(self.t) | np.isnan(self.y) | np.isnan(self.ye) |
            np.isinf(self.t) | np.isinf(self.y) | np.isinf(self.ye)
        )

    def get_mask(self, bitmask=None, flag_bad=True):
        """Returns logical array of all points which contain NaN values or
        which match the bitmask 'bitmask' (i.e. bitwise_and > 0) The
        intention is to return a boolean mask suitable for creating
        numpy.ma.masked_array objects. i.e. bitmask=(hcam.NOSKY |
        hcam.TARGET_SATURATED) would pick up all points said to have
        no sky or saturated target pixels. Thus the array returned
        flags bad data.

        Arguments::

           bitmask : None | numpy.uint32
              32-bit bitmask to select points according to the internal bmask
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
           attribute 'bmask' which contains the bitmask array itself.

        .. Note::

           Because of the way the mask is applied using
           'numpy.bitwise_and', it makes no sense to set
           bitmask=ALL_OK (it would actually be equivalent to setting
           bitmask=None in terms of operation) and so to do so raises
           a ValueErrror. In this instance you may want to set it to
           ANY_FLAG instead to flag up points flagged for any reason.

        """

        if bitmask == ALL_OK:
            raise ValueError(
                'bitmask=ALL_OK is invalid; ANY_FLAG may be what you want'
            )

        if flag_bad:
            # Flag bad data
            bad = self.get_bad()
        else:
            # Do not flag bad data
            bad = np.zeros_like(self.t, dtype=np.bool)

        if bitmask is not None:
            # Flag data matching the bitmask
            bad |= (self.bmask & bitmask) > 0

        return bad

    def get_data(self, bitmask=None, flagged=False):
        """
        Returns the data as (times, exposures, values, errors). 'bitmask' is a bitmask
        to define which the data to exclude in addition to any NaN values.

        Arguments::

           bitmask : int | None
             bitmask to identify points flagged in the internal bmask array. Use to
             exclude bad data (in addition to data with NaN). See get_mask for usage.

           flagged : bool
             If set True, rather than using 'bitmask' to exclude data, this uses it
             to include data (but still avoids bad NaN data)

        Returns::
           data : tuple
             times, exposures times [or None if they are not defined], y values and their 
             uncertainties
        """
        if flagged:
            bad = self.get_mask(bitmask,False)
            return (self.t[bad], self.te if self.te is None else self.te[bad], self.y[bad], self.ye[bad])
        else:
            good = ~self.get_mask(bitmask)
            return (self.t[good], self.te if self.te is None else self.te[good], self.y[good], self.ye[good])

    def percentile(self, q, bitmask=None):
        """Returns percentiles of the y-array and their -/+ 1 sigma to allow
        for error bars.  q = a percentile or array of percentiles
        (values from 0 to 100) The return value is a tuple of three
        lists containing the percentiles for y, y-ye, y+ye.

        Arguments::

           q : float | array
             Percentile or list of percentiles to calculate, in range 0-100.

           bitmask : None | int
             see get_mask for meaning. Bad data are always ignored.

        Returns::

           Three lists containing percentiles for y, y-ye and y+ye
        """
        good = ~self.get_mask(bitmask)
        y = np.percentile(self.y[good],q)
        ym = np.percentile(self.y[good]-self.ye[good],q)
        yp = np.percentile(self.y[good]+self.ye[good],q)
        return (y,ym,yp)

    def mplot(
            self,
            axes,
            color="b",
            fmt=".",
            bitmask=None,
            flagged=False,
            capsize=0,
            errx=False,
            erry=True,
            trange=None,
            mask=None,
            **kwargs
    ):
        """Plots a Tseries to a matplotlib Axes instance, only plotting points
        that match the bitmask `mask` and have positive errors.

        Arguments::

           axes : Axes
              the axis instance to plot to. Can just be matplotlib.pyplot

           color : matplotlib colour
              the colour to use (fed into kwargs to establish a fixed default)

           fmt : str
              marker to use for points or style for line plots, e.g. ',',
              '.', 'o' give different points while '-' and '--' give
              different lines. (fed into kwargs)

           bitmask : None | int
              bitmask to remove bad points. See 'get_mask' for usage.

           flagged : bool
              set True to plot points that match the bitmask rather than
              ones that don't.

           capsize : float
              if error bars are plotted with points, this sets
              the length of terminals (fed into kwargs)

           errx : bool
              True / False for bars indicating exposure length (i.e. +/- 1/2
              whatever is in the te array)

           erry : bool
              True / False for vertical error bars or not

           trange : None | (t1,t2)
              Two element tuple to limit the time range

           mask : None | logical array
              True for point to be plotted

           kwargs : keyword arguments
              These will be fed to the plot routine which is either
              matplotlib.pyplot.errorbar or matplotlib.pyplot.plot.
              e.g. 'ms=2' will set the markersize to 2.
        """

        # Generate boolean array of the points to plot
        if flagged:
            plot = self.get_mask(bitmask)
        else:
            plot = ~self.get_mask(bitmask)
            
        if mask is not None:
            plot &= mask

        if trange is not None:
            # add in time range limits
            t1, t2 = trange
            plot &= (self.t > t1) & (self.t < t2)

        errx &= (self.te is not None) and np.any(self.te > 0)
        erry &= np.any(self.ye[~np.isnan(self.ye)] > 0)

        # create masked arrays
        t = self.t[plot]
        y = self.y[plot]

        kwargs["color"] = color

        if errx and erry:
            te = self.te[plot]/2.
            ye = self.ye[plot]
            kwargs["fmt"] = fmt
            kwargs["capsize"] = capsize
            axes.errorbar(t, y, ye, te, **kwargs)

        elif errx:
            te = self.te[plot]/2.
            kwargs["fmt"] = fmt
            kwargs["capsize"] = capsize
            axes.errorbar(t, y, xerr=te, **kwargs)

        elif erry:
            ye = self.ye[plot]
            kwargs["fmt"] = fmt
            kwargs["capsize"] = capsize
            axes.errorbar(t, y, ye, **kwargs)

        else:
            axes.plot(t, y, fmt, **kwargs)

    def __repr__(self):
        return f"Tseries(t={self.t!r}, y={self.y!r}, ye={self.ye!r}, bmask={self.bmask!r}, te={self.te!r}"

    def __truediv__(self, other):
        """Divides the Tseries by 'other' returning the result as another
        Tseries. Propagates errors where possible, but just takes one of the
        sets of exposure times (if present). 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together.
        """
        if isinstance(other, Tseries):
            
            with np.errstate(divide='ignore', invalid='ignore'):
                y = self.y / other.y
                ye = np.sqrt(
                    self.ye**2 + (self.y*other.ye/other.y)**2
                ) / np.abs(other.y)

            bmask = self.bmask | other.bmask
            te = self.te.copy() if self.te is not None else \
                other.te.copy() if other.te is not None else None

        else:

            # division by a constant or an array of constants
            y = self.y / other
            ye = self.ye / np.abs(other)
            bmask = self.bmask.copy()
            te = self.te.copy() if self.te is not None else None

        return Tseries(self.t.copy(), y, ye, bmask, te)

    def __itruediv__(self, other):
        """
        Divides the Tseries by 'other' in place. See __truediv__ for more details.
        """
        if isinstance(other, Tseries):
            
            with np.errstate(divide='ignore', invalid='ignore'):
                self.ye = np.sqrt(
                    self.ye**2 + (self.y*other.ye/other.y)**2
                ) / np.abs(other.y)

                self.y /= other.y
                
            self.bmask |= other.bmask
            if self.te is None and other.te is not None:
                self.te = other.te.copy()

        else:

            # division by a constant or an array of constants
            self.y /= other
            self.ye /= np.abs(other)

        return self

    def __mul__(self, other):
        """Multiplies the Tseries by 'other' returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. bad data also propagates,
        i.e. if either input on a given pixel is bad, then so too is the output
        """
        if isinstance(other, Tseries):

            y = self.y * other.y
            ye = np.sqrt(
                (other.y*self.ye)**2 + (self.y*other.ye)**2
            )
            bmask = self.bmask | other.bmask
            te = self.te.copy() if self.te is not None else \
                other.te.copy() if other.te is not None else None

        else:
            # multiplication by a constant or an array of constants
            y = self.y * other
            ye = self.ye * np.abs(other)
            bmask = self.bmask.copy()
            te = self.te.copy() if self.te is not None else None

        return Tseries(self.t.copy(), y, ye, bmask, te)

    def __imul__(self, other):
        """Multiplies the Tseries by 'other' in place. See __mul__ for more.
        """
        if isinstance(other, Tseries):

            self.ye = np.sqrt(
                (other.y*self.ye)**2 + (self.y*other.ye)**2
            )
            self.y *= other.y
            self.bmask |= other.bmask
            if self.te is None and other.te is not None:
                self.te = other.te.copy()

        else:
            # multiplication by a constant or an array of constants
            self.y *= other
            self.ye *= np.abs(other)

        return self

    def __add__(self, other):
        """Add 'other' to the Tseries returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together.
        """
        if isinstance(other, Tseries):

            y = self.y + other.y
            ye = np.sqrt(self.ye**2+ other.ye**2)
            bmask = self.bmask | other.bmask
            te = self.te.copy() if self.te is not None else \
                other.te.copy() if other.te is not None else None

        else:
            # addition of a constant or an array of constants
            y = self.y + other
            ye = self.ye.copy()
            bmask = self.bmask.copy()
            te = self.te.copy() if self.te is not None else None

        return Tseries(self.t.copy(), y, ye, bmask, te)

    def __iadd__(self, other):
        """
        Add 'other' to the Tseries in place. See __add_ for more details.
        """
        if isinstance(other, Tseries):

            self.y += other.y
            self.ye = np.sqrt(self.ye**2+ other.ye**2)
            self.bmask |= other.bmask
            if self.te is None and other.te is not None:
                self.te = other.te.copy()

        else:
            # addition of a constant or an array of constants
            self.y += other

        return self

    def __sub__(self, other):
        """Subtracts 'other' from the Tseries returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together.
        """
        if isinstance(other, Tseries):

            y = self.y - other.y
            ye = np.sqrt(self.ye**2 + other.ye**2)
            bmask = self.bmask | other.bmask
            te = self.te.copy() if self.te is not None else \
                other.te.copy() if other.te is not None else None

        else:
            # subtraction of a constant or an array of constants
            y = self.y - other
            ye = self.ye.copy()
            bmask = self.bmask.copy()
            te = self.te.copy() if self.te is not None else None

        return Tseries(self.t.copy(), y, ye, bmask, te)

    def __isub__(self, other):
        """
        Subtracts 'other' from the Tseries in place. See __sub__ for more.
        """
        if isinstance(other, Tseries):

            self.y -= other.y
            self.ye = np.sqrt(self.ye**2 + other.ye**2)
            self.bmask |= other.bmask
            if self.te is None and other.te is not None:
                self.te = other.te.copy()

        else:
            # subtraction of a constant or an array of constants
            self.y -= other

        return self

    def __getitem__(self, key):
        copy_self = copy.deepcopy(self)
        copy_self.t = self.t[key]
        copy_self.y = self.y[key]
        copy_self.ye = self.ye[key]
        copy_self.bmask = self.bmask[key]
        copy_self.te = self.te[key] if self.te is not None else None
        return copy_self

    def tadd(self, other):
        """
        Adds 'other' to the times, in place.
        """
        self.t += other

    def tsub(self, other):
        """
        Subtracts 'other' to the times, in place.
        """
        self.t -= other

    def tmul(self, other):
        """
        Multiplies the times by 'other', in place.
        """
        self.t *= other
        if self.te is not None:
            self.te *= np.abs(other)

    def tdiv(self, other):
        """
        Divides the times by 'other', in place.
        """
        self.t /= other
        if self.te is not None:
            self.te /= np.abs(other)

    def bin(self, binsize, bitmask=None, inplace=True):
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

        inplace : bool
            If True, self is modified, else a new Tseries is created
            and self is untouched. A reference to a Tseries is always returned

        Returns
        -------
        TSeries : Tseries object
            Binned Timeseries. If any output bin has no allowed input points (all
            bad or flagged by bitmask

        .. Notes::

           (1) If the ratio between the Tseries length and the binsize is not
           a whole number, then the remainder of the data points will be
           ignored.

           (2) The binned TSeries will report the root-mean-square error.

           (3) The bitwise OR of the quality flags will be returned per bin.

           (4) Any bins with no points will have their data set bad to mask them

           (5) The exposure time will be set to span the first to last time contributing
               to the bin.
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
        bmask = self.bmask[:n_used].reshape(n_bins, binsize)
        bmask[mask] = 0

        # Define minimum and maximum times to cope with the exposure times
        if self.te is None:
            tmin = t
            tmax = t
        else:
            te = np.ma.masked_array(self.te[:n_used].reshape(n_bins, binsize), mask)/2
            tmin = t - te
            tmax = t + te

        # take mean / sum along x-axis of 2D reshaped arrays, convert back to ordinary arrays
        t = np.ma.getdata(np.mean(t,1))
        y = np.ma.getdata(np.mean(y, 1))
        ye = np.ma.getdata(np.sqrt(np.sum(ye*ye,1)))

        # note this next bit is not quite right since the time had
        # been set to the mean, but it will do
        tmin = np.ma.getdata(np.min(tmin,1))
        tmax = np.ma.getdata(np.max(tmax,1))
        te = tmax-tmin
        y[nbin == 0] = NaN
        ye[nbin == 0] = NaN
        ye[nbin > 0] /= nbin[nbin > 0]
        bmask = np.bitwise_or.reduce(bmask, 1)
        if inplace:
            self.t = t
            self.y = y
            self.ye = ye
            self.bmask = bmask
            self.te = te
            return self
        else:
            # Return the binned Tseries
            return Tseries(t, y, ye, bmask, te)

    def flag_outliers(self, sigma=4.0, inplace=True, **kwargs):
        """
        Flags outlier flux values using sigma-clipping according to
        `astropy.stats.sigma_clip` [see docs on that for details].

        Parameters
        ----------
        sigma : float
            The number of standard deviations to use for clipping outliers.
            Defaults to 5.

        inplace : bool
            Default action is to apply the mask in place. Otherwise returns
            a new Tseries, leaving the current one untouched.

        **kwargs : dict
            Dictionary of arguments to be passed to `astropy.stats.sigma_clip`.

        Returns
        -------
        clean_tseries : Tseries object [if not inplace]
            Modified Tseries
        """

        # identify outliers
        bad = self.get_bad()
        ymask = np.ma.masked_array(self.y, mask=bad)
        mask = sigma_clip(data=ymask, sigma=sigma, **kwargs).mask

        # remove any NaNs or Infs from this as they do not
        # count as "outliers"
        mask &= ~bad

        if inplace:
            self.set_bitmask(OUTLIER, mask)
        else:
            new_ts = copy.deepcopy(self)
            new_ts.set_bitmask(OUTLIER, mask)
            return new_ts

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

        .. Note::

           To end with any exposure times, all contributing series
           will need to have exposure times defined.
        """
        if not hasattr(others, "__iter__"):
            others = [others]

        new_lc = copy.deepcopy(self)
        edef = new_lc.te is not None
        for i in range(len(others)):
            new_lc.t = np.append(new_lc.t, others[i].t)
            new_lc.y = np.append(new_lc.y, others[i].y)
            new_lc.ye = np.append(new_lc.ye, others[i].ye)
            new_lc.bmask = np.append(new_lc.bmask, others[i].bmask)
            edef &= others[i].te is not None
            if edef:
                new_lc.te = np.append(new_lc.te, others[i].te)
            else:
                new_lc.te = None

        return new_lc

    def ttrans(self, func, efunc=None, inplace=True):
        """
        Applies "func" to transform a time scale. "func" should return
        an array of transformed times given an array of raw times. "efunc"
        if defined will transform the exposure times (if they exist).

        For instance:

           >>> ts.ttrans(lambda t : 1440*(t-T0), lambda te : 1440*te)

        would convert times in days to minutes offset from T0.
        """

        if inplace:
            self.t = func(self.t)
            if efunc is not None and self.te is not None:
                self.te = efunc(self.te)
        else:
            new_ts = copy.deepcopy(self)
            new_ts.t = func(new_ts.t)
            if efunc is not None and self.te is not None:
                new_ts.te = efunc(new_ts.te)
            return new_ts

    def phase(self, t0, period, fold=False, inplace=True):
        """Convert the time into phase.

        This method returns a new ``Tseries`` object in which the time
        values range between -0.5 to +0.5.  Data points which occur exactly
        at ``t0`` or an integer multiple of `t0 + n*period` have time
        value 0.0.

        Parameters
        ----------

        t0 : float
            Time reference point.

        period : float
            The period

        fold : bool
            If True, map the phases to -0.5 to 0.5, and sort them.

        inplace : bool
            Operation in place, else return a new Tseries

        Returns
        -------
        folded_tseries: Tseries object
            A new ``Tseries`` with the time converted to phase
            (if inplace=False)
        """
        phase = (self.t - t0) / period
        pexpose = self.te/np.abs(period) if self.te is not None else None

        if fold:
            phase = np.mod(phase, 1)
            phase[phase > 0.5] -= 1
            isort = np.argsort(phase)

            if inplace:
                self.t = phase[isort]
                self.te = pexpose[isort] if pexpose is not None else None
                self.y = self.y[isort]
                self.ye = self.ye[isort]
                self.bmask = self.bmask[isort]
            else:
                return Tseries(
                    phase[isort],
                    self.y[isort], self.ye[isort], self.bmask[isort],
                    pexpose[isort] if pexpose is not None else None,
                    True
                )

        else:
            if inplace:
                self.t = phase
                self.te = pexpose
            else:
                return Tseries(
                    phase, self.y, self.ye, self.bmask, pexpose, True
                )

    def normalise(self, bitmask=None, method='median', weighted=False, inplace=True):
        """Returns a normalized version of the time series.

        The normalized timeseries is obtained by dividing `y` and `ye`
        by the median (or mean) y value. Arguments::

          bitmask : None | int
             see get_mask for meaning. Bad data are always ignored.

          method : str
             method to apply, either 'median' or 'mean'

          weighted : bool
             if method == 'mean', compute an inverse-variance weighted mean

          inplace: bool
             if True, self is modfied, else a new Tseries is created. A
             reference to a Tseries is always returned.

        Returns
        -------
        normalized_tseries : Tseries object
            A new ``Tseries`` in which `y` and `ye` are divided
            by the median / mean y value.
        """
        if inplace:
            lc = self
        else:
            lc = copy.deepcopy(self)
        mask = lc.get_mask(bitmask)
        ymask = np.ma.masked_array(lc.y, mask)
        if method == 'median':
            norm_factor = np.ma.median(ymask)
        elif method == 'mean':
            if weighted:
                wgt = 1/np.ma.masked_array(lc.ye, mask)**2
                norm_factor = np.sum(wgt*ymask)/np.sum(wgt)
            else:
                norm_factor = np.mean(ymask)
        lc.y /= norm_factor
        lc.ye /= norm_factor
        return lc

    def downsize(self, other, bitmask=None, inplace=True):
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

          inplace : bool
            If True, overwrite the Tseries, else just return a new one.
            A reference to a Tseries is always returned

        Returns
        -------
          TSeries : Tseries
             Binned version of the Tseries. Every output bin will be present
             but some could be flagged as bad if there was no valid input. This
             allows the binned and coarse Tseries to stay 'in step'.
        """

        lself, lother = len(self), len(other)
        if lself % lother != 0:
            raise ValueError(
                f"length of other ({lother:d}) not a divisor of self ({lself:d})"
            )

        return self.bin(lself // lother, bitmask, inplace=inplace)

    def ymean(self, bitmask=None, weighted=True):
        """
        Computes the mean of the y-values. "weighted=True"
        means the mean is calculated with inverse variance
        weighting.
        """
        mask = self.get_mask(bitmask)
        ymask = np.ma.masked_array(self.y,mask)
        if weighted:
            wgt = 1/np.ma.masked_array(self.ye, mask)**2
            return np.sum(wgt*ymask)/np.sum(wgt)
        else:
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
            te = self.te[nstart:-nend] if self.te is not None else None
            return Tseries(
                self.t[nstart:-nend],
                self.y[nstart:-nend],
                self.ye[nstart:-nend],
                self.bmask[nstart:-nend],
                te
            )
        elif nstart:
            te = self.te[nstart:] if self.te is not None else None
            return Tseries(
                self.t[nstart:], self.y[nstart:], self.ye[nstart:],
                self.bmask[nstart:], te
            )
        else:
            return copy.deepcopy(self)

    def report(self):
        """Reports numbers and types of bad points and details which are problematical."""

        ntot = len(self)
        bad = self.get_mask()
        nbad = len(bad[bad])
        flagged = np.bitwise_and(self.bmask, ANY_FLAG) > 0
        nflag = len(flagged[flagged])

        print(f"There were {nbad} bad data points out of {ntot}")
        print(f"There were {nflag} flagged data points out of {ntot}")

        if nflag:
            print('\nFlags raised:\n')
            for fname, flag in FLAGS:
                if flag != ALL_OK and flag != ANY_FLAG:
                    match = np.bitwise_and(self.bmask, flag) > 0
                    nfl = len(match[match])
                    if nfl:
                        print(f"   Flag = {fname} was raised for {nfl} points out of {ntot}")

        if nbad or nflag:
            print('\nPoint-by-point (bad and/or flagged):\n')

            for i in range(len(self.t)):
                if bad[i] or flagged[i]:

                    if flagged[i]:
                        flags_raised = []
                        for fname, flag in FLAGS:
                            if flag != ANY_FLAG: 
                                if self.bmask[i] & flag:
                                    flags_raised.append(fname)
                        flags_raised = ', '.join(flags_raised)
                    else:
                        flags_raised = "none"

                    print(f"   Index {i}: (t,y,ye) = ({self.t[i]},{self.y[i]},{self.ye[i]}), bad = {bad[i]}, flags raised = {flags_raised}")

    def to_mag(self, inplace=True):
        """
        Convert to magnitudes, i.e. take -2.5*log10(y) (inverse of from_mag)
        """
        if inplace:
            self.ye *= (2.5 / np.log(10.)) / self.y
            self.y = -2.5 * np.log10(self.y)
        else:
            new_ts = copy.deepcopy(self)
            new_ts.y = -2.5 * np.log10(self.y)
            new_ts.ye = (2.5 / np.log(10.)) * self.ye / self.y
            return new_ts

    def from_mag(self, inplace=True):
        """
        Convert from magnitudes to a linear scale (inverse of to_mag)
        """
        if inplace:
            self.y = 10**(-self.y/2.5)
            self.ye *= (np.log(10.) / 2.5) * self.y

        else:
            new_ts = copy.deepcopy(self)
            new_ts.y = 10**(-self.y/2.5)
            new_ts.ye = (np.log(10.) / 2.5) * new_ts.y * self.ye
            return new_ts

    def write(self, fname, lcurve=False, bitmask=None, **kwargs):
        """Writes out the Tseries to an ASCII file using numpy.savetxt. Other
        arguments to savetxt can be passed via kwargs. The data are
        written with a default choice of "fmt" sent to savetxt unless
        overridden using kwargs.  The default "fmt" uses a high
        precision on the times based upon typical usage for ULTRACAM
        derived (MJD) data.

        Arguments::

           fname : str
             output file name

           lcurve : bool
             if True, tack on two extra columns to make it easy to read the data
             into the program lcurve. This does not write bad data or the bitmask
             but uses the supplied bitmask to set weights on bad points to zero.
             If False all data including bad data are written and the fill bitmask
             value is given. Explanotory header info is tacked onto to whatever header
             is passed via kwargs to savetxt.

           bitmask : None | int
             if lcurve=True, this is used to set the weight on matching points to 0.
             This allows bad data to be plotted but not fitted e.g. by "lroche"

        """

        if lcurve:
            if self.te is None:
                raise ValueError("Can only write to an lcurve file with exposure times")

            # Avoid bad data in the case of lcurve data
            good = ~self.get_mask()
            mask = self.get_mask(bitmask)
            wgt = np.ones_like(good[good])
            wgt[mask[good]] = 0.0

            # we don't supply bitmask flags in this case
            header = kwargs.get("header", "") + \
                """
The data are written in a format suitable for using in the lcurve
light-curve modelling program and include two final columns of weights
and integer sub-division factors. The weights might be zero for data
thought to be dubious.

-----------------------------------------------------------------
Data written by hipercam.hlog.Tseries.write.

Six columns: times exposures fluxes flux-errors weights sub-div-facs
                """
            fmt = kwargs.get("fmt", "%.16g %.3e %.6e %.3e %.1f 1")
            data = np.column_stack([self.t[good], self.te[good], self.y[good], self.ye[good], wgt])

        else:

            flags = []
            for flag, message in FLAG_MESSAGES.items():
                flags.append(f"   Flag = {flag:6d}: {message}")
            flags = '\n'.join(flags)

            header = kwargs.get("header", "") + \
                f"""
Unrecoverably bad data are indicated by "nan" (not-a-number)
values. There are also a series of flags which are combined into a
single integer "bitmask". The values & meanings of these bitmask flags
are as follows:

{flags}

Thus for instance a value of 20 (=16+4) would indicate that the
annulus used to measure the sky had crossed the edge of the data
window and that the saturation level had been breached. A value of 0
means no flags were set. The flags vary in their severity and do not
always indicate that data is to be thrown away. There will also be
some data affected by glitches such as cosmic rays, so don't expect
all bad data to be indicated as such.

-----------------------------------------------------------------
Data written by hipercam.hlog.Tseries.write.

"""

            if self.te is None:
                header += "Four columns: times y-values y-errors bitmask"
                fmt = kwargs.get("fmt", "%.16g %.6e %.3e %d")
                data = np.column_stack([self.t, self.y, self.ye, self.bmask])
            else:
                header += "Five columns: times exposures y-values y-errors bitmask"
                fmt = kwargs.get("fmt", "%.16g %.3e %.6e %.3e %d")
                data = np.column_stack([self.t, self.te, self.y, self.ye, self.bmask])

        kwargs["header"] = header

        # Finally we are ready for savetxt
        np.savetxt(fname, data, fmt=fmt, **kwargs)

    @classmethod
    def read(cls, fname, has_exposures=True):
        """
        Creates a Tseries from an ASCII column data file. This
        should either contain times, y-values, y-errors, bitmask
        or times, exposures, y-values, y-errors, bitmask
        """
        if has_exposures:
            t, te, y, ye, bmask = np.loadtxt(fname, unpack=True)
            return Tseries(t, y, ye, mask, te)
        else:
            t, y, ye, bmask = np.loadtxt(fname, unpack=True)
            return Tseries(t, y, ye, mask, None)


def scatter(
        axes, xts, yts, color="b", fmt=".", bitmask=None,
        flagged=False, capsize=0, errx=True, erry=True,
        **kwargs):
    """
    Plots data in one Tseries versus another.

    Arguments::

         axes : Axes
            An Axes instance to plot to

         xts : Tseries
            The X-axis data

         yts : Tseries
            The Y-axis data

         color: str | rgb tuple
            matplotlib colour

         fmt: str
            matplotlib marker symbol to use

         bitmask : None | int
            bitmask to remove bad points. See 'get_mask' for usage.

         flagged : bool
            set True to plot points that match the bitmask rather than
            ones that don't.

         capsize : float
            if error bars are plotted with points, this sets
            the length of terminals

         errx : bool
            True / False for bars indicating exposure length (i.e. +/- 1/2
            whatever is in the te array)

         erry : bool
            True / False for vertical error bars or not

         kwargs : extra arguments
            These will be fed to the plot routine which is either
            matplotlib.pyplot.errorbar or matplotlib.pyplot.plot.
            e.g. 'ms=2' will set the markersize to 2.
    """

    # Generate boolean array of the points to plot
    if flagged:
        plot = xts.get_mask(bitmask) & yts.get_mask(bitmask)
    else:
        plot = ~(xts.get_mask(bitmask) & yts.get_mask(bitmask))

    # create arrays to plot
    x = xts.y[plot]
    y = yts.y[plot]

    errx &= np.any(x.ye[~np.isnan(x.ye)] > 0)
    erry &= np.any(y.ye[~np.isnan(y.ye)] > 0)

    kwargs["color"] = color

    if errx and erry:
        xe = xts.ye[plot]
        ye = yts.ye[plot]
        kwargs["fmt"] = fmt
        kwargs["capsize"] = capsize
        axes.errorbar(x, y, ye, xe, **kwargs)

    elif errx:
        xe = xts.ye[plot]
        kwargs["fmt"] = fmt
        kwargs["capsize"] = capsize
        axes.errorbar(x, y, xerr=te, **kwargs)

    elif erry:
        ye = self.ye[plot]
        kwargs["fmt"] = fmt
        kwargs["capsize"] = capsize
        axes.errorbar(x, y, ye, **kwargs)

    else:
        axes.plot(x, y, fmt, **kwargs)
