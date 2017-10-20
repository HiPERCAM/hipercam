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
  >> hlog = hcam.hlog.Hlog.fromLog('ts.log')
  >> targ = hlog.tseries('2','t')
  >> comp = hlog.tseries('2','c')
  >>
  >> ratio = targ / comp
  >> ratio.mplot(plt, 'r')
  >> plt.show()


"""

import numpy as np

from .core import *
from . import utils

__all__ = ('Hlog', 'Tseries')

# Bit masks
NO_FWHM         = 0b1     # even though variable apertures are being used
NO_SKY          = 0b10    # no sky pixels
SKY_OFF_EDGE    = 0b100   # sky aperture off edge of window
TARGET_OFF_EDGE = 0b1000  # target aperture off edge of window
ANY = NO_FWHM | NO_SKY | SKY_OFF_EDGE | TARGET_OFF_EDGE

class Hlog(dict):
    """
    Class to represent a HiPERCAM log as produced by reduce.
    Based on dictionaries, Hlog files contain numpy structured
    arrays for each CCD which can be accessed by the CCD label
    name.
    """

    @classmethod
    def fromLog(cls, fname):
        """
        Loads a HiPERCAM log file written by reduce. Each CCD is loaded
        into a separate structured array, returned in a dictionary labelled
        by the CCD key. 
        """
        hlog = cls()
        read_cnames = False
        read_dtypes = False
        read_data = False
        cnames = {}
        dtype_defs = {}
        dtypes = {}
        with open(fname) as fin:
            for line in fin:

                if read_cnames:
                    # reading the column names
                    if line.find('End of column name definitions') > -1:
                        read_cnames = False
                    elif line.find('=') > -1:
                        # store the column names
                        cnam = line[1:line.find('=')].strip()
                        cnames[cnam] = line[line.find('=')+1:].strip().split()[1:]

                elif read_dtypes:
                    # reading the data types
                    if line.find('End of data type definitions') > -1:
                        read_dtypes = False
                        # can now create the record array dtype
                        for cnam in cnames:
                            dtypes[cnam] = np.dtype(
                                list(zip(cnames[cnam], dtype_defs[cnam]))
                            )
                        read_data = True

                    elif line.find('=') > -1:
                        # store the data types
                        cnam = line[1:line.find('=')].strip()
                        dtype_defs[cnam] = line[line.find('=')+1:].strip().split()[1:]

                        # get ready for the data
                        hlog[cnam] = []

                elif line.find('Start of column name definitions') > -1:
                    read_cnames = True

                elif line.find('Start of data type definitions') > -1:
                    read_dtypes = True

                elif read_data and not line.startswith('#'):
                    # read a data line
                    items = line.strip().split()
                    cnam = items.pop(0)
                    dts = dtype_defs[cnam]

                    # convert types
                    for n in range(len(items)):
                        dt = dts[n]
                        if dt.startswith('f'):
                            items[n] = float(items[n])
                        elif dt.startswith('i') or dt.startswith('u'):
                            items[n] = int(items[n])
                        elif dt == '?':
                            items[n] = bool(items[n])

                    hlog[cnam].append(tuple(items))

        for cnam in hlog:
            hlog[cnam] = np.array(hlog[cnam],dtypes[cnam])

        return hlog

    def tseries(self, cnam, apnam, name='counts'):
        """
        Returns with a Tseries corresponding to CCD cnam and
        aperture apnam. By default it accesses the 'counts',
        but 'x', 'y', 'fwhm', 'beta' and 'sky' are alternative
        choices.
        """
        ccd = self[cnam]
        return Tseries(
            ccd['MJD'], ccd['{:s}_{:s}'.format(name,apnam)],
            ccd['{:s}e_{:s}'.format(name,apnam)],
            ccd['flag_{:s}'.format(apnam)]
        )

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

    """

    def __init__(self, t, y, ye, mask):
        self.t = t
        self.y = y
        self.ye = ye
        self.mask = mask

    def mplot(self, axes, colour='b', fmt='.', mask=ANY, capsize=0,
              **kwargs):
        """
        Plots a Tseries to a matplotlib Axes instance, only
        plotting points without any of the same bits set as
        'mask'. Points with negative errors are plotted
        without error bars.
        """
        ok = np.bitwise_and(self.mask, mask) == 0

        err = ok & (self.ye > 0.)
        axes.errorbar(
            self.t[err],self.y[err],self.ye[err],
            fmt=fmt, color=colour, capsize=capsize, **kwargs)

        nerr = ok & (self.ye <= 0.)
        axes.plot(self.t[nerr], self.y[nerr], fmt, color=colour, **kwargs)

    def __repr__(self):
        return 'Tseries(t={!r}, y={!r}, ye={!r}, mask={!r}'.format(
            self.t, self.y, self.ye, self.mask)

    def __truediv__(self, other):
        """Divides the Tseries by 'other' returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if not all(self.t == other.t):
                raise ValueError('input times do not match')

            y = self.y / other.y
            ye = np.empty_like(y)
            ok = (self.ye > 0) & (other.ye > 0)
            ye[ok] = np.sqrt(
                self.ye[ok]**2 +
                (other.ye[ok]*self.y[ok]/other.y[ok])**2)/other.y[ok]
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # division by a constant or an array of constants
            y = self.y / other
            ye = np.empty_like(y)
            ok = (self.ye > 0)
            ye[ok] = self.ye[ok] / other
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t, y, ye, mask)

    def __mul__(self, other):
        """Multiplies the Tseries by 'other' returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if not all(self.t == other.t):
                raise ValueError('input times do not match')

            y = self.y * other.y
            ye = np.empty_like(y)
            ok = (self.ye > 0) & (other.ye > 0)
            ye[ok] = np.sqrt(
                (other.y[ok]*self.ye[ok])**2 +
                (self.ye[ok]*other.y[ok])**2
            )
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # multiplication by a constant or an array of constants
            y = self.y * other
            ye = np.empty_like(y)
            ok = (self.ye > 0)
            ye[ok] = self.ye[ok] * other
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t, y, ye, mask)

    def __add__(self, other):
        """Add 'other' to the Tseries returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if not all(self.t == other.t):
                raise ValueError('input times do not match')

            y = self.y + other.y
            ye = np.empty_like(y)
            ok = (self.ye > 0) & (other.ye > 0)
            ye[ok] = np.sqrt(self.ye[ok]**2 + other.y[ok]**2)
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # addition of a constant or an array of constants
            y = self.y + other
            ye = np.empty_like(y)
            ok = (self.ye > 0)
            ye[ok] = self.ye[ok]
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t, y, ye, mask)

    def __sub__(self, other):
        """Subtracts 'other' from the Tseries returning the result as another
        Tseries. Propagates errors where possible. 'other' can be a constant
        or another Tseries or a compatible numpy array. If it is a Tseries,
        the bit masks are bitwise_or-ed together. If any input errors are
        negative, the equivalent output errors are set = -1.
        """
        if isinstance(other, Tseries):
            if not all(self.t == other.t):
                raise ValueError('input times do not match')

            y = self.y - other.y
            ye = np.empty_like(y)
            ok = (self.ye > 0) & (other.ye > 0)
            ye[ok] = np.sqrt(self.ye[ok]**2 + other.y[ok]**2)
            ye[~ok] = -1
            mask = np.bitwise_or(self.mask, other.mask)

        else:
            # subtraction of a constant or an array of constants
            y = self.y - other
            ye = np.empty_like(y)
            ok = (self.ye > 0)
            ye[ok] = self.ye[ok]
            ye[~ok] = -1
            mask = self.mask.copy()

        return Tseries(self.t, y, ye, mask)

