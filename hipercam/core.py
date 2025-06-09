# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core data, functions and classes for the hipercam package
"""

from collections import OrderedDict as Odict
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
from importlib.metadata import version as _version, PackageNotFoundError

__all__ = (
    "FIELD",
    "HCAM",
    "LIST",
    "APER",
    "HRAW",
    "RED",
    "LOG",
    "DFCT",
    "SEP",
    "TBTS",
    "HipercamError",
    "HipercamWarning",
    "DMINS",
    "CNAMS",
    "CIS",
    "version",
    "gregorian_to_mjd",
    "mjd_to_gregorian",
    "fday_to_hms",
    "ALL_OK",
    "NO_FWHM",
    "NO_SKY",
    "SKY_AT_EDGE",
    "TARGET_AT_EDGE",
    "TARGET_NONLINEAR",
    "TARGET_SATURATED",
    "REDUCE_FILE_VERSION",
    "NO_EXTRACTION",
    "NO_DATA",
    "CLOUDS",
    "JUNK",
    "ANY_FLAG",
    "BAD_FLAT",
    "BAD_COLUMN",
    "BAD_TIME",
    "FLAGS",
    "OUTLIER",
    "FLAG_MESSAGES",
    "FRNG",
)

# Constants for general use

# Version of the reduce file in operation (used by 'reduce' and 'genred')
# Format: YYYYMMDD(.#) where the optional .# part is an integer to allow for
# multiple versions in a day, although that should be rare I hope.
REDUCE_FILE_VERSION = "20230222"

# Standard file extensions
FIELD = ".fld"
HCAM = ".hcm"
LIST = ".lis"
APER = ".ape"
HRAW = ".fits"
RED = ".red"
LOG = ".log"
DFCT = ".dft"
SEP = ".sep"
TBTS = ".tbts"
FRNG = ".frng"

# number of minutes in a day
DMINS = 1440.0

# For compatibility between PGPLOT and matplotlib, establish a uniform
# set of colours. These are fed directly to matplotlib, and used to override
# the colour indices used by PGPLOT.
CIS = (
    (1, 1, 1),  # 0 -- white
    (0, 0, 0),  # 1 -- black
    (0.898, 0, 0),  # 2 -- red
    (0.082, 0.690, 0.102),  # 3 -- green
    (0.012, 0.263, 0.875),  # 4 -- blue
    (1, 0.008, 0.553),  # 5 -- pink
    (0.494, 0.118, 0.612),  # 6 -- purple
    (0.024, 0.60, 0.675),  # 7 -- turquoise
    (0.976, 0.451, 0.024),  # 8 -- orange
    (0.808, 0.702, 0.004),  # 9 -- mustard
    (0.761, 0.0, 0.471),  # 10 -- magenta
    (0.682, 0.443, 0.505),  # 11 -- mauve
    (0.635, 0.643, 0.082),  # 12 -- vomit
    (0.7, 0.008, 0.0),  # 13 -- darkred
    (0.498, 0.372, 0),  # 14 -- mud
    (0.95, 0.95, 0),  # 15 -- yellow
)

CNAMS = {
    "white": 0,
    "black": 1,
    "red": 2,
    "green": 3,
    "blue": 4,
    "pink": 5,
    "purple": 6,
    "turquoise": 7,
    "orange": 8,
    "mustard": 9,
    "magenta": 10,
    "mauve": 11,
    "vomit": 12,
    "darkred": 13,
    "mud": 14,
    "yellow": 15,
}

# Bit masks (used in reduce.py and hlog.py)
bmask = np.uint32
ALL_OK = bmask(0)  # No problems detected (no bits set)
ANY_FLAG = ~ALL_OK  # Matches any point flagged in some way

# Now the flags; maximum number possible = 32
NO_FWHM = bmask(1 << 0)  # No FWHM, even though variable apertures are being used
NO_SKY = bmask(1 << 1)  # No sky pixels at all
SKY_AT_EDGE = bmask(1 << 2)  # Sky aperture edge of window
TARGET_AT_EDGE = bmask(1 << 3)  # Target aperture overlaps edge of window
TARGET_SATURATED = bmask(1 << 4)  # At least one pixel in target above saturation level
TARGET_NONLINEAR = bmask(1 << 5)  # At least one pixel in target above non-linear level
NO_EXTRACTION = bmask(1 << 6)  # No extraction possible
NO_DATA = bmask(1 << 7)  # No valid pixels in aperture
CLOUDS = bmask(1 << 8)  # Point affected by clouds
JUNK = bmask(1 << 9)  # Unspecified junk data (e.g. cosmic ray hit)
BAD_FLAT = bmask(1 << 10)  # Bad flat field pixel (e.g. deep dust speck)
BAD_COLUMN = bmask(1 << 11)  # Bad column (e.g. target aperture includes such data)
BAD_TIME = bmask(1 << 12)  # Bad time
OUTLIER = bmask(1 << 13)  # Point far from local trend (for Tseries)

# all flags for useful reference in other places
FLAGS = (
    ("ALL_OK", ALL_OK),
    ("ANY_FLAG", ANY_FLAG),
    ("NO_FWHM", NO_FWHM),
    ("NO_SKY", NO_SKY),
    ("SKY_AT_EDGE", SKY_AT_EDGE),
    ("TARGET_AT_EDGE", TARGET_AT_EDGE),
    ("TARGET_SATURATED", TARGET_SATURATED),
    ("TARGET_NONLINEAR", TARGET_NONLINEAR),
    ("NO_EXTRACTION", NO_EXTRACTION),
    ("NO_DATA", NO_DATA),
    ("CLOUDS", CLOUDS),
    ("JUNK", JUNK),
    ("BAD_FLAT", BAD_FLAT),
    ("BAD_COLUMN", BAD_COLUMN),
    ("BAD_TIME", BAD_TIME),
    ("OUTLIER", OUTLIER),
)

# messages if various bitflags are set
FLAG_MESSAGES = Odict(
    (
        (NO_FWHM, "no FWHM could be measured"),
        (NO_SKY, "zero sky pixels"),
        (SKY_AT_EDGE, "sky aperture overlapped the edge of the data window"),
        (TARGET_AT_EDGE, "target aperture overlapped the edge of the data window"),
        (TARGET_SATURATED, "target aperture had saturated pixels"),
        (TARGET_NONLINEAR, "target aperture had non-linear pixels"),
        (NO_EXTRACTION, "no extraction was possible"),
        (NO_DATA, "there were no valid pixels in the target aperture"),
        (CLOUDS, "marked as affected by clouds"),
        (JUNK, "junk data of unspecified nature"),
        (BAD_FLAT, "bad flat field feature in target aperture"),
        (BAD_COLUMN, "bad column in in target aperture"),
        (BAD_TIME, "the time was flagged as bad"),
        (OUTLIER, "identified as an outlier"),
    )
)


def version():
    """Returns version number of installed HiPERCAM pipeline"""
    try:
        return _version("hipercam")
    except PackageNotFoundError:
        return "not found"


def gregorian_to_mjd(year, month, day):
    """This returns the Modified Julian Day number corresponding to the
    start of the Gregorian date supplied. Algorithm taken from
    Duffet-Smith. Tested against astropy.time.Time. Designed to be
    simple and fast.

    Parameters
    ----------
    year : int
        Year (valid range: 1582 onwards)
    month : int
        Month of year, 1 to 12 inclusive.
    day : int
        day of month (1 to 31; will not check month by month)

    Returns
    -------
    mjd : int
        Modified Julian Day number
    """

    if year < 1582 or month < 1 or month > 12 or day < 1 or day > 31:
        raise ValueError(
            "Supplied date ({:d}-{:d}-{:d}) is not valid".format(year, month, day)
        )

    if month < 3:
        month += 12
        year -= 1
    A = year // 100
    B = 2 - A + A // 4
    C = int(365.25 * year)
    D = int(30.6001 * (month + 1))

    return B + C + D - 679006 + day


def mjd_to_gregorian(mjd):
    """This returns the year, month, and day corresponding to the supplied
    MJD. From Duffet-Smith. Checked against 'gregorian_to_mjd'.

    Parameters
    ----------
    mjd : int
       Modified Julian Day number

    Returns
    -------
    Returns a tuple of (year,month,day)

    """

    jd = mjd + 2400001
    if jd > 2299160:
        A = int((jd - 1867216.25) / 36524.25)
    else:
        A = jd
    B = jd + 1 + A - A // 4
    C = B + 1524
    D = int((C - 122.1) / 365.25)
    E = int(365.25 * D)
    G = int((C - E) / 30.6001)
    day = C - E - int(30.6001 * G)
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    return (year, month, day)


def fday_to_hms(fday):
    """Returns a tuple of (hours, minutes, seconds) given a fraction of a
    day.  The hours and minutes are integers, the seconds are floating
    point.

    fday should lie in the interval [0,1)
    """
    hours = int(24 * fday)
    minutes = int(60 * (24 * fday - hours))
    seconds = 86400 * fday - 3600 * hours - 60 * minutes
    return (hours, minutes, seconds)


class HipercamError(Exception):
    """
    Class for the hipercam package errors
    """


class HipercamWarning(AstropyUserWarning):
    """
    Class for hipercam package warnings. Use with warnings.warn
    """


# if __name__ == '__main__':
#    import doctest
#    doctest.testmod()
