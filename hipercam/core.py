# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core data and a few classes for the hipercam package
"""

from astropy.utils.exceptions import AstropyUserWarning
import numpy as np

__all__ = (
    'FIELD', 'HCAM', 'LIST', 'APER', 'HRAW', 'RED',
    'HipercamError', 'HipercamWarning', 'DMINS',
    'LOG', 'CNAMS', 'CIS', 'ALL_OK', 'NO_FWHM',
    'NO_SKY', 'SKY_OFF_EDGE', 'TARGET_OFF_EDGE',
    'EXTRA_OFF_EDGE', 'DATA_SATURATED', 'ANY'
)

# Constants for general use

# Standard file extensions
FIELD = '.fld'
HCAM = '.hcm'
LIST = '.lis'
APER = '.ape'
HRAW = '.fits'
RED  = '.red'
LOG  = '.log'

# number of minutes in a day
DMINS = 1440.

# For compatibility between PGPLOT and matplotlib, establish a uniform
# set of colours. These are fed directly to matplotlib, and used to override
# the colour indices used by PGPLOT. Could extend from 0 to 15 max.
CIS = (
    (1,1,1),              # 0 -- white
    (0,0,0),              # 1 -- black
    (0.898,0,0),          # 2 -- red
    (0.082,0.690,0.102),  # 3 -- green
    (0.012,0.263,0.875),  # 4 -- blue
    (1,0.008,0.553),      # 5 -- pink
    (0.494,0.118,0.612),  # 6 -- purple
    (0.024,0.60,0.675),   # 7 -- turquoise
    (0.976,0.451,0.024),  # 8 -- orange
    (0.808,0.702,0.004),  # 9 -- mustard
    (0.761,0.,0.471),     # 10 -- magenta
    (0.682,0.443,0.505),  # 11 -- mauve
    (0.635,0.643,0.082),  # 12 -- vomit
    (0.7,0.008,0.0),      # 13 -- darkred
    (0.498, 0.372, 0),    # 14 -- mud 
    (0.8, 0.8, 0),        # 15 -- yellow 
)

CNAMS = {
    'white'  : 0,
    'black'  : 1,
    'red'    : 2,
    'green'  : 3,
    'blue'   : 4,
    'pink'   : 5,
    'purple' : 6,
    'turquoise' : 7,
    'orange' : 8,
    'mustard' : 9,
    'magenta' : 10,
    'mauve' : 11,
    'vomit' : 12,
    'darkred' : 13,
    'mud' : 14,
}

# Bit masks (used in reduce.py and hlog.py)
ALL_OK          = 0       # no bits set
NO_FWHM         = 1 << 0  # even though variable apertures are being used
NO_SKY          = 1 << 1  # no sky pixels
SKY_OFF_EDGE    = 1 << 2  # sky aperture off edge of window
TARGET_OFF_EDGE = 1 << 3  # target aperture off edge of window
EXTRA_OFF_EDGE  = 1 << 4  # extra target aperture off edge of window
DATA_SATURATED  = 1 << 5  # At least one pixel in target aperture

ANY = NO_FWHM | NO_SKY | SKY_OFF_EDGE | TARGET_OFF_EDGE | \
      EXTRA_OFF_EDGE | DATA_SATURATED

class HipercamError (Exception):
    """
    Class for the hipercam package errors
    """

class HipercamWarning (AstropyUserWarning):
    """
    Class for hipercam package warnings. Use with warnings.warn
    """

#if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

