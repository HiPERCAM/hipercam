# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core classes and functions for the hipercam package
"""

from astropy.utils.exceptions import AstropyUserWarning
import numpy as np

__all__ = (
    'FIELD', 'HCAM', 'LIST', 'APER', 'HRAW', 'RED',
    'add_extension', 'HipercamError', 'HipercamWarning', 'CIS',
    'DMINS'
)

# Constants for general use

# Standard file extensions
FIELD = '.fld'
HCAM = '.hcm'
LIST = '.lis'
APER = '.ape'
HRAW = '.fits'
RED  = '.red'

# number of minutes in a day
DMINS = 1440.

# For compatibility between PGPLOT and matplotlib, establish a uniform
# set of colours. These are fed directly to matplotlib, and used to override
# the colour indices used by PGPLOT. Could extend from 0 to 15 max.
CIS = (
    (1,1,1),      # 0 -- white
    (0,0,0),      # 1 -- black
    (0.5,0,0),    # 2 -- red
    (0,0.4,0),    # 3 -- green
    (0,0,0.5),    # 4 -- blue
    (0.5,0.3,0),  # 5 -- pinkish
    (0.5,0,0.5),  # 6 -- purple
    (0,0.5,0.5),  # 7 -- turquoise
    )


def add_extension(fname, ext):
    """Add extension ext to a file name if it is not already there, and returns
    the revised name

    """
    if len(ext) and not fname.endswith(ext):
        return '{}{}'.format(fname, ext)
    else:
        return fname

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

