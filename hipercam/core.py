# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core classes and functions for the hipercam package
"""

from astropy.utils.exceptions import AstropyUserWarning

__all__ = ('FIELD', 'HCAM', 'add_extension',
           'HipercamError', 'HipercamWarning')

# Constants for general use

# Standard file extensions
FIELD = '.fld'
HCAM = '.fits'

def add_extension(fname, ext):
    """Add extension ext to a file name if it is not already there, and returns
    the revised name

    """
    if len(ext) and not fname.endswith(ext):
        return fname + ext
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

