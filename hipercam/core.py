# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core classes and functions for the hipercam package
"""

# Standard pre-amble from astropy
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

# Constants for general use

# Standard file extensions
FIELD = '.fld'
HCAM = '.hcm'

class HipercamError (Exception):
    """
    Error class for the hipercam package
    """

#if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

