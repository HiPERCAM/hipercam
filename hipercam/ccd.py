# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Class to represent a CCD and a multi-CCD
"""

# Standard pre-amble from astropy
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.extern import six

import numpy as np
from astropy.io import fits
from .core import *
from .group import *
from .window import *

class CCD(Group):
    """
    Class representing a CCD as a Group of Windata objects
    plus a few FITS header items
    """
    def __init__(self, winds, nxtot, nytot, head=None):
        """
        Constructs a :class:`CCD`

        Arguments::

          winds : (dict)
              dictionary of Windat objects. The dictionary keys should be integers.

          nxtot : (int)
              Unbinned X-dimension of CCD

          nytot : (int)
              Unbinned Y-dimension of CCD

          head : (astropy.io.fits.Header)
              a header which will be written along with each of the sub-images
              if writing to a file.
        """
        super(CCD,self).__init__(winds)
        self.nxtot = nxtot
        self.nytot = nytot
        self.head = head

    def __repr__(self):
        return 'CCD(winds=' + super(CCD,self).__repr__() + \
                            ', nxtot=' + repr(self.nxtot) + \
                            ', nytot=' + repr(self.nytot) + \
                            ', head=' + repr(self.head) + ')'
