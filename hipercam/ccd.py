# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Class to represent a CCD and a multi-CCD
"""

# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

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

    def min(self):
        """
        Returns the minimum value of the :class:`CCD`.
        """
        winds = self.values()
        vmin = winds[0].min()
        for wind in winds[1:]:
            vmin = min(vmin, wind.min())
        return vmin

    def max(self):
        """
        Returns the maximum value of the :class:`CCD`.
        """
        winds = self.values()
        vmax = winds[0].max()
        for wind in winds[1:]:
            vmax = min(vmax, wind.max())
        return vmax

    def mean(self):
        """
        Returns the mean value of the :class:`Windata`.
        """
        npix, total = 0, 0.
        for wind in self.values():
            npix += wind.size
            total += wind.sum()
        return total / float(npix)

    def percentile(self, q):
        """
        Computes percentile(s) of the :class:`CCD`.

        Arguments::

        q : (float or sequence of floats)
          Percentile(s) to use, in range [0,100]
        """

        # Flatten into a single 1D array
        arrs = []
        for wind in ccd.values():
            arrs.append(wind.data.flatten())
        arr = np.concatenate(arrs)

        # Then compute percentiles
        return np.percentile(arr, q)

    def __repr__(self):
        return 'CCD(winds=' + super(CCD,self).__repr__() + \
                            ', nxtot=' + repr(self.nxtot) + \
                            ', nytot=' + repr(self.nytot) + \
                            ', head=' + repr(self.head) + ')'
