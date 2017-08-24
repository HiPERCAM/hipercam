# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core classes and functions for the hipercam package
"""

from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft

__all__ = (
    'FIELD', 'HCAM', 'LIST', 'APER', 'HRAW', 'add_extension',
    'HipercamError', 'HipercamWarning', 'HCM_NXTOT',
    'HCM_NYTOT', 'detect'
)

# Constants for general use

# Standard file extensions
FIELD = '.fld'
HCAM = '.hcm'
LIST = '.lis'
APER = '.ape'
HRAW = '.fits'

# Maximum dimensions of HiPERCAM CCDs
HCM_NXTOT, HCM_NYTOT = 2048, 1040

def add_extension(fname, ext):
    """Add extension ext to a file name if it is not already there, and returns
    the revised name

    """
    if len(ext) and not fname.endswith(ext):
        return fname + ext
    else:
        return fname

def detect(img, fwhm, fft=True):
    """
    Detects the position of one star in an image. Works by convolving the
    image with a gaussian of FWHM = fwhm, and returning the location of the
    brightest pixel of the result. The convolution reduces the chance of
    problems being caused by cosmic rays, although if there is more overall
    flux in a cosmic ray than the star, it will go wrong, so some form of
    prior cleaning is advisable in such cases. It is up to the user to pass a
    small enough image so that the star of interest is the brightest object.

    This routine is intended to provide a first cut in position for more
    precise methods to polish.

    Arguments::

       img   : (2D numpy array)
          the image

       fwhm  : (float)
          Gaussian FWHM in pixels. If <= 0, there will be no convolution, although
          this is not advisable as a useful strategy.

       fft   : (bool)
          The astropy.convolution routines are used. By default FFT-based
          convolution is applied as it scales better with fwhm, especially for
          fwhm >> 1, however the direct method (fft=False) may be faster for
          small fwhm values.

    Returns: (x,y), the location of brightest pixel in the convolved image,
    with the lower-left pixel taken to lie at (0,0).
    """
    if fwhm > 0:
        kern = Gaussian2DKernel(fwhm/np.sqrt(8*np.log(2)))
        if fft:
            cimg = convolve_fft(img, kern)
        else:
            cimg = convolve(img, kern)
    else:
        cimg = img

    # locate coords of maximum pixel.
    y, x = np.unravel_index(cimg.argmax(), cimg.shape)

    return(x,y)

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

