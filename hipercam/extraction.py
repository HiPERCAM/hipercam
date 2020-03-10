# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes and functions used to detect and characterise
objects in photometric data.

Adds optional dependency on sep package ()
"""

import struct
import sep
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from numpy.lib import recfunctions
import numpy as np

# to save space
SMALL_TYPE = np.dtype(
    [
        ('thresh', '<f4'), ('npix', '<i4'), ('tnpix', '<i4'), ('xmin', '<i4'), ('xmax', '<i4'),
        ('ymin', '<i4'), ('ymax', '<i4'), ('x', '<f8'), ('y', '<f8'), ('x2', '<f4'), ('y2', '<f4'),
        ('xy', '<f4'), ('errx2', '<f4'), ('erry2', '<f4'), ('errxy', '<f4'), ('a', '<f4'), ('b', '<f4'),
        ('theta', '<f4'), ('cxx', '<f4'), ('cyy', '<f4'), ('cxy', '<f4'), ('cflux', '<f4'),
        ('flux', '<f4'), ('cpeak', '<f4'), ('peak', '<f4'), ('xcpeak', '<i4'), ('ycpeak', '<i4'),
        ('xpeak', '<i4'), ('ypeak', '<i4'), ('flag', '<i4'), ('fwhm', '<f4'), ('hfd', '<f4')
    ]
)

def findStars(wind, thresh, kernel_fwhm, return_bkg=False):
    """
    Use sep to find objects in image. Not sure the outputs returned
    are correct for xbin != ybin

    Parameters
    ----------
    wind : `~hipercam.window.Window`
        window in which to detect objects.

    thresh : float
        threshold for object detection, in muliples of background RMS.


    kernel_fwhm : float
        Image is convolved with a Gaussian kernel of this FWHM prior to object
        detection. Should be set to something similar to the typical FWHM in image.

    return_bkg : bool
        True to return the background as calculated by sep.Background

    Returns
    -------
    objects : np.ndarray
        Extracted object parameters (structured array).
        For available fields, see sep documentation.
        http://sep.readthedocs.io/en/v1.0.x/api/sep.extract.html#sep.extract

        Most quantities are converted to unbinned coordinates and pixels,
        apart from npix, and tnpix, the number of pixels in the objects.

    bkg : `~sep.Background` [if return_bkg]
        The background estimated by 'sep'. Use sep.subfrom
        to subtract from data.

    """
    # ensure float type, c-contiguous and native byte order
    data = wind.data.astype('float')

    bkg = sep.Background(data)
    bkg.subfrom(data)  # in-place background subtraction

    sigma = 2.5 * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    objects = sep.extract(data, thresh, err=bkg.globalrms, clean=True,
                          filter_kernel=kernel.array)

    # add crude FWHM estimate and HFD measurement
    fwhm = 2. * np.sqrt(np.log(2.) * (objects['a']**2 + objects['b']**2))
    hfr, flags = sep.flux_radius(data, objects['x'], objects['y'], 25*np.ones_like(objects['x']),
                                 0.5, normflux=objects['cflux'])
    bad_hfd = np.logical_or(flags != 0, hfr >= 25)
    hfr[bad_hfd] = np.nan
    objects = recfunctions.append_fields(objects, ('fwhm', 'hfd'), (fwhm, 2*hfr))

    # convert to un-binned pixels
    objects['xmin'] = (wind.x(objects['xmin']-(wind.xbin-1)/2)).astype(np.int32)
    objects['xmax'] = (wind.x(objects['xmin']+(wind.xbin-1)/2)).astype(np.int32)
    objects['ymin'] = (wind.x(objects['ymin']-(wind.ybin-1)/2)).astype(np.int32)
    objects['ymax'] = (wind.x(objects['ymin']+(wind.ybin-1)/2)).astype(np.int32)
    for key in ('x', 'xcpeak', 'xpeak'):
        objects[key] = wind.x(objects[key])
    for key in ('y', 'ycpeak', 'ypeak'):
        objects[key] = wind.y(objects[key])
    for key in ('xy', 'cxy'):
        objects[key] *= wind.xbin * wind.ybin
    for key in ('x2', 'a', 'errx2', 'cxx'):
        objects[key] *= wind.xbin
    for key in ('y2', 'b', 'erry2', 'cyy'):
        objects[key] *= wind.ybin
    for key in ('fwhm', 'hfd'):
        objects[key] *= np.sqrt(wind.xbin * wind.ybin)

    # change the data type to save on space
    objects = objects.astype(SMALL_TYPE)

    if return_bkg:
        return (objects, bkg)
    else:
        return objects


