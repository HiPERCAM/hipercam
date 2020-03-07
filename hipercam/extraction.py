# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes and functions used to detect and characterise
objects in photometric data.

Adds optional dependency on sep package ()
"""

import sep
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from numpy.lib import recfunctions
import numpy as np


def findStars(wind, thresh, kernel_fwhm, return_bkg=False):
    """
    Use sep to find objects in image.

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

        Quantities are in un-binned pixels

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
    for key in ('x', 'xmin', 'xmax', 'xcpeak', 'xpeak'):
        objects[key] = wind.x(objects[key])
    for key in ('y', 'ymin', 'ymax', 'ycpeak', 'ypeak'):
        objects[key] = wind.y(objects[key])
    for key in ('npix', 'tnpix', 'xy', 'cxy'):
        objects[key] *= wind.xbin * wind.ybin
    for key in ('x2', 'a', 'errx2', 'cxx'):
        objects[key] *= wind.xbin
    for key in ('y2', 'b', 'erry2', 'cyy'):
        objects[key] *= wind.ybin
    for key in ('fwhm', 'hfd'):
        objects[key] *= np.sqrt(wind.xbin * wind.ybin)

    if return_bkg:
        return (objects, bkg)
    else:
        return objects
