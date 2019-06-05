"""
Cython routines
"""

import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport sqrt

FTYPE = np.float32
ctypedef np.float32_t FTYPE_t

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# Next two decorators reduce time from 108ms to 75ms
# in the test I carried out.
@cython.cdivision(True)
@cython.boundscheck(False)
def avgstd(np.ndarray[FTYPE_t, ndim=3] cube, float sigma):
    """
    Given a 3D cube of dimensions (nf,ny,nx) which can be thought of as 'nf'
    frames each of dimension (ny,nx), avgstd computes the mean, standard
    deviation and number of frames used to calculate them for each pixel after
    rejection at 'sigma' standard deviations across the first dimension. It
    returns three 2D frames each of dimension (ny,nx). The rejection is
    carried out in a cycle until no more rejection occurs.

    This routine is heavily cythonised and runs much faster than a pure-python
    implementation, but that comes with a restriction on the data type of
    'cube' as listed below.

    Arguments::

      cube : (3D numpy array, dtype=np.float32)
           the set of frames to process.

      sigma : (float)
           the rejection threshold. 3 to 4 goodish value. Must be > 1.

    Returns (avg,std,num) the average, standard deviation and number
    of contributing frames.
    """

    # Dimensions and indices
    cdef unsigned int nf = cube.shape[0]
    cdef unsigned int ny = cube.shape[1]
    cdef unsigned int nx = cube.shape[2]
    cdef unsigned int ix, iy, iz

    # Sanity checks
    assert (cube.dtype == FTYPE) and (sigma > 1.)

    # Output arrays: average, standard deviation and numbers of frames
    cdef np.ndarray[FTYPE_t, ndim=2] avg = np.empty((ny,nx), dtype=FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=2] std = np.empty((ny,nx), dtype=FTYPE)
    cdef np.ndarray[int, ndim=2] nfm = np.empty((ny,nx), dtype=np.int32)

    # Temporary arrays for storing values and rejection flags
    cdef np.ndarray[FTYPE_t, ndim=1] vals = np.empty((nf,), dtype=FTYPE)
    cdef np.ndarray[int, ndim=1] ok = np.empty((nf,), dtype=np.int32)

    # Other temporaries
    cdef FTYPE_t tavg, tstd, thresh
    cdef unsigned int num, numt
    cdef double sum

    # The main loop starts. Anything inside this is done nx*ny
    # times and better be fast...
    for iy in xrange(ny):
        for ix in xrange(nx):

            # extract values of given pixel for all frames, set flags to 1
            for iz in xrange(nf):
                vals[iz] = cube[iz,iy,ix]
                ok[iz] = 1

            # calculate initial values of tavg, tstd and num
            num = nf
            sum = 0.
            for iz in xrange(nf):
                sum += vals[iz]
            tavg = sum / num
            sum = 0.
            for iz in xrange(nf):
                sum += (vals[iz]-tavg)**2
            if num > 1:
                tstd = sqrt(sum/(num-1))
            else:
                tstd = 0.
                break

            # Now the rejection loop. Keep rejecting if
            # number of rejections is positive
            while True:
                # pre-compute rejection threshold
                thresh = sigma*tstd
                sum = 0.
                numt = 0
                for iz in xrange(nf):
                    if abs(vals[iz]-tavg) > thresh:
                        ok[iz] = 0

                    if ok[iz]:
                        sum += vals[iz]
                        numt += 1

                tavg = sum / numt
                sum = 0.
                for iz in xrange(nf):
                    if ok[iz]:
                        sum += (vals[iz]-tavg)**2
                if numt > 1:
                    tstd = sqrt(sum/(numt-1))
                else:
                    tstd = 0.
                    break

                if numt == num:
                  break
                num = numt

            # store the final numbers
            avg[iy,ix] = tavg
            std[iy,ix] = tstd
            nfm[iy,ix] = num

    # return the frames
    return (avg,std,nfm)

@cython.cdivision(True)
@cython.boundscheck(False)
def gaussian(
        np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] y,
        double sky, double height, double xcen, double ycen, double fwhm,
        int xbin, int ybin, int ndiv):
    """Returns a numpy array corresponding to the ordinate grids in xy set to a
    symmetric 2D Gaussian plus a constant. The profile is essentially defined
    by sky + height*exp(-alpha*r**2) where r is the distance from the centre
    and alpha is set to give the desired FWHM, but account is taken of the
    finite size of the pixels by summing over multiple points in each one.

    Arguments::

      x  : 2D numpy array
         the X ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      y  : 2D numpy array
         the Y ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      sky : float
         sky background

      height : float
         height of central peak

      xcen : float
         X-ordinate of centre

      ycen : float
         Y-ordinate of centre

      fwhm : float
         FWHM of the profile (unbinned pixels)

      xbin : int
         X-size of pixels in terms of unbinned pixels, i.e. the binning factor in X

      ybin : int
         Y-size of pixels in terms of unbinned pixels, i.e. the binning factor in Y

      ndiv : int
         Parameter controlling treatment of sub-pixellation. If > 0, every
         pixel will be sub-divided first into unbinned pixels, and then each
         unbinned pixels will be split into a square array of ndiv by ndiv
         points. The profile will be evaluated and averaged over each of these
         points. Thus if the pixels are binned xbin by ybin, there will be
         xbin*ybin*ndiv**2 evaluations per pixel. This is to cope with the
         case where the seeing is becoming small compared to the pixels, but
         obviously will slow things. To simply evaluate the profile once at
         the centre of each pixel, set ndiv = 0.

    Returns:: 2D numpy array containing the Gaussian plus constant evaluated
    on the ordinate grids in xy.

    """

    # Dimensions and indices
    cdef unsigned int ny = x.shape[0]
    cdef unsigned int nx = x.shape[1]
    cdef double alpha = 4.*np.log(2.)/fwhm**2
    cdef np.ndarray[DTYPE_t, ndim=2] rsq = np.empty((ny,nx), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] prof = np.zeros((ny,nx), dtype=DTYPE)
    cdef unsigned int ix, iy, isx, isy
    cdef double xoff, yoff, xsoff, ysoff, soff

    if ndiv > 0:
        # Complicated case with sub-pixellation allowed for
        soff = (ndiv-1)/(2*ndiv)
        for iy in range(ybin):
            # loop over unbinned pixels in Y
            yoff = iy-(ybin-1)/2 - soff
            for ix in range(xbin):
                # loop over unbinned pixels in X
                xoff = ix-(xbin-1)/2 - soff
                for isy in range(ndiv):
                    # loop over sub-pixels in y
                    ysoff = yoff + isy/ndiv
                    for isx in range(ndiv):
                        # loop over sub-pixels in x
                        xsoff = xoff + isx/ndiv
                        rsq = (x+xsoff-xcen)**2+(y+ysoff-ycen)**2
                        prof += np.exp(-alpha*rsq)

        return sky+(height/xbin/ybin/ndiv**2)*prof

    else:
        # Fast as possible, compute profile at pixel centres only
        rsq = (x-xcen)**2+(y-ycen)**2
        return sky+height*np.exp(-alpha*rsq)
