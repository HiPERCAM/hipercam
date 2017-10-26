"""
Cython routines
"""

import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport sqrt

ITYPE = np.float32
ctypedef np.float32_t ITYPE_t

OTYPE = np.float32
ctypedef np.float32_t OTYPE_t

# Next two decorators reduce time from 108ms to 75ms
# in the test I carried out.
@cython.cdivision(True)
@cython.boundscheck(False)
def avgstd(np.ndarray[ITYPE_t, ndim=3] cube, float sigma):
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
    assert (cube.dtype == ITYPE) and (sigma > 1.)

    # Output arrays: average, standard deviation and numbers of frames
    cdef np.ndarray[OTYPE_t, ndim=2] avg = np.empty((ny,nx), dtype=OTYPE)
    cdef np.ndarray[OTYPE_t, ndim=2] std = np.empty((ny,nx), dtype=OTYPE)
    cdef np.ndarray[int, ndim=2] nfm = np.empty((ny,nx), dtype=np.int32)

    # Temporary arrays for storing values and rejection flags
    cdef np.ndarray[ITYPE_t, ndim=1] vals = np.empty((nf,), dtype=ITYPE)
    cdef np.ndarray[int, ndim=1] ok = np.empty((nf,), dtype=np.int32)

    # Other temporaries
    cdef OTYPE_t tavg, tstd, thresh
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

