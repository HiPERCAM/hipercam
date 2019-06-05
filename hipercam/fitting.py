# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Code for profile fitting. Currently supports symmetric 2D Gaussian and
Moffat profiles plus constants. 
"""

from numba import jit
import numpy as np
from scipy.optimize import leastsq
from .core import *
from .window import *
from . import support

__all__ = (
    'combFit', 'fitMoffat', 'fitGaussian', 'moffat', 'gaussian'
)

def combFit(wind, method, sky, height, x, y, fwhm, fwhm_min, fwhm_fix,
            beta, read, gain, thresh, ndiv=0):
    """Fits a stellar profile in a :class:Window using either a 2D Gaussian
    or Moffat profile. This is a convenient wrapper of fitMoffat and
    fitGaussian because one often wants both options.

    Arguments::

        wind : :class:`Window`
            the Window containing the stellar profile to fit.

        method    : string
            fitting method 'g' for Gaussian, 'm' for Moffat

        sky       : float
            initial sky level

        height    : float
            initial peak height

        x         : float
            initial central X value

        y         : float
            initial central Y value

        fwhm      : float
            initial FWHM, unbinned pixels.

        fwhm_min  : float
            minimum FWHM, unbinned pixels.

        fwhm_fix  : float
            fix the FWHM (i.e. don't fit it)

        beta      : float [if method == 'm']
            exponent of Moffat function

        read      : float | array
            readout noise, RMS ADU

        gain      : float | array
            gain, electrons per ADU

        thresh    : float
            threshold in terms of RMS for rejection. The RMS is obtained from
            read and gain but scaled by the sqrt(chi**2/d.o.f). Typical value
            = 4

        ndiv      : int
            Parameter controlling treatment of sub-pixellation. If ndiv > 0,
            every pixel in `wind` will be sub-divided first into unbinned
            pixels, and then each unbinned pixel will be split into a square
            array of ndiv by ndiv points. The profile will be evaluated and
            averaged over each of these points. Thus if the pixels in `wind`
            are binned xbin by ybin, there will be xbin*ybin*ndiv**2
            evaluations per pixel. This is to cope with the case where the
            seeing is becoming small compared to the pixels, but obviously
            will slow things. To simply evaluate the profile once at the
            centre of each pixel in `wind`, set ndiv = 0.

    Returns:: (pars, epars, extras)

    where::

       pars    : tuple
          (sky, height, x, y, fwhm, beta), the fitted parameters
          (`beta` == 0 if `method`=='g')

       epars   : tuple
          (esky, eheight, ex, ey, efwhm, ebeta), fitted standard errors.
          `ebeta` == -1 if `method`=='g'; `efwhm` == -1 if the FWHM was fixed
          or it hit the minimum value.

       extras : tuple
          (fit,x,y,sigma,chisq,nok,nrej,npar,message) -- `fit` is an Window
          containing the best fit; `x` and `x` are the X and Y coordinates of
          the pixels, sigma are the final uncertainties used for each pixel,
          any negative were rejected; `chisq` is the chi**2 of the fit; `nok`
          is the number of points fitted; `nrej` is the number rejected;
          `npar` is the number of parameters fitted; `message` summarises the
          fit values. `x`, `y` and `sigma` are all 2D numpy arrays matching
          the dimensions of the data.

    Raises a HipercamError if the leastsq fails.

    """

    if method == 'g':
        # gaussian fit
        (sky, height, x, y, fwhm), \
            (esky, eheight, ex, ey, efwhm), \
            (fit, X, Y, sigma, chisq, nok, nrej, npar) = fitGaussian(
                wind, sky, height, x, y, fwhm, fwhm_min, fwhm_fix,
                read, gain, thresh, ndiv
            )

    elif method == 'm':
        # moffat fit
        (sky, height, x, y, fwhm, beta), \
            (esky, eheight, ex, ey, efwhm, ebeta), \
            (fit, X, Y, sigma, chisq, nok, nrej, npar) = fitMoffat(
                wind, sky, height, x, y, fwhm, fwhm_min, fwhm_fix,
                beta, read, gain, thresh, ndiv
            )

    else:
        raise NotImplementedError(
            '{:s} fitting method not implemented'.format(method)
        )

    if method == 'g':
        message = (
            'x,y = {:.1f}({:.1f}),{:.1f}({:.1f}),'
            ' FWHM = {:.2f}({:.2f}), peak = {:.1f}({:.1f}),'
            ' sky = {:.1f}({:.1f}), counts = {:.0f}, chi**2 = {:.1f}, nok = {:d}, nrej = {:d}'
        ).format(x,ex,y,ey,fwhm,efwhm,height,eheight,sky,esky,fit.sum(),chisq,nok,nrej)
        beta, ebeta = 0., -1.
    elif method == 'm':
        message = (
            'x,y = {:.1f}({:.1f}),{:.1f}({:.1f}),'
            ' FWHM = {:.2f}({:.2f}), peak = {:.1f}({:.1f}),'
            ' sky = {:.1f}({:.1f}), counts = {:.0f}, beta = {:.2f}({:.2f}), chi**2 = {:.1f},'
            ' nok = {:d}, nrej = {:d}'
        ).format(x,ex,y,ey,fwhm,efwhm,height,eheight,sky,esky,fit.sum(),beta,ebeta,chisq,nok,nrej)

    return (
        (sky, height, x, y, fwhm, beta),
        (esky, eheight, ex, ey, efwhm, ebeta),
        (fit, X, Y, sigma, chisq, nok, nrej, npar, message)
    )

##########################################
#
# 2D Moffat + constant section
#
##########################################

def fitMoffat(wind, sky, height, xcen, ycen, fwhm, fwhm_min, fwhm_fix,
              beta, read, gain, thresh, ndiv):
    """Fits the profile of one target in a Window with a symmetric 2D Moffat
    profile plus a constant "c + h/(1+alpha**2)**beta" where r is the distance
    from the centre of the aperture. The constant alpha is determined by the
    parameters `fwhm` and `beta`. The routine returns the fitted parameters,
    the covariances and a few extras including the best fit itself. The FWHM
    can be constrained to lie above a lower limit which can be useful when
    data are heavily binned. If when `fwhm_fix` == False, the first fit with a
    free FWHM fails, the routine will try again with the fwhm held fixed at
    its input value `fwhm`.

    Arguments::

        wind : :class:`Window`
            the Window with the data to be fitted. Ideally should contain
            only one target, so chopping down to a small region around
            the target of interest is usually a good idea.

        sky : float
            initial value of the (assumed constant) sky background (counts
            per pixel)

        height : float
            initial peak height of profile (counts)

        xcen : float
            initial X value at centre of profile (unbinned absolute
            coordinates)

        ycen : float
            initial Y value at centre of profile (unbinned absolute
            coordinates)

        fwhm : float
            initial FWHM (unbinned pixels).

        fwhm_min : float
            minimum value to allow FWHM to go to. Useful for heavily binned
            data especially where all signal might be in a single pixel to
            prevent going to silly values. Routine first fits with a free
            FWHM. If it goes below this value it tries again, setting to
            this value. (unbinned pixels)

        fwhm_fix : bool
            fix the FWHM at the value `fwhm` during fits. Should be a bit
            faster and more robust.

        read : float | array
            readout noise, from which weights are generated (RMS counts)

        gain : float | array
            gain, from which weights are generated (e- per count)

        thresh : float
            threshold in terms of RMS for rejection. The RMS is obtained from
            read and gain but scaled by the sqrt(chi**2/d.o.f). Typical value
            = 4

        ndiv : int
            Parameter controlling treatment of sub-pixellation. If ndiv > 0,
            every pixel will be sub-divided first into unbinned pixels, and
            then each unbinned pixels will be split into a square array of
            ndiv by ndiv points. The profile will be evaluated and averaged
            over each of these points. Thus if the pixels in `wind` are binned
            xbin by ybin, there will be xbin*ybin*ndiv**2 evaluations per
            pixel. This is to cope with the case where the seeing is becoming
            small compared to the pixels, but obviously will slow things. To
            simply evaluate the profile once at the centre of each pixel in
            `wind`, set ndiv = 0.

    Returns:: tuple of tuples

        (pars, sigs, extras) where::

           pars : tuple
                (sky, height, xcen, yce, fwhm, beta), the fitted parameters.
                fwhm = fwhm_min if initial fit returns an fwhm < fwhm_min,
                or == initial fwhm if fwhm_fix == True.

           sigs : tuple
                (skye, heighte, xcene, ycene, fwhme, betae), standard errors
                on the fit parameters. If fwhm defaults to `fwhm_min` or
                `fwhm_fix` == True, then `fwhme` will come back as -1.

           extras : tuple
                (fit, x, y, sigma, chisq, nok, nrej, npar) where: `fit` is an
                :class:`Window` containing the best fit; `x` and `y` are the x
                and y positions of all pixels (2D numpy arrays); `sigma` is
                the final set of RMS uncertainties on each pixel (2D numpy
                array); `chisq` is the raw chi**2; `nok` is the number of
                points fitted; `nrej` is the number rejected; `npar` is the
                number of parameters fitted (5 or 6).

    The program re-scales the uncertainties on the fit parameters by
    sqrt(chi**2/ndof) where ndof is the number of degrees of freedom = number
    of points - 5 or 6, depending on the fit. This allows for poor values
    of `read` and `gain` to a certain extent, which happen especially if a
    sky background has been removed at the start. The value of `sigma` on any
    rejected pixel is < 0.

    Raises a HipercamError if the leastsq fails.

    """

    # construct function objects of both types from the first because during
    # rejection, we assume both exist. sigma arrays are constructed here all > 0
    mfit1 = Mfit1(wind, read, gain, ndiv)
    dmfit1 = Dmfit1(mfit1)
    mfit2 = Mfit2(wind, read, gain, fwhm, ndiv)
    dmfit2 = Dmfit2(mfit2)

    while True:
        # Rejection loop. There is a break statement at the end of the loop
        # if no pixels have been rejected.

        # Two pass strategy: fit with a free FWHM, and if it is
        # below fwhm_min, re-fit with the fwhm set = fwhm_min

        if fwhm_fix:
            # FWHM held fixed
            mfit2.fwhm = fwhm
            dmfit2.fwhm = fwhm
            param = (sky, height, xcen, ycen, beta)

            # carry out fit
            soln, covar, info, mesg, ier = leastsq(
                mfit2, param, Dfun=dmfit2, col_deriv=True,
                full_output=True
            )
            if ier < 1 or ier > 4:
                raise HipercamError(mesg)
            elif covar is None:
                raise HipercamError('leastsq failed returning covar = None')

            # process results
            fit = Window(wind, mfit2.model(soln))
            skyf, heightf, xf, yf, betaf = soln
            covs = np.diag(covar)
            if (covs < 0).any():
                raise HipercamError('Negative covariance in fitMoffat')
            skyfe, heightfe, xfe, yfe, betafe = np.sqrt(covs)
            fwhmf, fwhmfe = fwhm, -1

        else:
            # FWHM free to vary
            param = (sky, height, xcen, ycen, fwhm, beta)

            # carry out fit
            soln, covar, info, mesg, ier = leastsq(
                mfit1, param, Dfun=dmfit1, col_deriv=True,
                full_output=True
            )
            if ier < 1 or ier > 4:
                raise HipercamError(mesg)
            elif covar is None:
                raise HipercamError('leastsq failed returning covar = None')

            # process results
            skyf, heightf, xf, yf, fwhmf, betaf = soln
            fwhmf = abs(fwhmf)

            if fwhmf > fwhm_min:
                # Free fit is OK
                fit = Window(wind, mfit1.model(soln))
                covs = np.diag(covar)
                if (covs < 0).any():
                    raise HipercamError('Negative covariance in fitMoffat')
                skyfe, heightfe, xfe, yfe, \
                    fwhmfe, betafe = np.sqrt(covs)

            else:
                # fall back to fixed FWHM
                mfit2.fwhm = fwhm_min
                dmfit2.fwhm = fwhm_min
                param = (sky, height, xcen, ycen, beta)

                # carry out fit
                soln, covar, info, mesg, ier = leastsq(
                    mfit2, param, Dfun=dmfit2, col_deriv=True,
                    full_output=True
                )
                if ier < 1 or ier > 4:
                    raise HipercamError(mesg)
                elif covar is None:
                    raise HipercamError('leastsq failed returning covar = None')

                # process results
                fit = Window(wind, mfit2.model(soln))
                skyf, heightf, xf, yf, betaf = soln
                covs = np.diag(covar)
                if (covs < 0).any():
                    raise HipercamError('Negative covariance in fitMoffat')
                skyfe, heightfe, xfe, yfe, betafe = np.sqrt(covs)
                fwhmf, fwhmfe = fwhm_min, -1

        # now look for bad outliers
        ok = mfit1.sigma > 0
        resid = (wind.data - fit.data) / mfit1.sigma
        chisq = (resid[ok]**2).sum()
        nok1 = len(resid[ok])
        sfac = np.sqrt(chisq/(nok1-len(soln)))

        # reject any above the defined threshold
        mfit1.sigma[ok & (np.abs(resid)> sfac*thresh)] *= -1

        # check whether any have been rejected
        ok = mfit1.sigma > 0
        nok = len(mfit1.sigma[ok])
        if nok < nok1:
            # some more pixels have been rejected
            mfit2.sigma = mfit1.sigma
            dmfit1.sigma = mfit1.sigma
            dmfit2.sigma = mfit1.sigma
        else:
            # no more pixels have been rejected.  calculate how many have
            # been, re-scale the uncertainties to reflect the actual chi**2
            nrej = mfit1.sigma.size - nok
            skyfe *= sfac
            heightfe *= sfac
            xfe *= sfac
            yfe *= sfac
            if fwhmfe > 0:
                fwhmfe *= sfac
            betafe *= sfac
            break

    # OK we are done.
    return (
        (skyf,heightf,xf,yf,fwhmf,betaf),
        (skyfe,heightfe,xfe,yfe,fwhmfe,betafe),
        (fit,mfit1.x,mfit1.y,mfit1.sigma,chisq,nok,nrej,len(soln))
    )

@jit(nopython=True,cache=True)
def moffat(x, y, sky, height, xcen, ycen, fwhm, beta, xbin, ybin, ndiv):
    """
    Returns a numpy array corresponding to the ordinate grids in xy
    set to a Moffat profile plus a constant. Defined by
    sky + height/(1+alpha*r**2)**beta where r is the distance from the
    centre and alpha is set to give the desired FWHM given the exponent
    beta. As beta becomes large, this tends to a Gaussian shape but has
    more extended wings at low beta.

    Parameters:

      x  : 2D numpy array
         the X ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      y  : 2D numpy array
         the Y ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      sky    : float
         sky background

      height : float
         height of central peak

      xcen   : float
         X-ordinate of centre

      ycen   : float
         Y-ordinate of centre

      fwhm   : float
         FWHM of the profile

      beta   : float
         Moffat exponent

      xbin   : int
         X-size of pixels in terms of unbinned pixels, i.e. the binning factor in X

      ybin   : int
         Y-size of pixels in terms of unbinned pixels, i.e. the binning factor in Y

      ndiv   : int
         Parameter controlling treatment of sub-pixellation. If > 0, every
         pixel will be sub-divided first into unbinned pixels, and then each
         unbinned pixels will be split into a square array of ndiv by ndiv
         points. The profile will be evaluated and averaged over each of these
         points. Thus if the pixels are binned xbin by ybin, there will be
         xbin*ybin*ndiv**2 evaluations per pixel. This is to cope with the
         case where the seeing is becoming small compared to the pixels, but
         obviously will slow things. To simply evaluate the profile once at
         the centre of each pixel, set ndiv = 0.

    Returns:: 2D numpy array containg the Moffat profile plus constant evaluated
    on the ordinate grids in xy.

    """
    tbeta = max(0.01,beta)
    alpha = 4*(2**(1./tbeta)-1)/fwhm**2

    if ndiv > 0:
        # Complicated case with sub-pixellation allowed for

        # mean offset within sub-pixels
        soff = (ndiv-1)/(2*ndiv)

        prof = np.zeros_like(x)
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
                        prof += (height/xbin/ybin/ndiv**2)*(1+alpha*rsq)**(-tbeta)

        return sky+prof

    else:
        # Fast as possible, compute profile at pixel centres only
        rsq = (x-xcen)**2+(y-ycen)**2
        return height*(1+alpha*rsq)**(-tbeta) + sky

@jit(nopython=True,cache=True)
def dmoffat(x, y, sky, height, xcen, ycen, fwhm, beta, xbin, ybin, ndiv,
            comp_dfwhm, comp_dbeta):
    """Returns a list of numpy arrays corresponding to the ordinate grids
    in xy set to the partial derivatives of a Moffat profile plus a
    constant. Defined by sky + height/(1+alpha*r**2)**beta where r is
    the distance from the centre and alpha is set to give the desired
    FWHM given the exponent beta. As beta becomes large, this tends to
    a Gaussian shape but has more extended wings at low beta. The
    partial derivatives are in the order of the parameters.

    Parameters:

      x  : 2D numpy array
         the X ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      y  : 2D numpy array
         the Y ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      sky    : float
         sky background

      height : float
         height of central peak

      xcen   : float
         X-ordinate of centre

      ycen   : float
         Y-ordinate of centre

      fwhm   : float
         FWHM of the profile

      beta   : float
         Moffat exponent

      xbin   : int
         X-size of pixels in terms of unbinned pixels, i.e. the binning factor in X

      ybin   : int
         Y-size of pixels in terms of unbinned pixels, i.e. the binning factor in Y

      ndiv   : int
         Parameter controlling treatment of sub-pixellation. If > 0, every
         pixel will be sub-divided first into unbinned pixels, and then each
         unbinned pixels will be split into a square array of ndiv by ndiv
         points. The profile will be evaluated and averaged over each of these
         points. Thus if the pixels are binned xbin by ybin, there will be
         xbin*ybin*ndiv**2 evaluations per pixel. This is to cope with the
         case where the seeing is becoming small compared to the pixels, but
         obviously will slow things. To simply evaluate the profile once at
         the centre of each pixel, set ndiv = 0.

      comp_dfwhm : bool
         compute derivative wrt FWHM (or not). If True, then derivatives with respect
         to (sky, height, xcen, ycen, fwhm, beta) are returned in a 6-element tuple.

      comp_dbeta : bool [only operative if comp_dfwhm=False]
         compute derivative wrt beta (or not). Only has an effect if comp_dfwhm=False.
         If True, 5 derivatives (sky, height, xcen, ycen, beta) are returned, if False
         4 derivatives (sky, height, xcen, ycen) come back.

    Returns:: in all cases a 6-element tuple is returned, but depending upon the comp flags
    only the first 4, 5 or 6 elements may be useful. 6-elements always come back because this
    helps the numba just-in-time compiler function better.

    """
    tbeta = max(0.01, beta)
    alpha = 4*(2**(1/tbeta)-1)/fwhm**2

    dsky = np.ones_like(x)

    if ndiv > 0:
        # complicated sub-pixellation case
        dheight = np.zeros_like(x)
        dxcen = np.zeros_like(x)
        dycen = np.zeros_like(x)
        if comp_dfwhm:
            dfwhm = np.zeros_like(x)
            dbeta = np.zeros_like(x)
        elif comp_dbeta:
            dbeta = np.zeros_like(x)

        # mean offset within sub-pixels
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

                        # finally compute stuff
                        dx = x + xsoff - xcen
                        dy = y + ysoff - ycen
                        rsq = dx**2 + dy**2

                        denom = 1+alpha*rsq
                        save1 = height*denom**(-tbeta-1)
                        save2 = save1*rsq

                        # derivatives. beta is a bit complicated because it appears directly
                        # through the exponent but also indirectly through alpha
                        dh = denom**(-tbeta)
                        dheight += dh
                        dxcen += (2*alpha*tbeta)*dx*save1
                        dycen += (2*alpha*tbeta)*dy*save1

                        if comp_dfwhm:
                            dfwhm += (2*alpha*tbeta/fwhm)*save2
                            dbeta += -np.log(denom)*height*dh + (4.*np.log(2)*2**(1/tbeta)/tbeta/fwhm**2)*save2
                        elif comp_dbeta:
                            dbeta += -np.log(denom)*height*dh + (4.*np.log(2)*2**(1/tbeta)/tbeta/fwhm**2)*save2

        # Normalise by number of evaluations
        nadd = xbin*ybin*ndiv**2
        dheight /= nadd
        dxcen /= nadd
        dycen /= nadd

        # in the next lines we return the same number of arrays in all cases
        # to help 'numba'
        if comp_dfwhm:
            # full set of derivs
            dfwhm /= nadd
            dbeta /= nadd
            return (dsky, dheight, dxcen, dycen, dfwhm, dbeta)
        elif comp_dbeta:
            # miss out dfwhm
            dbeta /= nadd
            return (dsky, dheight, dxcen, dycen, dbeta, dbeta)
        else:
            # miss out dbeta well
            return (dsky, dheight, dxcen, dycen, dycen, dycen)

    else:
        # fast as possible, only compute at centre of pixels
        dx = x - xcen
        dy = y - ycen
        rsq = dx**2 + dy**2

        denom = 1+alpha*rsq
        save1 = height*denom**(-tbeta-1)
        save2 = save1*rsq

        # derivatives. beta is a bit complicated because it appears directly
        # through the exponent but also indirectly through alpha
        dsky = np.ones_like(x)
        dheight = denom**(-tbeta)
        dxcen = (2*alpha*tbeta)*dx*save1
        dycen = (2*alpha*tbeta)*dy*save1

        if comp_dfwhm:
            dfwhm = (2*alpha*tbeta/fwhm)*save2
            dbeta = -np.log(denom)*height*dheight + (4.*np.log(2)*2**(1/tbeta)/tbeta/fwhm**2)*save2
            return (dsky, dheight, dxcen, dycen, dfwhm, dbeta)
        elif comp_dbeta:
            dbeta = -np.log(denom)*height*dheight + (4.*np.log(2)*2**(1/tbeta)/tbeta/fwhm**2)*save2
            return (dsky, dheight, dxcen, dycen, dbeta, dbeta)
        else:
            return (dsky, dheight, dxcen, dycen, dycen, dycen)

class Mfit1:
    """
    Function object to pass to leastsq for Moffat + constant
    background model with free FWHM.
    """

    def __init__(self, wind, read, gain, ndiv):
        """
        Arguments::

          wind  : Window
             the Window containing data to fit

          read  : float | array
             readout noise in RMS counts. Can be a 2D array if dimensions
             same as the data in wind

          gain  : float | array
             gain in electrons / count. Can be a 2D array if dimensions
             same as the data in wind

          ndiv  : int
             pixel sub-division factor. See comments in fitMoffat
        """
        self.sigma = np.sqrt(read**2+np.maximum(0,wind.data)/gain)
        x = wind.x(np.arange(wind.nx))
        y = wind.y(np.arange(wind.ny))
        self.x, self.y = np.meshgrid(x, y)
        self.data = wind.data
        self.xbin = wind.xbin
        self.ybin = wind.ybin
        self.ndiv = ndiv

    def __call__(self, param):
        """
        Returns 1D array of normalised residuals. See the model
        method for a description of the argument 'param'
        """
        mod = self.model(param)
        diff = (self.data-mod)/self.sigma
        ok = self.sigma > 0
        return diff[ok].ravel()

    def model(self, param):
        """
        Returns 2D array with model given a parameter vector.

        Argument::

           param : 1D array
              unpacks to (sky, height, xcen, ycen, fwhm, beta) where:
              'sky' is the background per pixel; 'height' is the central
              height of the Moffat function; 'xcen' and 'ycen' are the ordinates
              of its centre in unbinned CCD pixels with (1,1) at the left
              corner of the physical imagine area; 'fwhm' is the FWHM in unbinned
              pixels; 'beta' is the Moffat exponent.
        """
        sky, height, xcen, ycen, fwhm, beta = param
        return moffat(
            self.x, self.y, sky, height, xcen, ycen, fwhm, beta,
            self.xbin, self.ybin, self.ndiv
        )

class Dmfit1:
    """
    Function object to pass to leastsq to calculate the Jacobian equivalent
    to Mfit1
    """

    def __init__(self, mfit1):
        """
        Arguments::

          mfit1 : Mfit1
             the Mfit1 object passed as 'func' to leastsq
        """
        self.sigma = mfit1.sigma
        self.x, self.y = mfit1.x, mfit1.y
        self.xbin = mfit1.xbin
        self.ybin = mfit1.ybin
        self.ndiv = mfit1.ndiv

    def __call__(self, param):
        """
        Returns list of 1D arrays of the partial derivatives of
        the normalised residuals with respect to the variable
        parameters.
        """
        sky, height, xcen, ycen, fwhm, beta = param
        derivs = dmoffat(
            self.x, self.y, sky, height, xcen, ycen, fwhm, beta,
            self.xbin, self.ybin, self.ndiv, True, True
        )
        ok = self.sigma > 0
        return [(-deriv[ok]/self.sigma[ok]).ravel() for deriv in derivs]

class Mfit2:
    """
    Function object to pass to leastsq for Moffat + constant
    background model with fixed FWHM
    """

    def __init__(self, wind, read, gain, fwhm, ndiv):
        """
        Arguments:

          wind  : Window
             the Window containing data to fit

          read  : float / array
             readout noise in RMS counts. Can be a 2D array if dimensions
             same as the data in wind

          gain  : float / array
             gain in electrons / count. Can be a 2D array if dimensions
             same as the data in wind

          fwhm  : float
             the fixed FWHM to use

          ndiv  : int
             pixel sub-division factor. See comments in fitMoffat
        """
        self.sigma = np.sqrt(read**2+np.maximum(0,wind.data)/gain)
        x = wind.x(np.arange(wind.nx))
        y = wind.y(np.arange(wind.ny))
        self.x, self.y = np.meshgrid(x, y)
        self.data = wind.data
        self.fwhm = fwhm
        self.xbin = wind.xbin
        self.ybin = wind.ybin
        self.ndiv = ndiv

    def __call__(self, param):
        """
        Returns 1D array of normalised residuals
        """
        mod = self.model(param)
        diff = (self.data-mod)/self.sigma
        ok = self.sigma > 0
        return diff[ok].ravel()

    def model(self, param):
        """
        Returns 2D array with model given a parameter vector.

        Argument::

           param : 1D array
              unpacks to (sky, height, xcen, ycen, beta) where:
              'sky' is the background per pixel; 'height' is the central
              height of the Moffat function; 'xcen' and 'ycen' are the ordinates
              of its centre in unbinned CCD pixels with (1,1) at the left
              corner of the physical imagine area; 'beta' is the Moffat exponent.
        """
        sky, height, xcen, ycen, beta = param
        return moffat(
            self.x, self.y, sky, height, xcen, ycen, self.fwhm,
            beta, self.xbin, self.ybin, self.ndiv
        )

class Dmfit2:
    """
    Function object to pass to leastsq to calculate the Jacobian equivalent
    to Mfit2
    """

    def __init__(self, mfit2):
        """
        Arguments::

          mfit2 : Mfit2
             the Mfit2 object passed as 'func' to leastsq
        """
        self.sigma = mfit2.sigma
        self.x, self.y = mfit2.x, mfit2.y
        self.fwhm = mfit2.fwhm
        self.xbin = mfit2.xbin
        self.ybin = mfit2.ybin
        self.ndiv = mfit2.ndiv

    def __call__(self, param):
        """
        Returns list of 1D arrays of the partial derivatives of
        the normalised residuals with respect to the variable
        parameters.
        """
        sky, height, xcen, ycen, beta = param
        derivs = dmoffat(
            self.x, self.y, sky, height, xcen, ycen, self.fwhm, beta,
            self.xbin, self.ybin, self.ndiv, False, True
        )
        ok = self.sigma > 0
        return [(-deriv[ok]/self.sigma[ok]).ravel() for deriv in derivs[:-1]]

class Mfit3:
    """
    Function object to pass to leastsq for Moffat + constant
    background model with fixed FWHM and beta
    """

    def __init__(self, wind, read, gain, fwhm, beta, ndiv):
        """
        Arguments:

          wind  : Window
             the Window containing data to fit

          read  : float / array
             readout noise in RMS counts. Can be a 2D array if dimensions
             same as the data in wind

          gain  : float / array
             gain in electrons / count. Can be a 2D array if dimensions
             same as the data in wind

          fwhm  : float
             the fixed FWHM to use

          beta  : float
             the fixed FWHM to use

          ndiv  : int
             pixel sub-division factor. See comments in fitMoffat
        """
        self.sigma = np.sqrt(read**2+np.maximum(0,wind.data)/gain)
        x = wind.x(np.arange(wind.nx))
        y = wind.y(np.arange(wind.ny))
        self.x, self.y = np.meshgrid(x, y)
        self.data = wind.data
        self.fwhm = fwhm
        self.beta = beta
        self.xbin = wind.xbin
        self.ybin = wind.ybin
        self.ndiv = ndiv

    def __call__(self, param):
        """
        Returns 1D array of normalised residuals
        """
        mod = self.model(param)
        diff = (self.data-mod)/self.sigma
        ok = self.sigma > 0
        return diff[ok].ravel()

    def model(self, param):
        """
        Returns 2D array with model given a parameter vector.

        Argument::

           param : 1D array
              unpacks to (sky, height, xcen, ycen) where:
              'sky' is the background per pixel; 'height' is the central
              height of the Moffat function; 'xcen' and 'ycen' are the ordinates
              of its centre in unbinned CCD pixels with (1,1) at the left
              corner of the physical imagine area.
        """
        sky, height, xcen, ycen = param
        return moffat(
            self.x, self.y, sky, height, xcen, ycen, self.fwhm,
            self.beta, self.xbin, self.ybin, self.ndiv
        )

class Dmfit3:
    """
    Function object to pass to leastsq to calculate the Jacobian equivalent
    to Mfit3
    """

    def __init__(self, mfit3):
        """
        Arguments::

          mfit3 : Mfit3
             the Mfit3 object passed as 'func' to leastsq
        """
        self.sigma = mfit3.sigma
        self.x, self.y = mfit3.x, mfit3.y
        self.fwhm = mfit3.fwhm
        self.beta = mfit3.beta
        self.xbin = mfit3.xbin
        self.ybin = mfit3.ybin
        self.ndiv = mfit3.ndiv

    def __call__(self, param):
        """
        Returns list of 1D arrays of the partial derivatives of
        the normalised residuals with respect to the variable
        parameters.
        """
        sky, height, xcen, ycen = param
        derivs = dmoffat(
            self.x, self.y, sky, height, xcen, ycen, self.fwhm, self.beta,
            self.xbin, self.ybin, self.ndiv, False, False
        )
        ok = self.sigma > 0
        return [(-deriv[ok]/self.sigma[ok]).ravel() for deriv in derivs[:-2]]

##########################################
#
# 2D Gaussian + constant section
#
##########################################

def fitGaussian(
        wind, sky, height, xcen, ycen, fwhm, fwhm_min, fwhm_fix,
        read, gain, thresh, ndiv):
    """Fits the profile of one target in an Window with a 2D symmetric Gaussian
    profile "c + h*exp(-alpha*r**2)" where r is the distance from the centre
    of the aperture. The constant alpha is fixed by the FWHM. The function
    returns the fitted parameters, covariances and a few extras; see below.
    The FWHM can be constrained to lie above a lower limit which can be useful
    when data are heavily binned. If when fwhm_fix == False, the first fit
    with a free FWHM fails, the routine will try again with the fwhm held
    fixed at its input value `fwhm`.

    Arguments::

        wind : Window
            the Window with the data to be fitted. Ideally should contain only
            one target, so chopping down to a small region around the target
            of interest is usually a good idea.

        sky : float
            initial value of the (assumed constant) sky background (counts per
            pixel)

        height : float
            initial peak height of profile (counts)

        xcen : float
            initial X value at centre of profile (unbinned absolute
            coordinates)

        ycen : float
            initial Y value at centre of profile (unbinned absolute
            coordinates)

        fwhm : float
            initial FWHM (unbinned pixels).

        fwhm_min : float
            minimum value to allow FWHM to go to. Useful for heavily binned
            data especially where all signal might be in a single pixel to
            prevent going to silly values. Routine first fits with a free
            FWHM. If it goes below this value it tries again, setting to this
            value. (unbinned pixels)

        fwhm_fix : bool
            fix the FWHM at the value 'fwhm' during fits. Should be a bit
            faster and more robust.

        read : float / array
            readout noise, from which weights are generated. If an
            array it must match the dimensions of the Window. (RMS counts)

        gain : float / array
            gain, from which weights are generated. If an array it
            must match the dimensions of the Window (e- per count)

        thresh : float
            threshold in terms of RMS for rejection. The RMS is obtained from
            read and gain but scaled by the sqrt(chi**2/d.o.f). Typical value
            = 4

        ndiv : int
            Parameter controlling treatment of sub-pixellation. If ndiv > 0,
            every pixel will be sub-divided first into unbinned pixels, and
            then each unbinned pixels will be split into a square array of
            ndiv by ndiv points. The profile will be evaluated and averaged
            over each of these points. Thus if the pixels in `wind` are binned
            xbin by ybin, there will be xbin*ybin*ndiv**2 evaluations per
            pixel. This is to cope with the case where the seeing is becoming
            small compared to the pixels, but obviously will slow things. To
            simply evaluate the profile once at the centre of each pixel in
            `wind`, set ndiv = 0.

    Returns:: tuple of tuples

        (pars, sigs, extras) where::

           pars : tuple
                (sky, height, xcen, yce, fwhm). Fitted parameters.
                `fwhm` == `fwhm_min` if fwhm after first fit is < `fwhm_min`,
                or == initial fwhm if `fwhm_fix` == True.

           sigs : tuple
                (skye, heighte, xcene, ycene, fwhme). Standard errors on
                the fit parameters. If `fwhm` defaults to `fwhm_min` or
                `fwhm_fix` == True, then `fwhme` will come back as -1.

           extras : tuple
                (fit, x, y, sigma) where `fit` is a :class:`Window` containing
                the best fit, `x` and `y` are the x and y positions of all
                pixels (2D numpy arrays) and `sigma` is the final set of RMS
                uncertainties on each pixels [option for future: some may come
                back < 0 indicating rejection].

    Raises a HipercamError if the leastsq fails.

    """

    # construct function objects of both types from the first because during
    # rejection, we assume both exist. sigma arrays are constructed here all > 0
    gfit1 = Gfit1(wind, read, gain, ndiv)
    dgfit1 = Dgfit1(gfit1)
    gfit2 = Gfit2(wind, read, gain, fwhm, ndiv)
    dgfit2 = Dgfit2(gfit2)

    while True:
        # Rejection loop. There is a break statement at the end of the loop
        # if no pixels have been rejected.

        # Two pass strategy: fit with a free FWHM, and if it is
        # below fwhm_min, re-fit with the fwhm set = fwhm_min

        if fwhm_fix:
            # FWHM held fixed
            gfit2.fwhm = fwhm
            dgfit2.fwhm = fwhm
            param = (sky, height, xcen, ycen)

            # carry out fit
            soln, covar, info, mesg, ier = leastsq(
                gfit2, param, Dfun=dgfit2, col_deriv=True,
                full_output=True
            )
            if ier < 1 or ier > 4:
                raise HipercamError(mesg)
            elif covar is None:
                raise HipercamError('leastsq failed returning covar = None')

            # process results
            fit = Window(wind, gfit2.model(soln))
            skyf, heightf, xf, yf = soln
            covs = np.diag(covar)
            if (covs < 0).any():
                raise HipercamError('Negative covariance in fitGaussian')
            skyfe, heightfe, xfe, yfe = np.sqrt(covs)
            fwhmf, fwhmfe = fwhm, -1

        else:
            # FWHM free to vary
            param = (sky, height, xcen, ycen, fwhm)

            # carry out fit
            soln, covar, info, mesg, ier = leastsq(
                gfit1, param, Dfun=dgfit1, col_deriv=True,
                full_output=True
            )
            if ier < 1 or ier > 4:
                raise HipercamError(mesg)
            elif covar is None:
                raise HipercamError('leastsq failed returning covar = None')

            # process results
            skyf, heightf, xf, yf, fwhmf = soln
            fwhmf = abs(fwhmf)

            if fwhmf > fwhm_min:
                # Free fit is OK
                covs = np.diag(covar)
                if (covs < 0).any():
                    raise HipercamError('Negative covariance in fitGaussian')
                skyfe, heightfe, xfe, yfe, fwhmfe = np.sqrt(covs)
                fit = Window(wind, gfit1.model(soln))

            else:

                # fall back to fixed FWHM
                gfit2.fwhm = fwhm_min
                dgfit2.fwhm = fwhm_min
                param = (sky, height, xcen, ycen)

                # carry out fit
                soln, covar, info, mesg, ier = leastsq(
                    gfit2, param, Dfun=dgfit2, col_deriv=True,
                    full_output=True
                )
                if ier < 1 or ier > 4:
                    raise HipercamError(mesg)
                elif covar is None:
                    raise HipercamError('leastsq failed returning covar = None')

                # process results
                fit = Window(wind, gfit2.model(soln))
                skyf, heightf, xf, yf = soln
                covs = np.diag(covar)
                if (covs < 0).any():
                    raise HipercamError('Negative covariance in fitGaussian')
                skyfe, heightfe, xfe, yfe = np.sqrt(covs)
                fwhmf, fwhmfe = fwhm_min, -1

        # now look for bad outliers
        ok = gfit1.sigma > 0
        resid = (wind.data - fit.data) / gfit1.sigma
        chisq = (resid[ok]**2).sum()
        nok1 = len(resid[ok])
        sfac = np.sqrt(chisq/(nok1-len(soln)))

        # reject any above the defined threshold
        gfit1.sigma[ok & (np.abs(resid)> sfac*thresh)] *= -1

        # check whether any have been rejected
        ok = gfit1.sigma > 0
        nok = len(gfit1.sigma[ok])
        if nok < nok1:
            # some more pixels have been rejected
            gfit2.sigma = gfit1.sigma
            dgfit1.sigma = gfit1.sigma
            dgfit2.sigma = gfit1.sigma
        else:
            # no more pixels have been rejected.  calculate how many have
            # been, re-scale the uncertainties to reflect the actual chi**2
            nrej = gfit1.sigma.size - nok
            skyfe *= sfac
            heightfe *= sfac
            xfe *= sfac
            yfe *= sfac
            if fwhmfe > 0:
                fwhmfe *= sfac
            break

    # OK we are done.
    return (
        (skyf,heightf,xf,yf,fwhmf),
        (skyfe,heightfe,xfe,yfe,fwhmfe),
        (fit,gfit1.x,gfit1.y,gfit1.sigma,chisq,nok,nrej,len(soln))
    )

@jit(nopython=True,cache=True)
def gaussian(x, y, sky, height, xcen, ycen, fwhm, xbin, ybin, ndiv):
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

    alpha = 4.*np.log(2.)/fwhm**2

    if ndiv > 0:
        # Complicated case with sub-pixellation allowed for

        # mean offset within sub-pixels
        soff = (ndiv-1)/(2*ndiv)

        prof = np.zeros_like(x)
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

@jit(nopython=True,cache=True)
def dgaussian(x, y, sky, height, xcen, ycen, fwhm, xbin, ybin, ndiv):
    """Returns a list of four or five numpy arrays corresponding to the ordinate
    grids in xy set to the partial derivatives of a symmetric 2D Gaussian plus
    a constant. Defined by sky + height*exp(-alpha*r**2) where r is the
    distance from the centre and alpha is set to give the desired FWHM.  The
    partial derivatives are in the order of the parameters.

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
         FWHM of the profile

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

    Returns:: a list of five 2D numpy arrays containing the partial
    derivatives of a symmetric 2D Gaussian plus constant evaluated on
    the ordinate grids in xy. They come in the same order as they
    appear in the function call.

    """
    alpha = 4.*np.log(2.)/fwhm**2

    dsky = np.ones_like(x)

    if ndiv > 0:
        # complicated sub-pixellation case
        dheight = np.zeros_like(x)
        dxcen = np.zeros_like(x)
        dycen = np.zeros_like(x)
        dfwhm = np.zeros_like(x)

        # mean offset within sub-pixels
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

                        # finally compute stuff
                        dx = x + xsoff - xcen
                        dy = y + ysoff - ycen
                        rsq = dx**2 + dy**2

                        dh = np.exp(-alpha*rsq)
                        dheight += dh
                        dxcen += (2*alpha*height)*dh*dx
                        dycen += (2*alpha*height)*dh*dy
                        dfwhm += (2*alpha*height/fwhm)*dh*rsq

        # Normalise by number of evaluations
        nadd = xbin*ybin*ndiv**2
        dheight /= nadd
        dxcen /= nadd
        dycen /= nadd

        dfwhm /= nadd
        #            return np.dstack((dsky, dheight, dxcen, dycen, dfwhm))
        return (dsky, dheight, dxcen, dycen, dfwhm)

    else:
        # fast as possible, only compute at centre of pixels
        dx = x - xcen
        dy = y - ycen
        rsq = dx**2 + dy**2

        dheight = np.exp(-alpha*rsq)
        dxcen = (2*alpha*height)*dheight*dx
        dycen = (2*alpha*height)*dheight*dy
        dfwhm = (2*alpha*height/fwhm)*dheight*rsq
        #            return np.dstack((dsky, dheight, dxcen, dycen, dfwhm))
        return (dsky, dheight, dxcen, dycen, dfwhm)

class Gfit1:
    """
    Function object to pass to leastsq for Gaussian + constant
    background model with free FWHM.
    """

    def __init__(self, wind, read, gain, ndiv):
        """
        Arguments::

          wind : Window
             the Window containing data to fit

          read : float / array
             readout noise in RMS counts. Can be a 2D array if dimensions
             same as the data in wind

          gain : float / array
             gain in electrons / count. Can be a 2D array if dimensions
             same as the data in wind

          ndiv : int
             pixel sub-division factor. See comments in fitGaussian.
        """
        self.sigma = np.sqrt(read**2+np.maximum(0,wind.data)/gain)
        x = wind.x(np.arange(wind.nx))
        y = wind.y(np.arange(wind.ny))
        self.x, self.y = np.meshgrid(x, y)
        self.data = wind.data
        self.xbin = wind.xbin
        self.ybin = wind.ybin
        self.ndiv = ndiv

    def __call__(self, param):
        """
        Returns 1D array of normalised residuals. See the model
        method for a description of the argument 'param'
        """
        mod = self.model(param)
        diff = (self.data-mod)/self.sigma
        ok = self.sigma > 0
        return diff[ok].ravel()

    def model(self, param):
        """Returns 2D array with model given a parameter vector.

        Argument::

           param : 1D array
              unpacks to (sky, height, xcen, ycen, fwhm) where: 'sky' is the
              background per pixel; 'height' is the central height of the
              Gaussian; 'xcen' and 'ycen' are the ordinates of its centre in
              unbinned CCD pixels with (1,1) at the left corner of the
              physical imagine area; 'fwhm' is the FWHM in unbinned pixels.

        """
        sky, height, xcen, ycen, fwhm = param
        return support.gaussian(
            self.x, self.y, sky, height, xcen, ycen, fwhm,
            self.xbin, self.ybin, self.ndiv
        )

class Dgfit1:
    """
    Function object to pass to leastsq to calculate the Jacobian equivalent
    to Gfit1
    """

    def __init__(self, gfit1):
        """
        Arguments::

          gfit1 : Gfit1
             the Gfit1 object passed as 'func' to leastsq
        """
        self.sigma = gfit1.sigma
        self.x, self.y = gfit1.x, gfit1.y
        self.xbin = gfit1.xbin
        self.ybin = gfit1.ybin
        self.ndiv = gfit1.ndiv

    def __call__(self, param):
        """Returns list of 1D arrays of the partial derivatives of the normalised
        residuals with respect to the variable parameters.

        """
        sky, height, xcen, ycen, fwhm = param
        derivs = dgaussian(
            self.x, self.y, sky, height, xcen, ycen, fwhm,
            self.xbin, self.ybin, self.ndiv
        )
        ok = self.sigma > 0
        return [(-deriv[ok]/self.sigma[ok]).ravel() for deriv in derivs]

class Gfit2:
    """
    Function object to pass to leastsq for Gaussian + constant
    background model with a fixed FWHM.
    """

    def __init__(self, wind, read, gain, fwhm, ndiv):
        """
        Arguments:

          wind : Window
             the Window containing data to fit

          read : float / array
             readout noise in RMS counts. Can be a 2D array if dimensions
             same as the data in wind

          gain : float / array
             gain in electrons / count. Can be a 2D array if dimensions
             same as the data in wind

          fwhm : float
             the fixed FWHM to use

          ndiv : int
             sub-pixellation factor. See fitGaussian for more.

        """
        self.sigma = np.sqrt(read**2+np.maximum(0,wind.data)/gain)
        x = wind.x(np.arange(wind.nx))
        y = wind.y(np.arange(wind.ny))
        self.x, self.y = np.meshgrid(x, y)
        self.data = wind.data
        self.fwhm = fwhm
        self.xbin = wind.xbin
        self.ybin = wind.ybin
        self.ndiv = ndiv

    def __call__(self, param):
        """
        Returns 1D array of normalised residuals
        """

        mod = self.model(param)
        diff = (self.data-mod)/self.sigma
        return diff[self.sigma > 0].ravel()

    def model(self, param):
        """Returns 2D array with model given a parameter vector.

        Argument::

           param : 1D array
              unpacks to (sky, height, xcen, ycen) where: 'sky' is the
              background per pixel; 'height' is the central height of the
              Moffat function; 'xcen' and 'ycen' are the ordinates of its
              centre in unbinned CCD pixels with (1,1) at the left corner of
              the physical imagine area.

        """
        sky, height, xcen, ycen = param
        return gaussian(
            self.x, self.y, sky, height, xcen, ycen, self.fwhm,
            self.xbin, self.ybin, self.ndiv
        )

class Dgfit2:
    """
    Function object to pass to leastsq to calculate the Jacobian equivalent
    to Gfit2
    """

    def __init__(self, gfit2):
        """
        Arguments::

          gfit2 : Gfit2
             the Gfit2 object passed as 'func' to leastsq
        """
        self.sigma = gfit2.sigma
        self.x, self.y = gfit2.x, gfit2.y
        self.fwhm = gfit2.fwhm
        self.xbin = gfit2.xbin
        self.ybin = gfit2.ybin
        self.ndiv = gfit2.ndiv

    def __call__(self, param):
        """
        Returns list of 1D arrays of the partial derivatives of
        the normalised residuals with respect to the variable
        parameters.
        """
        sky, height, xcen, ycen = param
        derivs = dgaussian(
            self.x, self.y, sky, height, xcen, ycen, self.fwhm,
            self.xbin, self.ybin, self.ndiv
        )

        # note one more derivative is resturned than we need, so
        # we cut it out (the last one)
        ok = self.sigma > 0
        sigs = self.sigma[ok]
        return [(-deriv[ok]/sigs).ravel() for deriv in derivs[:-1]]
