# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Code for profile fitting. Currently supports symmetric 2D Gaussian and
Moffat profiles plus constants.
"""

from numba import jit
import numpy as np
from scipy.optimize import least_squares
from .core import *
from .window import *
from . import support

__all__ = ("combFit", "fitMoffat", "fitGaussian", "moffat", "gaussian")


def combFit(
        wind,
        sigma,
        method,
        sky,
        height,
        x,
        y,
        fwhm,
        fwhm_min,
        fwhm_fix,
        beta,
        beta_max,
        beta_fix,
        thresh,
        ndiv=0,
        max_nfev=None
):
    """Fits a stellar profile in a :class:Window using either a 2D Gaussian
    or Moffat profile. This is a convenience wrapper of fitMoffat and
    fitGaussian because one often wants both options.

    Arguments::

        wind : :class:`Window`
            the Window containing the stellar profile to fit.

        sigma : np.array
            array of sigma estimates matching dimensions of wind. Modified
            on exit. Rejected points will be negative.

        method : string
            fitting method 'g' for Gaussian, 'm' for Moffat

        sky : float
            initial sky level

        height : float
            initial peak height

        x : float
            initial central X value

        y : float
            initial central Y value

        fwhm : float
            initial FWHM, unbinned pixels.

        fwhm_min : float
            minimum FWHM, unbinned pixels.

        fwhm_fix : bool
            fix the FWHM (i.e. don't fit it)

        beta : float [if method == 'm']
            exponent of Moffat function

        beta_max : float [if method == 'm']
            maximum beta

        beta_fix : bool [if method == 'm']
            fix beta (i.e. don't fit it)

        thresh : float
            threshold in terms of RMS for rejection. The RMS is obtained from
            read and gain but scaled by the sqrt(chi**2/d.o.f). Typical value
            = 4. The first fit is carried out with 2*thresh as the rejection
            threshold to avoid spurious rejections.

        ndiv : int
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

       pars : tuple
          (sky, height, x, y, fwhm, beta), the fitted parameters
          (`beta` == 0 if `method`=='g')

       epars : tuple
          (esky, eheight, ex, ey, efwhm, ebeta), fitted standard errors.
          `ebeta` == -1 if `method`=='g'; `efwhm` == -1 if the FWHM was fixed
          or it hit the minimum value.

       extras : tuple
          (fit,x,y,chisq,nok,nrej,npar,message) -- `fit` is an Window
          containing the best fit; `x` and `x` are the X and Y coordinates of
          the pixels; `chisq` is the chi**2 of the fit; `nok`
          is the number of points fitted; `nrej` is the number rejected;
          `npar` is the number of parameters fitted; `message` summarises the
          fit values. `x`  and `y` are 2D numpy arrays matching
          the dimensions of the data.

    Raises a HipercamError if least_squares fails.

    """

    if method == "g":
        # gaussian fit
        (
            (sky, height, x, y, fwhm),
            (esky, eheight, ex, ey, efwhm),
            (fit, X, Y, chisq, nok, nrej, npar, nfev)
        ) = fitGaussian(
            wind, sigma, sky, height, x, y, fwhm, fwhm_min, fwhm_fix,
            thresh, ndiv, max_nfev
        )

    elif method == "m":
        # moffat fit
        (
            (sky, height, x, y, fwhm, beta),
            (esky, eheight, ex, ey, efwhm, ebeta),
            (fit, X, Y, chisq, nok, nrej, npar, nfev)
        ) = fitMoffat(
            wind, sigma, sky, height, x, y, fwhm, fwhm_min, fwhm_fix, beta,
            beta_max, beta_fix, thresh, ndiv, max_nfev
        )

    else:
        raise NotImplementedError(f"Fitting method = {method} method not implemented")

    if method == "g":
        message = (
            f"x,y= {x:.2f}({ex:.2f}), {y:.2f}({ey:.2f}),"
            f" FWHM= {fwhm:.2f}({efwhm:.2f}),"
            f" peak= {height:.1f}({eheight:.1f}),"
            f" sky= {sky:.1f}({esky:.1f}),"
            f" counts= {(fit-sky).sum():.0f}, chi**2= {chisq:.1f},"
            f" nok= {nok}, nrej= {nrej}, nfev= {nfev}"
        )
        beta, ebeta = 0.0, -1.0
    elif method == "m":
        message = (
            f"x,y= {x:.2f}({ex:.2f}), {y:.2f}({ey:.2f}),"
            f" FWHM= {fwhm:.2f}({efwhm:.2f}),"
            f" peak= {height:.1f}({eheight:.1f}),"
            f" sky= {sky:.1f}({esky:.1f}),"
            f" counts= {(fit-sky).sum():.0f}, "
            f" beta= {beta:.2f}({ebeta:.2f}),"
            f" chi**2= {chisq:.1f},"
            f" nok= {nok}, nrej= {nrej}, nfev= {nfev}"
        )

    return (
        (sky, height, x, y, fwhm, beta),
        (esky, eheight, ex, ey, efwhm, ebeta),
        (fit, X, Y, chisq, nok, nrej, npar, nfev, message),
    )


################################
#
# Symmetric 2D Moffat + constant
#
################################


def fitMoffat(
        wind,
        sigma,
        sky,
        height,
        xcen,
        ycen,
        fwhm,
        fwhm_min,
        fwhm_fix,
        beta,
        beta_max,
        beta_fix,
        thresh,
        ndiv,
        max_nfev=None,
):
    """Fits the profile of one target in a Window with a symmetric 2D Moffat
    profile plus a constant "c + h/(1+alpha**2)**beta" where r is the distance
    from the centre of the aperture. The constant alpha is determined by the
    parameters `fwhm` and `beta`. The routine returns the fitted parameters,
    the covariances and a few extras including the best fit itself. The FWHM
    can be constrained to lie above a lower limit which can be useful when
    data are heavily binned. If when `fwhm_fix` == False, the first fit with a
    free FWHM fails, the routine will try again with the fwhm held fixed at
    its input value `fwhm`. Similarly the Moffat exponent beta can be limited
    to a maximum value as Moffat profiles are degenerate with Gaussians for
    large beta.

    Arguments::

        wind : :class:`Window`
            the Window with the data to be fitted. Ideally should contain
            only one target, so chopping down to a small region around
            the target of interest is usually a good idea.

        sigma : np.array
            array of sigma estimates matching dimensions of wind. Modified on
            exit; rejected points are multiplied by -1

        sky : float or None
            initial value of the (assumed constant) sky background (counts
            per pixel). Set None to ignore (i.e. regard as = 0)

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

        beta : float
            initial value of Moffat exponent beta

        beta_max : float
            maximum value to allow beta to go to. For large beta,
            Moffat profiles degenerate into guassians. The routines
            tries with a free beta, but restarts if beta goes above
            this value.

        beta_fix : bool
            fix beta at its initial value during fits. Should be a bit
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

        max_nfev : int or None
           maximum number of function evaluations during fits. Passed
           direct to least_squares.

    Returns:: tuple

        (pars, sigs, extras) where::

           pars : tuple
                (sky, height, xcen, yce, fwhm, beta), the fitted parameters.
                fwhm = fwhm_min if initial fit returns an fwhm < fwhm_min,
                or == initial fwhm if fwhm_fix == True. sky = 0 if sky
                not fitted.

           sigs : tuple
                (skye, heighte, xcene, ycene, fwhme, betae), standard errors
                on the fit parameters. If fwhm defaults to `fwhm_min` or
                `fwhm_fix` == True, then `fwhme` will come back as -1. If
                the sky is not fitted, skye comes back as -1.

           extras : tuple
                (fit, x, y, chisq, nok, nrej, npar, nfev) where: `fit`
                is an :class:`Window` containing the best fit; `x` and
                `y` are the x and y positions of all pixels (2D numpy
                arrays); `chisq` is the raw chi**2; `nok` is the
                number of points fitted; `nrej` is the number
                rejected; `npar` is the number of parameters fitted (3
                to 6), `nfev` is the number of functions evaluations.

    The program re-scales the uncertainties on the fit parameters by
    sqrt(chi**2/ndof) where ndof is the number of degrees of freedom =
    number of points - 5 or 6, depending on the fit. This allows for
    poor values of `read` and `gain` to a certain extent, which happen
    especially if a sky background has been removed at the start. The
    value of `sigma` on any rejected pixel is < 0.

    Raises a HipercamError if the least_squares fails.

    """

    if fwhm < fwhm_min:
        raise HipercamError(f"fwhm < fwhm_min ({fwhm} < {fwhm_min})")

    if beta > beta_max:
        raise HipercamError(f"beta > beta_max ({beta} > {beta_max})")

    # construct function object for Moffat fit. First the mode.
    mode = "" if sky is None else "s"
    mode = mode if fwhm_fix else mode + "f"
    mode = mode if beta_fix else mode + "b"
    mfit = Mfit(wind, sigma, ndiv, mode, fwhm, beta)

    mode_switch = False

    first_fit = True
    nfev = 0
    while True:

        # retrieve current fit mode
        mode = mfit.mode

        # Rejection loop. There is a break statement at the end of the
        # loop if no pixels have been rejected.

        # set parameters
        param = mfit.set_par(sky, height, xcen, ycen, fwhm, beta)

        # carry out fit
        res = least_squares(
            mfit.fun, param, jac=mfit.jac, method="lm", max_nfev=max_nfev
        )
        if not res.success:
            raise HipercamError(res.message)
        nfev += res.nfev

        # get Jacobian
        J = np.matrix(res.jac)
        try:
            # try / except for singularity
            covar = (J.T * J).I
        except (np.linalg.LinAlgError, ValueError) as err:
            raise HipercamError(err)

        param = res.x
        skyf, heightf, xf, yf, fwhmf, betaf = mfit.get_par(param)
        fwhmf = abs(fwhmf)

        if mode.find("b") > -1 and betaf > beta_max:
            # switch to fixed beta fit
            beta = beta_max
            mfit.set_mode(mode.replace("b", ""), fwhm, beta)

        elif mode.find("f") > -1 and fwhmf < fwhm_min:
            # switch to fixed fwhm fit
            mfit.set_mode(mode.replace("f", ""), fwhm_min, beta)

        else:

            # no mode switch, compute fit
            fit = Window(wind, mfit.model(param))

            covs = np.diag(covar)
            if (covs < 0).any():
                raise HipercamError("Negative covariance in fitMoffat")

            # compute chi**2 and number of OK points
            ok = mfit.mask & (sigma > 0)
            resid = (wind.data - fit.data) / sigma
            chisq = (resid[ok] ** 2).sum()
            nok1 = len(resid[ok])
            sfac = np.sqrt(chisq / nok1)

            # reject any above the defined threshold
            if first_fit:
                # first fit carried out with higher threshold for safety
                sigma[ok & (np.abs(resid) > 2*sfac*thresh)] *= -1
            else:
                sigma[ok & (np.abs(resid) > 2*sfac*thresh)] *= -1

            # check whether any have been rejected
            ok = mfit.mask & (sigma > 0)
            nok = len(sigma[ok])

            if nok == nok1:
                # no more pixels have been rejected.  calculate how
                # many have been, re-scale the uncertainties to
                # reflect the actual chi**2
                nrej = len(sigma[mfit.mask]) - nok
                skyfe, heightfe, xfe, yfe, fwhmfe, betafe = mfit.get_epar(
                    sfac * np.sqrt(covs)
                )
                if not first_fit:
                    # always carry out at least 1
                    break
            first_fit = False

    # Make sure all regions masked off have sigma < 0
    sigma[~mfit.mask] = - np.abs(sigma[~mfit.mask])

    # OK we are done.
    extras = (
        fit,
        mfit.x,
        mfit.y,
        chisq,
        nok,
        nrej,
        len(param),
        nfev
    )

    if sky is None:
        return (
            (0., heightf, xf, yf, fwhmf, betaf),
            (-1., heightfe, xfe, yfe, fwhmfe, betafe), extras
        )
    else:
        return (
            (skyf, heightf, xf, yf, fwhmf, betaf),
            (skyfe, heightfe, xfe, yfe, fwhmfe, betafe), extras
        )


@jit(nopython=True, cache=True)
def moffat(x, y, sky, height, xcen, ycen, fwhm, beta, xbin, ybin, ndiv):
    """
    Returns a numpy array corresponding to the ordinate grids in xy
    set to a Moffat profile plus a constant. Defined by
    sky + height/(1+alpha*r**2)**beta where r is the distance from the
    centre and alpha is set to give the desired FWHM given the exponent
    beta. As beta becomes large, this tends to a Gaussian shape but has
    more extended wings at low beta.

    Parameters:

      x : 2D numpy array
         the X ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      y : 2D numpy array
         the Y ordinates over which you want to compute the
         profile. They should be measured in term of unbinned pixels.

      sky : float
         sky background.

      height : float
         height of central peak

      xcen : float
         X-ordinate of centre

      ycen : float
         Y-ordinate of centre

      fwhm : float
         FWHM of the profile

      beta : float
         Moffat exponent

      xbin : int
         X-size of pixels in terms of unbinned pixels, i.e. the binning factor in X

      ybin : int
         Y-size of pixels in terms of unbinned pixels, i.e. the binning factor in Y

      ndiv : int
         Parameter controlling treatment of sub-pixellation. If > 0, every
         pixel will be sub-divided first into unbinned pixels, and then each
         unbinned pixel will be split into a square array of ndiv by ndiv
         points. The profile will be evaluated and averaged over each of these
         points. Thus if the pixels are binned xbin by ybin, there will be
         xbin*ybin*ndiv**2 evaluations per pixel. This is to cope with the
         case where the seeing is becoming small compared to the pixels, but
         obviously will slow things. To simply evaluate the profile once at
         the centre of each pixel, set ndiv = 0.

    Returns:: 2D numpy array containg the Moffat profile plus constant evaluated
    on the ordinate grids in xy.

    """
    tbeta = max(0.01, beta)
    alpha = 4 * (2 ** (1.0 / tbeta) - 1) / fwhm ** 2

    if ndiv > 0:
        # Complicated case with sub-pixellation allowed for

        # mean offset within sub-pixels
        soff = (ndiv - 1) / (2 * ndiv)

        prof = np.zeros_like(x)
        for iy in range(ybin):
            # loop over unbinned pixels in Y
            yoff = iy - (ybin - 1) / 2 - soff
            for ix in range(xbin):
                # loop over unbinned pixels in X
                xoff = ix - (xbin - 1) / 2 - soff
                for isy in range(ndiv):
                    # loop over sub-pixels in y
                    ysoff = yoff + isy / ndiv
                    for isx in range(ndiv):
                        # loop over sub-pixels in x
                        xsoff = xoff + isx / ndiv
                        rsq = (x + xsoff - xcen) ** 2 + (y + ysoff - ycen) ** 2
                        prof += (height / xbin / ybin / ndiv ** 2) * (
                            1 + alpha * rsq
                        ) ** (-tbeta)

        return sky + prof

    else:
        # Fast as possible, compute profile at pixel centres only
        rsq = (x - xcen) ** 2 + (y - ycen) ** 2
        return height * (1 + alpha * rsq) ** (-tbeta) + sky


@jit(nopython=True, cache=True)
def dmoffat(
    x, y, sky, height, xcen, ycen, fwhm, beta, xbin, ybin, ndiv, comp_dfwhm, comp_dbeta
):
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

      sky : float or None
         sky background. If `None` it will be ignored.

      height : float
         height of central peak

      xcen : float
         X-ordinate of centre

      ycen : float
         Y-ordinate of centre

      fwhm : float
         FWHM of the profile

      beta : float
         Moffat exponent

      xbin : int
         X-size of pixels in terms of unbinned pixels, i.e. the
         binning factor in X

      ybin : int
         Y-size of pixels in terms of unbinned pixels, i.e. the
         binning factor in Y

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

      comp_dfwhm : bool
         compute derivative wrt FWHM (or not). If False, derivative with respect
         to FWHM will not be computed

      comp_dbeta : bool
         compute derivative wrt beta (or not). If False, derivative with respect
         to beta will not be computed.

    Returns:: in all cases a 6-element tuple is returned, but
    depending upon the comp flags only the first 4, 5 or 6 elements
    may be useful. 6-elements always come back because this helps the
    numba just-in-time compiler function better.

    """
    tbeta = max(0.01, beta)
    alpha = 4 * (2 ** (1 / tbeta) - 1) / fwhm ** 2

    dsky = np.ones_like(x)

    if ndiv > 0:
        # complicated sub-pixellation case
        dheight = np.zeros_like(x)
        dxcen = np.zeros_like(x)
        dycen = np.zeros_like(x)
        if comp_dfwhm:
            dfwhm = np.zeros_like(x)
        if comp_dbeta:
            dbeta = np.zeros_like(x)

        # mean offset within sub-pixels
        soff = (ndiv - 1) / (2 * ndiv)

        for iy in range(ybin):
            # loop over unbinned pixels in Y
            yoff = iy - (ybin - 1) / 2 - soff
            for ix in range(xbin):
                # loop over unbinned pixels in X
                xoff = ix - (xbin - 1) / 2 - soff
                for isy in range(ndiv):
                    # loop over sub-pixels in y
                    ysoff = yoff + isy / ndiv
                    for isx in range(ndiv):
                        # loop over sub-pixels in x
                        xsoff = xoff + isx / ndiv

                        # finally compute stuff
                        dx = x + xsoff - xcen
                        dy = y + ysoff - ycen
                        rsq = dx ** 2 + dy ** 2

                        denom = 1 + alpha * rsq
                        save1 = height * denom ** (-tbeta - 1)
                        save2 = save1 * rsq

                        # derivatives. beta is a bit complicated
                        # because it appears directly through the
                        # exponent but also indirectly through alpha
                        dh = denom ** (-tbeta)
                        dheight += dh
                        dxcen += (2 * alpha * tbeta) * dx * save1
                        dycen += (2 * alpha * tbeta) * dy * save1

                        if comp_dfwhm:
                            dfwhm += (2 * alpha * tbeta / fwhm) * save2

                        if comp_dbeta:
                            dbeta += (
                                -np.log(denom) * height * dh
                                + (
                                    4.0
                                    * np.log(2)
                                    * 2 ** (1 / tbeta)
                                    / tbeta
                                    / fwhm ** 2
                                )
                                * save2
                            )

        # Normalise by number of evaluations
        nadd = xbin * ybin * ndiv ** 2
        dheight /= nadd
        dxcen /= nadd
        dycen /= nadd

        # in the next lines we return the same number of arrays in all cases
        # to help 'numba'
        if comp_dfwhm and comp_dbeta:
            # full set of derivs
            dfwhm /= nadd
            dbeta /= nadd
            return (dsky, dheight, dxcen, dycen, dfwhm, dbeta)

        elif comp_dfwhm:
            dfwhm /= nadd
            return (dsky, dheight, dxcen, dycen, dfwhm, dfwhm)

        elif comp_dbeta:
            dbeta /= nadd
            return (dsky, dheight, dxcen, dycen, dbeta, dbeta)

        else:
            return (dsky, dheight, dxcen, dycen, dycen, dycen)

    else:
        # fast as possible, only compute at centre of pixels
        dx = x - xcen
        dy = y - ycen
        rsq = dx ** 2 + dy ** 2

        denom = 1 + alpha * rsq
        save1 = height * denom ** (-tbeta - 1)
        save2 = save1 * rsq

        # derivatives. beta is a bit complicated because it appears directly
        # through the exponent but also indirectly through alpha
        dheight = denom ** (-tbeta)
        dxcen = (2 * alpha * tbeta) * dx * save1
        dycen = (2 * alpha * tbeta) * dy * save1

        if comp_dfwhm and comp_dbeta:
            dfwhm = (2 * alpha * tbeta / fwhm) * save2
            dbeta = (
                -np.log(denom) * height * dheight
                + (4.0 * np.log(2) * 2 ** (1 / tbeta) / tbeta / fwhm ** 2) * save2
            )
            return (dsky, dheight, dxcen, dycen, dfwhm, dbeta)

        elif comp_dfwhm:
            dfwhm = (2 * alpha * tbeta / fwhm) * save2
            return (dsky, dheight, dxcen, dycen, dfwhm, dfwhm)

        elif comp_dbeta:
            dbeta = (
                -np.log(denom) * height * dheight
                + (4.0 * np.log(2) * 2 ** (1 / tbeta) / tbeta / fwhm ** 2) * save2
            )
            return (dsky, dheight, dxcen, dycen, dbeta, dbeta)

        else:
            return (dsky, dheight, dxcen, dycen, dycen, dycen)


def _mask(wind, x, y):
    """Returns circular mask centred on wind, extending to nearest side"""
    x1, x2 = wind.x(0), wind.x(wind.nx - 1)
    y1, y2 = wind.y(0), wind.y(wind.ny - 1)
    xc, yc = (x1 + x2) / 2, (y1 + y2) / 2.0
    rad = 1.01 * min((x2 - x1) / 2, (y2 - y1) / 2)
    return (x - xc) ** 2 + (y - yc) ** 2 < rad ** 2


class Mfit:
    """Object providing 'fun' and 'jac' methods for least_squares for
    Mofffat models. Eight operating modes each of which allows the
    following to be free [else not]. Sky assumed = 0 when not free.

       mode == 'sfb' : sky, FWHM, beta
       mode == 'sb' : sky, beta
       mode == 'sf' : sky, FWHM
       mode == 's' : sky
       mode == 'fb' : FWHM, beta
       mode == 'b' : beta
       mode == 'f' : FWHM
       mode == '' : None of the above.

    It is up to the user to set the values of fwhm and beta
    appropriate to the mode. A method 'set_mode' provides some basic
    checks.

    """

    def __init__(self, wind, sigma, ndiv, mode, fwhm=None, beta=None):
        """
        Arguments::

          wind : Window
             the Window containing data to fit

          sigma : np.array
             matching array of uncertainties

          ndiv  : int
             pixel sub-division factor. See comments in fitMoffat
        """
        self.sigma = sigma
        x = wind.x(np.arange(wind.nx))
        y = wind.y(np.arange(wind.ny))
        self.x, self.y = np.meshgrid(x, y)
        self.data = wind.data
        self.xbin = wind.xbin
        self.ybin = wind.ybin
        self.ndiv = ndiv
        self.mask = _mask(wind, self.x, self.y)
        self.set_mode(mode, fwhm, beta)

    def set_mode(self, mode, fwhm, beta):
        """Set the operation mode with some light checks"""
        if mode not in ("sfb", "sb", "sf", "s", "fb", "b", "f", ""):
            raise HipercamError("invalid mode = {:s}".format(mode))

        if (mode.find("f") == -1 and fwhm is None) or (
            mode.find("b") == -1 and beta is None
        ):
            raise HipercamError("invalid mode / fwhm / beta combination")
        self.mode = mode
        self.fwhm = fwhm
        self.beta = beta

    def get_par(self, param):
        """Gets parameters (sky, height, xcen, ycen, fwhm, beta) according to
        current least_squares parameter vector, accounting for the
        operating mode. All six are returned in each case, even when this
        is meaningless.

        """
        if self.mode == "sfb":
            sky, height, xcen, ycen, fwhm, beta = param
        elif self.mode == "sb":
            sky, height, xcen, ycen, beta = param
            fwhm = self.fwhm
        elif self.mode == "sf":
            sky, height, xcen, ycen, fwhm = param
            beta = self.beta
        elif self.mode == "s":
            sky, height, xcen, ycen = param
            fwhm = self.fwhm
            beta = self.beta
        elif self.mode == "fb":
            height, xcen, ycen, fwhm, beta = param
            sky = 0.0
        elif self.mode == "b":
            height, xcen, ycen, beta = param
            sky = 0.0
            fwhm = self.fwhm
        elif self.mode == "f":
            height, xcen, ycen, fwhm = param
            sky = 0.0
            beta = self.beta
        elif self.mode == "":
            height, xcen, ycen = param
            sky = 0.0
            fwhm = self.fwhm
            beta = self.beta
        else:
            raise HipercamError("invalid mode")

        return (sky, height, xcen, ycen, fwhm, beta)

    def set_par(self, sky, height, xcen, ycen, fwhm, beta):
        """Sets parameters (sky, height, xcen, ycen, fwhm, beta) into correct
        vector for least_squares according to the operating mode

        Argument::

           param : 1D array
              unpacks to (sky, height, xcen, ycen, fwhm, beta) where:
              'sky' is the background per pixel; 'height' is the
              central height of the Moffat function; 'xcen' and 'ycen'
              are the ordinates of its centre in unbinned CCD pixels
              with (1,1) at the left corner of the physical imagine
              area; 'fwhm' is the FWHM in unbinned

        """
        if self.mode == "sfb":
            return (sky, height, xcen, ycen, fwhm, beta)
        elif self.mode == "sb":
            return (sky, height, xcen, ycen, beta)
        elif self.mode == "sf":
            return (sky, height, xcen, ycen, fwhm)
        elif self.mode == "s":
            return (sky, height, xcen, ycen)
        elif self.mode == "fb":
            return (height, xcen, ycen, fwhm, beta)
        elif self.mode == "b":
            reurn(height, xcen, ycen, beta)
        elif self.mode == "f":
            return (height, xcen, ycen, fwhm)
        elif self.mode == "":
            return (height, xcen, ycen)
        else:
            raise HipercamError("invalid mode")

    def get_epar(self, param):
        """Gets parameter errors (skye, heighte, xcene, ycene, fwhme, betae)
        passed a vector of errors, accounting for the operating
        mode. All six are returned in each case, with fixed ones set =
        -1

        """
        if self.mode == "sfb":
            skye, heighte, xcene, ycene, fwhme, betae = param
        elif self.mode == "sb":
            skye, heighte, xcene, ycene, betae = param
            fwhme = -1.0
        elif self.mode == "sf":
            skye, heighte, xcene, ycene, fwhme = param
            betae = -1.0
        elif self.mode == "s":
            skye, heighte, xcene, ycene = param
            fwhme = -1.0
            betae = -1.0
        elif self.mode == "fb":
            heighte, xcene, ycene, fwhme, betae = param
            skye = -1.0
        elif self.mode == "b":
            heighte, xcene, ycene, betae = param
            skye = -1.0
            fwhme = -1.0
        elif self.mode == "f":
            heighte, xcene, ycene, fwhme = param
            skye = -1.0
            betae = -1.0
        elif self.mode == "":
            heighte, xcene, ycene = param
            skye = -1.0
            fwhme = -1.0
            betae = -1.0

        return (skye, heighte, xcene, ycene, fwhme, betae)

    def fun(self, param):
        """
        Returns 1D array of normalised residuals. See the model
        method for a description of the argument 'param'
        """
        mod = self.model(param)
        diff = (self.data - mod) / self.sigma
        ok = self.mask & (self.sigma > 0)
        return diff[ok].ravel()

    def jac(self, param):
        """
        Returns list of 1D arrays of the partial derivatives of
        the normalised residuals with respect to the variable
        parameters.
        """
        sky, height, xcen, ycen, fwhm, beta = self.get_par(param)

        comp_fwhm = self.mode.find("f") > -1
        comp_beta = self.mode.find("b") > -1

        # work out which derivatives to bother with
        if self.mode == "sfb":
            inds = (0, 1, 2, 3, 4, 5)
        elif self.mode == "sb" or self.mode == "sf":
            inds = (0, 1, 2, 3, 4)
        elif self.mode == "s":
            inds = (0, 1, 2, 3)
        elif self.mode == "fb":
            inds = (1, 2, 3, 4, 5)
        elif self.mode == "b" or self.mode == "f":
            inds = (1, 2, 3, 4)
        elif self.mode == "":
            inds = (1, 2, 3)

        derivs = dmoffat(
            self.x,
            self.y,
            sky,
            height,
            xcen,
            ycen,
            fwhm,
            beta,
            self.xbin,
            self.ybin,
            self.ndiv,
            comp_fwhm,
            comp_beta,
        )

        ok = self.mask & (self.sigma > 0)
        return np.column_stack(
            [(-derivs[ind][ok] / self.sigma[ok]).ravel() for ind in inds]
        )

    def model(self, param):
        """
        Returns 2D array with model given a parameter vector.

        Argument::

           param : 1D array
              parameter vector, with values that depend upon the mode,
              but could include some or all of (sky, height, xcen, ycen,
              fwhm, beta) where 'sky' is the background per pixel;
              'height' is the central height of the Moffat function;
              'xcen' and 'ycen' are the ordinates of its centre in
              unbinned CCD pixels with (1,1) at the left corner of the
              physical imaging area; 'fwhm' is the FWHM in unbinned
              pixels; 'beta' is the Moffat exponent.
        """
        sky, height, xcen, ycen, fwhm, beta = self.get_par(param)
        return moffat(
            self.x,
            self.y,
            sky,
            height,
            xcen,
            ycen,
            fwhm,
            beta,
            self.xbin,
            self.ybin,
            self.ndiv,
        )


##########################################
#
# 2D Gaussian + constant section
#
##########################################


def fitGaussian(
        wind,
        sigma,
        sky,
        height,
        xcen,
        ycen,
        fwhm,
        fwhm_min,
        fwhm_fix,
        thresh,
        ndiv,
        max_nfev=0,
):
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

        sigma : np.array
            array of uncertainties matching wind. Modified on exit. Rejected
            points multiplied by -1

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

        max_nfev : int
            maximum number of function evaluations during fits. Passed directly
            to leastsq.

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
                (fit, x, y, chisq, nok, nrej, npar, nfev) where: `fit`
                is an :class:`Window` containing the best fit; `x` and
                `y` are the x and y positions of all pixels (2D numpy
                arrays); `chisq` is the raw chi**2; `nok` is the
                number of points fitted; `nrej` is the number
                rejected; `npar` is the number of parameters fitted (3
                to 5), `nfev` is the number of functions evaluations.

    Raises a HipercamError least_squares fails.

    """

    if fwhm < fwhm_min:
        raise HipercamError("fwhm out of range")

    # construct function object for gaussian fit. First the mode.
    mode = "" if sky is None else "s"
    mode = mode if fwhm_fix else mode + "f"
    gfit = Gfit(wind, sigma, ndiv, mode, fwhm)

    mode_switch = False

    first_fit = True
    nfev = 0
    while True:

        # retrieve current fit mode
        mode = gfit.mode

        # Rejection loop. There is a break statement at the end of the
        # loop if no pixels have been rejected.

        # set parameters
        param = gfit.set_par(sky, height, xcen, ycen, fwhm)

        # carry out fit
        res = least_squares(
            gfit.fun, param, jac=gfit.jac, method="lm", max_nfev=max_nfev
        )
        if not res.success:
            raise HipercamError(res.message)
        nfev += nfev

        # get Jacobian
        J = np.matrix(res.jac)
        try:
            # try / except for singularity
            covar = (J.T * J).I
        except (np.linalg.LinAlgError, ValueError) as err:
            raise HipercamError(err)

        param = res.x
        skyf, heightf, xf, yf, fwhmf = gfit.get_par(param)
        fwhmf = abs(fwhmf)

        if mode.find("f") > -1 and fwhmf < fwhm_min:
            # switch to fixed fwhm fit
            gfit.set_mode(mode.replace("f", ""), fwhm_min)

        else:

            # no mode switch, compute fit
            fit = Window(wind, gfit.model(param))

            covs = np.diag(covar)
            if (covs < 0).any():
                raise HipercamError("Negative covariance in fitGaussian")

            # compute chi**2 and number of OK points
            ok = gfit.mask & (sigma > 0)
            resid = (wind.data - fit.data) / sigma
            chisq = (resid[ok] ** 2).sum()
            nok1 = len(resid[ok])
            sfac = np.sqrt(chisq / nok1)

            # reject any above the defined threshold
            if first_fit:
                # first fit carried out with higher threshold for safety
                sigma[ok & (np.abs(resid) > 2*sfac*thresh)] *= -1
            else:
                sigma[ok & (np.abs(resid) > sfac*thresh)] *= -1

            # check whether any have been rejected
            ok = gfit.mask & (sigma > 0)
            nok = len(sigma[ok])

            if nok == nok1:
                # no more pixels have been rejected.  calculate how
                # many have been, re-scale the uncertainties to
                # reflect the actual chi**2
                nrej = len(sigma[gfit.mask]) - nok
                skyfe, heightfe, xfe, yfe, fwhmfe = gfit.get_epar(
                    sfac * np.sqrt(covs)
                )
                if not first_fit:
                    # always carry out at least 1
                    break
            first_fit = False

    # Make sure all regions masked off have sigma < 0
    sigma[~gfit.mask] = - np.abs(sigma[~gfit.mask])

    # OK we are done.
    extras = (
        fit,
        gfit.x,
        gfit.y,
        chisq,
        nok,
        nrej,
        len(param),
        nfev
    )

    if sky is None:
        return (
            (0., heightf, xf, yf, fwhmf),
            (-1., heightfe, xfe, yfe, fwhmfe),
            extras
        )
    else:
        return (
            (skyf, heightf, xf, yf, fwhmf),
            (skyfe, heightfe, xfe, yfe, fwhmfe),
            extras
        )

@jit(nopython=True, cache=True)
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

    alpha = 4.0 * np.log(2.0) / fwhm ** 2

    if ndiv > 0:
        # Complicated case with sub-pixellation allowed for

        # mean offset within sub-pixels
        soff = (ndiv - 1) / (2 * ndiv)

        prof = np.zeros_like(x)
        for iy in range(ybin):
            # loop over unbinned pixels in Y
            yoff = iy - (ybin - 1) / 2 - soff
            for ix in range(xbin):
                # loop over unbinned pixels in X
                xoff = ix - (xbin - 1) / 2 - soff
                for isy in range(ndiv):
                    # loop over sub-pixels in y
                    ysoff = yoff + isy / ndiv
                    for isx in range(ndiv):
                        # loop over sub-pixels in x
                        xsoff = xoff + isx / ndiv
                        rsq = (x + xsoff - xcen) ** 2 + (y + ysoff - ycen) ** 2
                        prof += np.exp(-alpha * rsq)

        return sky + (height / xbin / ybin / ndiv ** 2) * prof

    else:
        # Fast as possible, compute profile at pixel centres only
        rsq = (x - xcen) ** 2 + (y - ycen) ** 2
        return sky + height * np.exp(-alpha * rsq)


@jit(nopython=True, cache=True)
def dgaussian(x, y, sky, height, xcen, ycen, fwhm, xbin, ybin, ndiv, comp_dfwhm):
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

    comp_dfwhm : bool
         compute derivative wrt FWHM (or not). If False, derivative with respect
         to FWHM will not be computed

    Returns:: a list of five 2D numpy arrays containing the partial
    derivatives of a symmetric 2D Gaussian plus constant evaluated on
    the ordinate grids in xy. They come in the same order as they
    appear in the function call.

    """
    alpha = 4.0 * np.log(2.0) / fwhm ** 2

    dsky = np.ones_like(x)

    if ndiv > 0:
        # complicated sub-pixellation case
        dheight = np.zeros_like(x)
        dxcen = np.zeros_like(x)
        dycen = np.zeros_like(x)
        if comp_dfwhm:
            dfwhm = np.zeros_like(x)

        # mean offset within sub-pixels
        soff = (ndiv - 1) / (2 * ndiv)

        for iy in range(ybin):
            # loop over unbinned pixels in Y
            yoff = iy - (ybin - 1) / 2 - soff
            for ix in range(xbin):
                # loop over unbinned pixels in X
                xoff = ix - (xbin - 1) / 2 - soff
                for isy in range(ndiv):
                    # loop over sub-pixels in y
                    ysoff = yoff + isy / ndiv
                    for isx in range(ndiv):
                        # loop over sub-pixels in x
                        xsoff = xoff + isx / ndiv

                        # finally compute stuff
                        dx = x + xsoff - xcen
                        dy = y + ysoff - ycen
                        rsq = dx ** 2 + dy ** 2

                        dh = np.exp(-alpha * rsq)
                        dheight += dh
                        dxcen += (2 * alpha * height) * dh * dx
                        dycen += (2 * alpha * height) * dh * dy
                        if comp_dfwhm:
                            dfwhm += (2 * alpha * height / fwhm) * dh * rsq

        # Normalise by number of evaluations
        nadd = xbin * ybin * ndiv ** 2
        dheight /= nadd
        dxcen /= nadd
        dycen /= nadd
        if comp_dfwhm:
            dfwhm /= nadd
        #            return np.dstack((dsky, dheight, dxcen, dycen, dfwhm))
        if comp_dfwhm:
            return (dsky, dheight, dxcen, dycen, dfwhm)
        else:
            return (dsky, dheight, dxcen, dycen, dycen)

    else:
        # fast as possible, only compute at centre of pixels
        dx = x - xcen
        dy = y - ycen
        rsq = dx ** 2 + dy ** 2

        dheight = np.exp(-alpha * rsq)
        dxcen = (2 * alpha * height) * dheight * dx
        dycen = (2 * alpha * height) * dheight * dy
        if comp_dfwhm:
            dfwhm = (2 * alpha * height / fwhm) * dheight * rsq
            return (dsky, dheight, dxcen, dycen, dfwhm)
        else:
            return (dsky, dheight, dxcen, dycen, dycen)

class Gfit:
    """Object providing 'fun' and 'jac' methods for least_squares for
    Gaussian models. Four operating modes each of which allows the
    following to be free [else not]. Sky assumed = 0 when not free.

       mode == 'sf' : sky, FWHM
       mode == 's' : sky
       mode == 'f' : FWHM
       mode == '' : None of the above.

    It is up to the user to set the value of fwhm appropriate to the
    mode. A method 'set_mode' provides some basic checks.

    """

    def __init__(self, wind, sigma, ndiv, mode, fwhm=None):
        """
        Arguments::

          wind  : Window
             the Window containing data to fit

          sigma : np.array
             array of uncertainties matching wind

          ndiv : int
             pixel sub-division factor. See comments in fitMoffat

          mode : str
             see comments about operating mode in class description

          fwhm : float
             FWHM in unbinned pixels

        """
        self.sigma = sigma
        x = wind.x(np.arange(wind.nx))
        y = wind.y(np.arange(wind.ny))
        self.x, self.y = np.meshgrid(x, y)
        self.data = wind.data
        self.xbin = wind.xbin
        self.ybin = wind.ybin
        self.ndiv = ndiv
        self.mask = _mask(wind, self.x, self.y)
        self.set_mode(mode, fwhm)

    def set_mode(self, mode, fwhm):
        """Set the operation mode with some light checks"""
        if mode not in ("sf", "s", "f", ""):
            raise HipercamError("invalid mode = {}".format(mode))

        if mode.find("f") == -1 and fwhm is None:
            raise HipercamError("invalid mode / fwhm combination")
        self.mode = mode
        self.fwhm = fwhm

    def get_par(self, param):
        """Gets parameters (sky, height, xcen, ycen, fwhm) according to
        current least_squares parameter vector, accounting for the
        operating mode. All five are returned in each case, even when this
        is meaningless.

        """
        if self.mode == "sf":
            sky, height, xcen, ycen, fwhm = param
        elif self.mode == "s":
            sky, height, xcen, ycen = param
            fwhm = self.fwhm
        elif self.mode == "f":
            height, xcen, ycen, fwhm = param
            sky = 0.0
        elif self.mode == "":
            height, xcen, ycen = param
            sky = 0.0
            fwhm = self.fwhm
        else:
            raise HipercamError("invalid mode")

        return (sky, height, xcen, ycen, fwhm)

    def set_par(self, sky, height, xcen, ycen, fwhm):
        """Sets parameters (sky, height, xcen, ycen, fwhm) into correct
        vector for least_squares according to the operating mode. The
        vector is returned. e.g. if self.mode == "s", then (sky, height, xcen, ycen)
        will be returned.
        """
        if self.mode == "sf":
            return (sky, height, xcen, ycen, fwhm)
        elif self.mode == "s":
            return (sky, height, xcen, ycen)
        elif self.mode == "f":
            return (height, xcen, ycen, fwhm)
        elif self.mode == "":
            return (height, xcen, ycen)
        else:
            raise HipercamError("invalid mode")

    def get_epar(self, param):
        """Gets parameter errors (skye, heighte, xcene, ycene, fwhme)
        passed a vector of errors, accounting for the operating
        mode. All five are returned in each case, with fixed ones set =
        -1.
        """
        if self.mode == "sf":
            skye, heighte, xcene, ycene, fwhme = param
        elif self.mode == "s":
            skye, heighte, xcene, ycene = param
            fwhme = -1.0
        elif self.mode == "f":
            heighte, xcene, ycene, fwhme = param
            skye = -1.0
        elif self.mode == "":
            heighte, xcene, ycene = param
            skye = -1.0
            fwhme = -1.0
        else:
            raise HipercamError("invalid mode")

        return (skye, heighte, xcene, ycene, fwhme)

    def fun(self, param):
        """
        Returns 1D array of normalised residuals. See the model
        method for a description of the argument 'param'
        """
        mod = self.model(param)
        diff = (self.data - mod) / self.sigma
        ok = self.mask & (self.sigma > 0)
        return diff[ok].ravel()

    def jac(self, param):
        """
        Returns list of 1D arrays of the partial derivatives of
        the normalised residuals with respect to the variable
        parameters.
        """
        sky, height, xcen, ycen, fwhm = self.get_par(param)

        comp_fwhm = self.mode.find("f") > -1

        # work out which derivatives to bother with
        if self.mode == "sf":
            inds = (0, 1, 2, 3, 4)
        elif self.mode == "s":
            inds = (0, 1, 2, 3)
        elif self.mode == "f":
            inds = (1, 2, 3, 4)
        elif self.mode == "":
            inds = (1, 2, 3)
        else:
            raise HipercamError("invalid mode")

        derivs = dgaussian(
            self.x,
            self.y,
            sky,
            height,
            xcen,
            ycen,
            fwhm,
            self.xbin,
            self.ybin,
            self.ndiv,
            comp_fwhm,
        )

        ok = self.mask & (self.sigma > 0)
        return np.column_stack(
            [(-derivs[ind][ok] / self.sigma[ok]).ravel() for ind in inds]
        )

    def model(self, param):
        """
        Returns 2D array with model given a parameter vector.

        Argument::

           param : 1D array
              parameter vector, with values that depend upon the mode,
              but could include some or all of (sky, height, xcen, ycen,
              fwhm) where 'sky' is the background per pixel;
              'height' is the central height of the gaussian function;
              'xcen' and 'ycen' are the ordinates of its centre in
              unbinned CCD pixels with (1,1) at the left corner of the
              physical imaging area; 'fwhm' is the FWHM in unbinned
              pixels.
        """
        sky, height, xcen, ycen, fwhm = self.get_par(param)
        return gaussian(
            self.x,
            self.y,
            sky,
            height,
            xcen,
            ycen,
            fwhm,
            self.xbin,
            self.ybin,
            self.ndiv,
        )

