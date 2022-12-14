import sys
import warnings

import numpy as np
from astropy import units as u
from astropy.modeling import Fittable2DModel, Parameter
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import SigmaClip, gaussian_fwhm_to_sigma, gaussian_sigma_to_fwhm
from astropy.table import Table
from photutils.background import MMMBackground
from photutils.psf import DAOGroup, BasicPSFPhotometry
from scipy.special import erf

import hipercam as hcam
from hipercam.reduction import moveApers


# Stuff below here are helper routines that are not exported
class IntegratedGaussianPRF2(Fittable2DModel):
    r"""
    Circular Gaussian model integrated over pixels.

    Because it is integrated, this model is considered a PRF, *not* a
    PSF (see :ref:`psf-terminology` for more about the terminology used
    here.)

    This model is a Gaussian *integrated* over an area of
    ``1`` (in units of the model input coordinates, e.g., 1
    pixel). This is in contrast to the apparently similar
    `astropy.modeling.functional_models.Gaussian2D`, which is the value
    of a 2D Gaussian *at* the input coordinates, with no integration.
    So this model is equivalent to assuming the PSF is Gaussian at a
    *sub-pixel* level.

    Based on the same named model in photutils, but with different
    std deviations in x and y.

    Parameters
    ----------
    sigma_x : float
        Width of the Gaussian PSF in x.
    sigma_y : float
        Width of the Gaussian PSF in y.
    flux : float, optional
        Total integrated flux over the entire PSF
    x_0 : float, optional
        Position of the peak in x direction.
    y_0 : float, optional
        Position of the peak in y direction.

    Notes
    -----
    This model is evaluated according to the following formula:

        .. math::

            f(x, y) =
                \frac{F}{4}
                \left[
                {\rm erf} \left(\frac{x - x_0 + 0.5}
                {\sqrt{2} \sigma_x} \right) -
                {\rm erf} \left(\frac{x - x_0 - 0.5}
                {\sqrt{2} \sigma_x} \right)
                \right]
                \left[
                {\rm erf} \left(\frac{y - y_0 + 0.5}
                {\sqrt{2} \sigma_y} \right) -
                {\rm erf} \left(\frac{y - y_0 - 0.5}
                {\sqrt{2} \sigma_y} \right)
                \right]

    where ``erf`` denotes the error function and ``F`` the total
    integrated flux.
    """

    flux = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    sigma_x = Parameter(default=1, fixed=True)
    sigma_y = Parameter(default=1, fixed=True)

    _erf = None

    @property
    def bounding_box(self):
        halfwidth = 2 * (self.sigma_x + self.sigma_y)
        return (
            (int(self.y_0 - halfwidth), int(self.y_0 + halfwidth)),
            (int(self.x_0 - halfwidth), int(self.x_0 + halfwidth)),
        )

    @property
    def fwhm(self):
        return u.Quantity(0.5 * (self.sigma_x + self.sigma_y)) * gaussian_sigma_to_fwhm

    def evaluate(self, x, y, flux, x_0, y_0, sigma_x, sigma_y):
        """Model function Gaussian PSF model."""
        return (
            flux
            / 4
            * (
                (
                    erf((x - x_0 + 0.5) / (np.sqrt(2) * sigma_x))
                    - erf((x - x_0 - 0.5) / (np.sqrt(2) * sigma_x))
                )
                * (
                    erf((y - y_0 + 0.5) / (np.sqrt(2) * sigma_y))
                    - erf((y - y_0 - 0.5) / (np.sqrt(2) * sigma_y))
                )
            )
        )


class MoffatPSF(Fittable2DModel):
    """
    Standard Moffat model, but with parameter names that work with photutils

    Different FWHM parameters apply in x and y, the fwhm property gives
    the average of both.

    Parameters
    ----------
    fwhm_x : float
        Full-Width at Half-Maximum of profile in x
    fwhm_y : float
        Full-Width at Half-Maximum of profile in y
    beta : float
        Beta parameter of Moffat profile
    flux : float (default 1)
        Total integrated flux over the entire PSF
    x_0 : float (default 0)
        Position of peak in x direction
    y_0 : float (default 0)
        Position of peak in y direction

    Notes
    -----
    This model is evaluated according to the following formula:

        .. math::
            f(x, y) = F [1 + ((x-x_o)/\alpha_x)**2 + ((y-y_o)/\alpha_y)**2]**(-beta)

    where ``F`` is the total integrated flux, and ``alpha_i`` is chosen
    so the profile has the appropriate FWHM.
    """

    flux = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    beta = Parameter(default=2.5, fixed=True)
    fwhm_x = Parameter(default=12, fixed=True)
    fwhm_y = Parameter(default=12, fixed=True)

    fit_deriv = None

    @property
    def fwhm(self):
        return u.Quantity(0.5 * (self.fwhm_x + self.fwhm_y))

    @property
    def bounding_box(self):
        halfwidth = 3 * self.fwhm
        return (
            (int(self.y_0 - halfwidth), int(self.y_0 + halfwidth)),
            (int(self.x_0 - halfwidth), int(self.x_0 + halfwidth)),
        )

    def evaluate(self, x, y, flux, x_0, y_0, beta, fwhm_x, fwhm_y):
        alpha_sq_x = fwhm_x**2 / 4 / (2 ** (1 / beta) - 1)
        alpha_sq_y = fwhm_y**2 / 4 / (2 ** (1 / beta) - 1)
        x_term = (x - x_0) ** 2 / alpha_sq_x
        y_term = (y - y_0) ** 2 / alpha_sq_y
        prof = flux / (1 + x_term + y_term) ** beta
        return (beta - 1) * prof / np.pi / np.sqrt(alpha_sq_x * alpha_sq_y)


def ccdproc_psf(cnam, ccds, bccds, rccds, nframes, read, gain, ccdwin, rfile, store):
    """Processing steps for a sequential set of images from the same
    CCD. This is designed for parallelising the processing across CCDs
    of multiarm cameras like ULTRACAM and HiPERCAM using
    multiprocessing. To be called *after* checking that any processing
    is needed.

    Arguments::

       cnam : string
          name of CCD, for information purposes (e.g. 'red', '3', etc)

       ccds : List of CCDs
          the CCDs for processing which should have been debiassed, flat
          fielded and multiplied by the gain to get into electrons.

       bccds : List of CCDs
          the CCDs for processing which should just have been debiassed

       rccds : List of CCDs
          unprocessed CCDs, one-to-one correspondence with 'ccds', used
          to measure saturation

       nframes : List of ints
          frame numbers for each CCD

       read : CCD
          readnoise frame divided by the flatfield

       gain : CCD
          gain frame dmultiplied by the flatfield

       ccdwin : dict
          label of the Window enclosing each aperture

       rfile : Rfile object
          reduction control parameters. rfile.aper used to store the aperture
          parameters.

       store : dict
          dictionary of results

    Returns: (cnam, list[res]) where 'res' represents a tuple
    of results for each input CCD and contains the following:

    (nframe, store, ccdaper, results, mjdint, mjdfrac, mjdok, expose)

    """
    # At this point 'ccds' contains a list of CCD each of which
    # contains all the Windows of a CCD, 'ccdaper' all of its
    # apertures, 'ccdwin' the label of the Window enclosing each
    # aperture, 'rfile' contains control parameters, 'rflat' contains
    # the readout noise in electrons and divided by the flat as a CCD,
    # 'store' is a dictionary initially with jus 'mfwhm' and 'mbeta'
    # set = -1, but will pick up extra stuff from moveApers for use by
    # extractFlux along with revised values of mfwhm and mbeta which
    # are used to initialise profile fits next time.

    res = []
    for ccd, bccd, rccd, nframe in zip(ccds, bccds, rccds, nframes):
        # Loop through the CCDs supplied

        # move the apertures.
        moveApers(cnam, ccd, bccd, read, gain, ccdwin, rfile, store)

        # extract flux from all apertures of each CCD. Return with the CCD
        # name, the store dictionary, ccdaper and then the results from
        # extractFlux for compatibility with multiprocessing. Note
        results = extractFluxPSF(
            cnam, ccd, bccd, rccd, read, gain, ccdwin, rfile, store
        )

        # Save the essentials
        res.append(
            (
                nframe,
                store,
                rfile.aper[cnam],
                results,
                ccd.head["MJDINT"],
                ccd.head["MJDFRAC"],
                ccd.head.get("GOODTIME", True),
                ccd.head.get("EXPTIME", 1.0),
            )
        )

    return (cnam, res)


def extractFluxPSF(cnam, ccd, bccd, rccd, read, gain, ccdwin, rfile, store):
    """This extracts the flux of all apertures of a given CCD.

    The steps are (1) creation of PSF model, (2) PSF fitting, (3)
    flux extraction. The apertures are assumed to be correctly positioned.

    It returns the results as a dictionary keyed on the aperture label. Each
    entry returns a list:

    [x, ex, y, ey, fwhm, efwhm, beta, ebeta, counts, countse, sky, esky,
    nsky, nrej, flag, cmax]

    flag = bitmask. See hipercam.core to see all the options which are
    referred to by name in the code e.g. ALL_OK. The various flags can
    signal that there no sky pixels (NO_SKY), the sky aperture was off
    the edge of the window (SKY_AT_EDGE), etc.

    This code::

       >> bset = flag & TARGET_SATURATED

    determines whether the data saturation flag is set for example.

    .. note::
        not all arguments are used by this routine, but the call signature is
        fixed to be a drop-in replacement for the normal and optimal photometry
        function called by reduce

    Arguments::

       cnam : string
          CCD identifier label

       ccd : CCD
           the fully processed CCD (debiassed, dark-corrected, flat-fielded, fringe-corrected)

       bccd : CCD
           debiassed-only CCD. Used in variance computation.

       rccd : CCD
          corresponding raw CCD, used to work out whether data are
          saturated in target aperture.

       read : CCD
           readnoise divided by the flat-field

       gain : CCD
           gain multiplied by the flat field

       ccdwin : dictionary of strings
           the Window label corresponding to each Aperture

       rfile : Rfile
           reduce file configuration parameters

       store : dict of dicts
           see moveApers for what this contains.
    """

    # initialise flag
    flag = hcam.ALL_OK

    ccdaper = rfile.aper[cnam]

    results = {}

    # get profile params from aperture store
    mfwhm = store["mfwhm"]
    mbeta = store["mbeta"]
    method = "m" if mbeta > 0.0 else "g"

    if mfwhm <= 0:
        # die hard, die soon as there's nothing we can do.
        print(
            (
                " *** WARNING: CCD {:s}: no measured FWHM to create PSF model"
                "; no extraction possible"
            ).format(cnam)
        )
        # set flag to indicate no FWHM
        flag = hcam.NO_FWHM

        for apnam, aper in ccdaper.items():
            info = store[apnam]
            results[apnam] = {
                "x": aper.x,
                "xe": info["xe"],
                "y": aper.y,
                "ye": info["ye"],
                "fwhm": info["fwhm"],
                "fwhme": info["fwhme"],
                "beta": info["beta"],
                "betae": info["betae"],
                "counts": 0.0,
                "countse": -1,
                "sky": 0.0,
                "skye": 0.0,
                "nsky": 0,
                "nrej": 0,
                "flag": flag,
                "cmax": 0,
            }
        return results

    # all apertures have to be in the same window, or we can't easily make a
    # postage stamp of the data
    wnames = set(ccdwin.values())
    if len(wnames) != 1:
        print(
            (
                " *** WARNING: CCD {:s}: not all apertures"
                " lie within the same window; no extraction possible"
            ).format(cnam)
        )

        # set flag to indicate no extraction
        flag = hcam.NO_EXTRACTION

        # return empty results
        for apnam, aper in ccdaper.items():
            info = store[apnam]
            results[apnam] = {
                "x": aper.x,
                "xe": info["xe"],
                "y": aper.y,
                "ye": info["ye"],
                "fwhm": info["fwhm"],
                "fwhme": info["fwhme"],
                "beta": info["beta"],
                "betae": info["betae"],
                "counts": 0.0,
                "countse": -1,
                "sky": 0.0,
                "skye": 0.0,
                "nsky": 0,
                "nrej": 0,
                "flag": flag,
                "cmax": 0,
            }
            return results
    wnam = wnames.pop()

    # PSF params are in binned pixels, so find binning
    bin_fac = ccd[wnam].xbin

    # create PSF model
    if method == "m":
        psf_model = MoffatPSF(
            beta=mbeta, fwhm_x=mfwhm / bin_fac, fwhm_y=mfwhm / bin_fac
        )
        psf_model.fwhm_x.min = 0.5 * mfwhm / bin_fac
        psf_model.fwhm_x.max = 1.5 * mfwhm / bin_fac
        psf_model.fwhm_y.min = 0.5 * mfwhm / bin_fac
        psf_model.fwhm_y.max = 1.5 * mfwhm / bin_fac
        psf_model.fwhm_x.fixed = True
        psf_model.fwhm_y.fixed = True
        extra_output_cols = ["fwhm_x", "fwhm_y"]
    else:
        psf_model = IntegratedGaussianPRF2(
            sigma_x=mfwhm * gaussian_fwhm_to_sigma / bin_fac,
            sigma_y=mfwhm * gaussian_fwhm_to_sigma / bin_fac,
        )
        extra_output_cols = ["sigma_x", "sigma_y"]

    # force photometry only at aperture positions
    # this means PSF shape and positions are fixed, we are only fitting flux
    if rfile["psf_photom"]["positions"] == "fixed":
        psf_model.x_0.fixed = True
        psf_model.y_0.fixed = True

    # create instances for PSF photometry
    gfac = float(rfile["psf_photom"]["gfac"])
    sclip = float(rfile["sky"]["thresh"])
    daogroup = DAOGroup(gfac * mfwhm / bin_fac)
    mmm_bkg = MMMBackground(sigma_clip=SigmaClip(sclip))
    fitter = LevMarLSQFitter()
    fitshape_box_size = int(2 * int(float(rfile["psf_photom"]["fit_half_width"])) + 1)
    fitshape = (fitshape_box_size, fitshape_box_size)

    photometry_task = BasicPSFPhotometry(
        group_maker=daogroup,
        bkg_estimator=mmm_bkg,
        psf_model=psf_model,
        fitter=fitter,
        fitshape=fitshape,
        extra_output_cols=extra_output_cols,
    )

    # initialise flag
    flag = hcam.ALL_OK

    # extract Windows relevant for these apertures
    wdata = ccd[wnam]
    wraw = rccd[wnam]

    # extract sub-windows that include all of the apertures, plus a little
    # extra around the edges.
    x1 = min([ap.x - ap.rsky2 - wdata.xbin for ap in ccdaper.values()])
    x2 = max([ap.x + ap.rsky2 + wdata.xbin for ap in ccdaper.values()])
    y1 = min([ap.y - ap.rsky2 - wdata.ybin for ap in ccdaper.values()])
    y2 = max([ap.y + ap.rsky2 + wdata.ybin for ap in ccdaper.values()])

    # extract sub-Windows
    swdata = wdata.window(x1, x2, y1, y2)
    swraw = wraw.window(x1, x2, y1, y2)

    # compute pixel positions of apertures in windows
    xpos, ypos = zip(
        *((swdata.x_pixel(ap.x), swdata.y_pixel(ap.y)) for ap in ccdaper.values())
    )
    positions = Table(names=["x_0", "y_0"], data=(xpos, ypos))

    # do the PSF photometry
    photom_results = photometry_task(swdata.data, init_guesses=positions)
    slevel = mmm_bkg(swdata.data)

    # unpack the results and check apertures
    for apnam, aper in ccdaper.items():
        try:
            # reset flag
            flag = hcam.ALL_OK

            result_row = photom_results[photom_results["id"] == int(apnam)]
            if len(result_row) == 0:
                flag |= hcam.NO_DATA
                raise hcam.HipercamError(
                    "no source in PSF photometry for this aperture"
                )
            elif len(result_row) > 1:
                flag |= hcam.NO_EXTRACTION
                raise hcam.HipercamError(
                    "ambiguous lookup for this aperture in PSF photometry"
                )
            else:
                result_row = result_row[0]

            # compute X, Y arrays over the sub-window relative to the centre
            # of the aperture and the distance squared from the centre (Rsq)
            # to save a little effort.
            x = swdata.x(np.arange(swdata.nx)) - aper.x
            y = swdata.y(np.arange(swdata.ny)) - aper.y
            X, Y = np.meshgrid(x, y)
            Rsq = X**2 + Y**2

            # size of a pixel which is used to taper pixels as they approach
            # the edge of the aperture to reduce pixellation noise
            size = np.sqrt(wdata.xbin * wdata.ybin)

            # target selection, accounting for extra apertures and allowing
            # pixels to contribute if their centres are as far as size/2 beyond
            # the edge of the circle (but with a tapered weight)
            dok = Rsq < (aper.rtarg + size / 2.0) ** 2
            if not dok.any():
                # check there are some valid pixels
                flag |= hcam.NO_DATA
                raise hcam.HipercamError("no valid pixels in aperture")

            # check for saturation and nonlinearity
            cmax = int(swraw.data[dok].max())
            if cnam in rfile.warn:
                if cmax >= rfile.warn[cnam]["saturation"]:
                    flag |= hcam.TARGET_SATURATED

                if cmax >= rfile.warn[cnam]["nonlinear"]:
                    flag |= hcam.TARGET_NONLINEAR
            else:
                warnings.warn("CCD {:s} has no nonlinearity or saturation levels set")

            counts = result_row["flux_fit"]
            try:
                countse = result_row["flux_unc"]
            except KeyError:
                raise hcam.HipercamError(
                    "unable to find errors on solution, model does not depend on params"
                )
            info = store[apnam]

            results[apnam] = {
                "x": aper.x,
                "xe": info["xe"],
                "y": aper.y,
                "ye": info["ye"],
                "fwhm": info["fwhm"],
                "fwhme": info["fwhme"],
                "beta": info["beta"],
                "betae": info["betae"],
                "counts": counts,
                "countse": countse,
                "sky": slevel,
                "skye": 0,
                "nsky": 0,
                "nrej": 0,
                "flag": flag,
                "cmax": cmax,
            }

        except hcam.HipercamError as err:
            print(
                f"CCD {cnam}, reference aperture {apnam},"
                f" fit failed; extraction aborted. Error = {err}",
                file=sys.stderr,
            )
            info = store[apnam]
            flag |= hcam.NO_EXTRACTION

            results[apnam] = {
                "x": aper.x,
                "xe": info["xe"],
                "y": aper.y,
                "ye": info["ye"],
                "fwhm": info["fwhm"],
                "fwhme": info["fwhme"],
                "beta": info["beta"],
                "betae": info["betae"],
                "counts": 0.0,
                "countse": -1,
                "sky": 0.0,
                "skye": 0.0,
                "nsky": 0,
                "nrej": 0,
                "flag": flag,
                "cmax": 0,
            }

    # finally, we are done
    return results
