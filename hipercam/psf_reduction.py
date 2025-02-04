import sys
import warnings

import numpy as np
from astropy import units as u
from astropy.modeling import Fittable2DModel, Parameter
from astropy.modeling.utils import ellipse_extent
from astropy.stats import SigmaClip
from astropy.table import Table
from matplotlib import pyplot as plt
from photutils.background import LocalBackground, MedianBackground
from photutils.psf import GaussianPRF, PSFPhotometry
from photutils.psf.functional_models import FLOAT_EPSILON, GAUSSIAN_FWHM_TO_SIGMA

import hipercam as hcam

# Stuff below here are helper routines that are not exported


class MoffatPSF(Fittable2DModel):
    r"""
    Standard Moffat model, but with parameter names that work with photutils

    Different FWHM parameters apply in x and y, the fwhm property gives
    the average of both.

    Parameters
    ----------
    flux : float, optional
        Total integrated flux over the entire PSF.

    x_0 : float, optional
        Position of the peak along the x axis.

    y_0 : float, optional
        Position of the peak along the y axis.

    x_fwhm : float, optional
        Full-Width at Half-Maximum of profile in x

    y_fwhm : float, optional
        Full-Width at Half-Maximum of profile in y

    beta : float, optional
        Beta parameter of Moffat profile

    theta : float, optional
        The counterclockwise rotation angle either as a float (in
        degrees) or a `~astropy.units.Quantity` angle (optional).

    bbox_factor : float, optional
        The multiple of the FWHM used to define the bounding box limits.

    Notes
    -----
    This model is evaluated according to the following formula:

        .. math::
            f(x, y) =  F \frac{beta - 1}{\pi \alpha^2}
                \left[1 + a(x-x_0)**2 + +b(x-x_0)(y-y_0) + c(y-y_0)**2 \right]**(-beta)

    where :math:`F` is the total integrated flux, :math:`(x_{0},
    y_{0})` is the position of the peak. Note that :math:`beta` must be
    greater than 1.

    The parameters :math:`a`, :math:`b` and :math:`c` are given by:

        .. math::
            a = \cos(\theta)^2 / \alpha_x^2 + \sin(\theta)^2 / \alpha_y^2
            b = \sin(2\theta) / \alpha_x^2 - \sin(2\theta) / \alpha_y^2
            c = \sin(\theta)^2 / \alpha_x^2 + \cos(\theta)^2 / \alpha_y^2

    where :math:`\alpha_x` and :math:`\alpha_y` are chosen to given the correct
    FWHM in the x and y directions respectively. :math:`\theta` is the anti-clockwise
    rotation of the PSF. :math:`alpha` is the geometric mean of :math:`\alpha_x` and
    :math:`\alpha_y`. The FWHM of the Moffat profile is given by:

    .. math::

        \rm{FWHM} = 2 \alpha \sqrt{2^{1 / beta} - 1}

    The model is normalized such that, for :math:`beta > 1`:

    .. math::

        \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} f(x, y)
            \,dx \,dy = F
    """

    flux = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    x_fwhm = Parameter(default=12, bounds=(FLOAT_EPSILON, None), fixed=False)
    y_fwhm = Parameter(default=12, bounds=(FLOAT_EPSILON, None), fixed=False)
    beta = Parameter(default=2.5, bounds=(1.0, None), fixed=False)
    theta = Parameter(default=0, bounds=(-90, 90), fixed=False)

    def __init__(
        self,
        *,
        flux=flux.default,
        x_0=x_0.default,
        y_0=y_0.default,
        x_fwhm=x_fwhm.default,
        y_fwhm=y_fwhm.default,
        beta=beta.default,
        theta=theta.default,
        bbox_factor=15.5,
        **kwargs,
    ):
        super().__init__(
            flux=flux,
            x_0=x_0,
            y_0=y_0,
            x_fwhm=x_fwhm,
            y_fwhm=y_fwhm,
            beta=beta,
            theta=theta,
            **kwargs,
        )
        self.bbox_factor = bbox_factor

    @property
    def fwhm(self):
        return u.Quantity(0.5 * (self.x_fwhm + self.y_fwhm))

    def _calc_bounding_box(self, factor=5.5):
        """
        Calculate a bounding box defining the limits of the model.

        The limits are adjusted for rotation.

        Parameters
        ----------
        factor : float, optional
            The multiple of the x and y standard deviations (sigma) used
            to define the limits.

        Returns
        -------
        bbox : tuple
            A bounding box defining the ((y_min, y_max), (x_min, x_max))
            limits of the model.
        """
        a = factor * self.x_fwhm * GAUSSIAN_FWHM_TO_SIGMA
        b = factor * self.y_fwhm * GAUSSIAN_FWHM_TO_SIGMA
        theta = self.theta
        if not isinstance(theta, u.Quantity):
            theta = np.deg2rad(theta)

        dx, dy = ellipse_extent(a, b, theta)
        return ((self.y_0 - dy, self.y_0 + dy), (self.x_0 - dx, self.x_0 + dx))

    @property
    def bounding_box(self):
        """
        The bounding box of the model.

        Examples
        --------
        >>> from hipercam.psf_reduction import MoffatPSF
        >>> model = MoffatPSF(x_0=0, y_0=0, x_fwhm=5, y_fwhm=3, theta=20)
        >>> model.bounding_box  # doctest: +FLOAT_CMP
        ModelBoundingBox(
            intervals={
                x: Interval(lower=-11.232523683597986, upper=11.232523683597986)
                y: Interval(lower=-7.701096862895421, upper=7.701096862895421)
            }
            model=MoffatPSF(inputs=('x', 'y'))
            order='C'
        )
        >>> model.bbox_factor = 7
        >>> model.bounding_box  # doctest: +FLOAT_CMP
        ModelBoundingBox(
            intervals={
                x: Interval(lower=-14.295939233670163, upper=14.295939233670163)
                y: Interval(lower=-9.801396007321445, upper=9.801396007321445)
            }
            model=MoffatPSF(inputs=('x', 'y'))
            order='C'
        )
        """
        return self._calc_bounding_box(factor=self.bbox_factor)

    def evaluate(self, x, y, flux, x_0, y_0, x_fwhm, y_fwhm, beta, theta):
        """
        Calculate the value of the 2D Moffat model at the input
        coordinates for the given model parameters.

        Parameters
        ----------
        x, y : float or array_like
            The x and y coordinates at which to evaluate the model.

        flux : float
            Total integrated flux over the entire PSF.

        x_0, y_0 : float
            Position of the peak along the x and y axes.

        x_fwhm, y_fwhm : float
            FWHM of the Gaussian along the x and y axes.

        beta : float, optional
            Beta parameter of Moffat profile

        theta : float
            The counterclockwise rotation angle either as a float (in
            degrees) or a `~astropy.units.Quantity` angle (optional).

        Returns
        -------
        result : `~numpy.ndarray`
            The value of the model evaluated at the input coordinates.
        """
        # find useful derived properties of theta
        if not isinstance(theta, u.Quantity):
            theta = np.deg2rad(theta)
        cost2 = np.cos(theta) ** 2
        sint2 = np.sin(theta) ** 2
        sin2t = np.sin(2.0 * theta)

        # alpha_sq terms for x and y
        alpha_sq_x = x_fwhm**2 / 4 / (2 ** (1 / beta) - 1)
        alpha_sq_y = y_fwhm**2 / 4 / (2 ** (1 / beta) - 1)

        # terms a, b, c
        a = cost2 / alpha_sq_x + sint2 / alpha_sq_y
        b = sin2t / alpha_sq_x - sin2t / alpha_sq_y
        c = sint2 / alpha_sq_x + cost2 / alpha_sq_y
        x_term = a * (x - x_0) ** 2
        xy_term = b * (x - x_0) * (y - y_0)
        y_term = c * (y - y_0) ** 2
        profile = (1 + x_term + y_term + xy_term) ** (-beta)
        normalisation = (beta - 1) / np.pi / np.sqrt(alpha_sq_x * alpha_sq_y)
        return normalisation * flux * profile

    @staticmethod
    def fit_deriv(x, y, flux, x_0, y_0, x_fwhm, y_fwhm, beta, theta):
        """
        Calculate the partial derivatives of the 2D Moffat function
        with respect to the parameters.

        Parameters
        ----------
        x, y : float or array_like
            The x and y coordinates at which to evaluate the model.

        flux : float
            Total integrated flux over the entire PSF.

        x_0, y_0 : float
            Position of the peak along the x and y axes.

        x_fwhm, y_fwhm : float
            FWHM of the Gaussian along the x and y axes.

        beta : float, optional
            Beta parameter of Moffat profile

        theta : float
            The counterclockwise rotation angle either as a float (in
            degrees) or a `~astropy.units.Quantity` angle (optional).

        Returns
        -------
        result : list of `~numpy.ndarray`
            The list of partial derivatives with respect to each
            parameter.
        """
        if not isinstance(theta, u.Quantity):
            theta = np.deg2rad(theta)
        else:
            theta = theta.to_value(u.rad)

        alpha_sq_x = x_fwhm**2 / 4 / (2 ** (1 / beta) - 1)
        alpha_sq_y = y_fwhm**2 / 4 / (2 ** (1 / beta) - 1)
        alpha_x = np.sqrt(alpha_sq_x)
        alpha_y = np.sqrt(alpha_sq_y)

        x1 = beta - 1
        x2 = 1 / np.pi
        x3 = 1 / alpha_x
        x4 = 1 / alpha_y
        x5 = x - x_0
        x6 = x5**2
        x7 = alpha_x ** (-2)
        x8 = np.cos(theta)
        x9 = x8**2
        x10 = alpha_y ** (-2)
        x11 = np.sin(theta)
        x12 = x11**2
        x13 = x10 * x12 + x7 * x9
        x14 = y - y_0
        x15 = x14**2
        x16 = x10 * x9 + x12 * x7
        x17 = 2 * theta
        x18 = np.sin(x17)
        x19 = -x10 * x18 + x18 * x7
        x20 = x14 * x19
        x21 = x13 * x6 + x15 * x16 + x20 * x5 + 1
        x22 = x21 ** (-beta)
        x23 = x1 * x2 * x22 * x3 * x4
        x24 = 1 / x21
        x25 = 2 / alpha_x**3
        x26 = x14 * x5
        x27 = alpha_y ** (-3)
        x28 = 2 * x27
        x29 = 2 * x11 * x8
        x30 = x10 * x29 - x29 * x7
        x31 = np.cos(x17)

        dg_dflux = x23
        dg_dx_0 = -flux * beta * x23 * x24 * (x13 * (-2 * x + 2 * x_0) - x20)
        dg_dy_0 = -flux * beta * x23 * x24 * (x16 * (-2 * y + 2 * y_0) - x19 * x5)
        dg_dalpha_x = (
            -flux
            * beta
            * x23
            * x24
            * (-x12 * x15 * x25 - x18 * x25 * x26 - x25 * x6 * x9)
            - flux * x1 * x2 * x22 * x4 * x7
        )
        dg_dalpha_y = (
            -flux
            * beta
            * x23
            * x24
            * (-x12 * x28 * x6 + 2 * x14 * x18 * x27 * x5 - x15 * x28 * x9)
            - flux * x1 * x10 * x2 * x22 * x3
        )
        partial_dg_dbeta = flux * x2 * x22 * x3 * x4 - flux * x23 * np.log(x21)
        dg_dtheta = (
            (
                -flux
                * beta
                * x23
                * x24
                * (-x15 * x30 + x26 * (-2 * x10 * x31 + 2 * x31 * x7) + x30 * x6)
            )
            * np.pi
            / 180
        )

        # jacobians
        dg_dx_fwhm = alpha_x * dg_dalpha_x / x_fwhm
        dg_dy_fwhm = alpha_y * dg_dalpha_y / y_fwhm
        dalpha_x_dbeta = (
            2 ** (1 / beta) * np.sqrt(x_fwhm**2 / (2 ** (1 / beta) - 1)) * np.log(2)
        )
        dalpha_x_dbeta /= 4 * beta**2 * (2 ** (1 / beta) - 1)
        dalpha_y_dbeta = (
            2 ** (1 / beta) * np.sqrt(y_fwhm**2 / (2 ** (1 / beta) - 1)) * np.log(2)
        )
        dalpha_y_dbeta /= 4 * beta**2 * (2 ** (1 / beta) - 1)
        dg_dbeta = (
            partial_dg_dbeta
            + dg_dalpha_x * dalpha_x_dbeta
            + dg_dalpha_y * dalpha_y_dbeta
        )

        return dg_dflux, dg_dx_0, dg_dy_0, dg_dx_fwhm, dg_dy_fwhm, dg_dbeta, dg_dtheta


def create_psf_model(photom_results, method, fixed_positions=False):
    if method == "moffat":
        psf_model = MoffatPSF(flux=1)
        for param in ["x_fwhm", "y_fwhm", "theta", "beta"]:
            colname = param + "_fit"
            if colname in photom_results.colnames:
                setattr(psf_model, param, np.median(photom_results[colname]))
        psf_model.beta.fixed = True
    else:
        psf_model = GaussianPRF(flux=1)

    if fixed_positions:
        psf_model.x_0.fixed = True
        psf_model.y_0.fixed = True

    psf_model.x_fwhm.fixed = True
    psf_model.y_fwhm.fixed = True
    psf_model.theta.fixed = True
    return psf_model


def display(data, phot, psfphot):
    from astropy.visualization import simple_norm

    colnames = [
        "id",
        "group_id",
        "group_size",
        "local_bkg",
        "x_fit",
        "y_fit",
        "flux_fit",
        "x_err",
        "y_err",
        "flux_err",
        "x_fwhm_fit",
        "y_fwhm_fit",
        "theta_fit",
        "beta_fit",
        "npixfit",
        "qfit",
        "cfit",
        "flags",
    ]

    display_colnames = []
    for colname in colnames:
        if colname in phot.colnames:
            if colname not in ["id", "group_id", "group_size", "flags"]:
                phot[colname].info.format = ".4f"
            display_colnames.append(colname)

    print(phot[display_colnames])
    resid = psfphot.make_residual_image(data)
    norm = simple_norm(data, "sqrt", percent=95)
    bkg = np.median(phot["local_bkg"])
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
    plt.tight_layout()
    ax[0].imshow(data, origin="lower", norm=norm)
    ax[1].imshow(data - resid + bkg, origin="lower", norm=norm)
    ax[2].imshow(resid, origin="lower", norm=norm)
    ax[0].set_title("Data")
    ax[1].set_title("Model")
    ax[2].set_title("Residual Image")

    plt.show()


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

    # we need at least one reference aperture to proceed
    nref = sum(1 for ap in ccdaper.values() if ap.ref)
    if nref == 0:
        print(
            (
                " *** WARNING: CCD {:s}: no reference apertures"
                " defined within window; no extraction possible"
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
    method = rfile["psf_photom"]["psf_model"]
    if method == "moffat":
        # fix beta if initial aperture tweak used Gaussian profiles
        beta = mbeta if mbeta > 0.0 else 2.5
        psf_model = MoffatPSF(beta=beta, x_fwhm=mfwhm / bin_fac, y_fwhm=mfwhm / bin_fac)
        # let the shape vary when we constrain psf_model based on reference stars
        psf_model.x_fwhm.fixed = False
        psf_model.y_fwhm.fixed = False
        psf_model.beta.fixed = False
        psf_model.theta.fixed = False
        psf_model.theta.bounds = (-90, 90)
    else:
        psf_model = GaussianPRF(
            x_fwhm=mfwhm / bin_fac,
            y_fwhm=mfwhm / bin_fac,
            theta=0.0,
        )
        # let the shape vary when we constrain psf_model based on reference stars
        psf_model.x_fwhm.fixed = False
        psf_model.y_fwhm.fixed = False
        psf_model.theta.fixed = False
        psf_model.theta.bounds = (-90, 90)

    # force photometry only at aperture positions
    # this means PSF shape and positions are fixed, we are only fitting flux
    if rfile["psf_photom"]["positions"] == "fixed":
        psf_model.x_0.fixed = True
        psf_model.y_0.fixed = True

    # initialise flag
    flag = hcam.ALL_OK

    # extract Windows relevant for these apertures
    wdata = ccd[wnam]
    wraw = rccd[wnam]
    wbias = bccd[wnam]
    wread = read[wnam]
    wgain = gain[wnam]

    # extract sub-windows that include all of the apertures, plus a little
    # extra around the edges.
    x1 = min([ap.x - ap.rsky2 - wdata.xbin for ap in ccdaper.values()])
    x2 = max([ap.x + ap.rsky2 + wdata.xbin for ap in ccdaper.values()])
    y1 = min([ap.y - ap.rsky2 - wdata.ybin for ap in ccdaper.values()])
    y2 = max([ap.y + ap.rsky2 + wdata.ybin for ap in ccdaper.values()])

    # extract sub-Windows
    swdata = wdata.window(x1, x2, y1, y2)
    swraw = wraw.window(x1, x2, y1, y2)
    swbias = wbias.window(x1, x2, y1, y2)
    swread = wread.window(x1, x2, y1, y2)
    swgain = wgain.window(x1, x2, y1, y2)
    # This is where the pure-debiassed data are used
    sigma = np.sqrt(swread.data**2 + np.maximum(0, swbias.data) / swgain.data)

    # compute pixel positions of apertures in windows (twice, once for reference aps only)
    xpos, ypos = zip(
        *((swdata.x_pixel(ap.x), swdata.y_pixel(ap.y)) for ap in ccdaper.values())
    )
    xpos = np.array(xpos)
    ypos = np.array(ypos)

    positions = Table(names=["x_0", "y_0"], data=(xpos, ypos))
    # and positions of PSF stars in the sub-window
    ref_idx = [i for i, ap in enumerate(ccdaper.values()) if ap.ref]
    psf_positions = Table(names=["x_0", "y_0"], data=(xpos[ref_idx], ypos[ref_idx]))
    # assign all PSF stars to the same group, so the same PSF params are used for all
    psf_positions["id"] = 1

    # create instances for PSF photometry
    sclip = float(rfile["sky"]["thresh"])
    bkg = LocalBackground(
        rfile["extraction"][cnam][7],  # inner sky
        rfile["extraction"][cnam][10],  # outer sky
        bkg_estimator=MedianBackground(sigma_clip=SigmaClip(sclip)),
    )
    fitshape_box_size = int(2 * int(float(rfile["psf_photom"]["fit_half_width"])) + 1)
    fit_shape = (fitshape_box_size, fitshape_box_size)

    # aperture radius is only used to guess initial flux
    aperture_radius = 0.5 * (
        rfile["extraction"][cnam][3] + rfile["extraction"][cnam][4]
    )

    photometry_task = PSFPhotometry(
        psf_model=psf_model,
        aperture_radius=aperture_radius,
        fit_shape=fit_shape,
        localbkg_estimator=bkg,
        xy_bounds=rfile["apertures"]["fit_max_shift"],
    )
    # do the PSF photometry
    photom_results = photometry_task(
        swdata.data, error=sigma, init_params=psf_positions
    )

    # now create the PSF model from fitting the reference stars
    psf_model = create_psf_model(
        photom_results, method, rfile["psf_photom"]["positions"] == "fixed"
    )

    # re-do the PSF photometry with the fixed PSF shape and all stars
    photometry_task = PSFPhotometry(
        psf_model=psf_model,
        aperture_radius=aperture_radius,
        fit_shape=fit_shape,
        localbkg_estimator=bkg,
        xy_bounds=rfile["apertures"]["fit_max_shift"],
    )
    photom_results = photometry_task(swdata.data, error=sigma, init_params=positions)

    #
    # if cnam == "4":
    #    display(swdata.data, photom_results, photometry_task)

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
                countse = result_row["flux_err"]
            except KeyError:
                raise hcam.HipercamError(
                    "unable to find errors on solution, model does not depend on params"
                )
            slevel = result_row["local_bkg"]

            # check flags from PSF photom
            if result_row["flags"] == 2:
                flag = hcam.TARGET_AT_EDGE
            elif result_row["flags"] > 8:
                # fit failed
                flag = hcam.NO_EXTRACTION

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
