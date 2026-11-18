"""Tests of Gaussian and Moffat profile evaluation and fitting."""

import numpy as np
import pytest

from hipercam import HipercamError, Window, Winhead
from hipercam import fitting


def analytic_profile(method, x, y, sky, height, xcen, ycen, fwhm, beta=2.5):
    """Evaluate the profiles independently of the compiled fitting functions."""
    radius2 = (x - xcen) ** 2 + (y - ycen) ** 2
    if method == "g":
        return sky + height * np.exp(-4 * np.log(2) * radius2 / fwhm**2)
    alpha = 4 * (2 ** (1 / beta) - 1) / fwhm**2
    return sky + height / (1 + alpha * radius2) ** beta


@pytest.mark.parametrize("method", ["g", "m"])
@pytest.mark.parametrize("xbin,ybin,ndiv", [(1, 1, 0), (2, 3, 0), (2, 3, 2)])
def test_profile_evaluation(method, xbin, ybin, ndiv):
    x, y = np.meshgrid(np.array([2.0, 3.2, 4.7]), np.array([4.1, 6.0]))
    pars = (12.0, 150.0, 3.1, 4.7, 3.8)
    beta = 2.5
    # ndiv=0 samples pixel centres, even when binned; otherwise average
    # ndiv samples per unbinned pixel in each direction.
    if ndiv:
        xoffsets = (np.arange(xbin * ndiv) + 0.5) / ndiv - xbin / 2
        yoffsets = (np.arange(ybin * ndiv) + 0.5) / ndiv - ybin / 2
    else:
        xoffsets = yoffsets = [0.0]
    expected = np.mean(
        [
            analytic_profile(method, x + dx, y + dy, *pars, beta)
            for dx in xoffsets
            for dy in yoffsets
        ],
        axis=0,
    )
    if method == "g":
        actual = fitting.gaussian(x, y, *pars, xbin, ybin, ndiv)
    else:
        actual = fitting.moffat(x, y, *pars, beta, xbin, ybin, ndiv)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


def make_window(method, sky=12.0):
    # The non-unit origin and off-centre star exercise absolute CCD coordinates.
    win = Winhead(11, 21, 21, 21, 1, 1)
    x, y = np.meshgrid(win.x(np.arange(win.nx)), win.y(np.arange(win.ny)))
    # Parameter order: sky, height, x, y, FWHM, Moffat beta.
    pars = (sky, 150.0, 21.3, 30.7, 4.2, 2.5)
    data = analytic_profile(method, x, y, *pars)
    # A small deterministic residual avoids a perfect, zero-variance fit.
    data += 0.01 * np.sin(x * 0.7 + y * 0.4)
    return Window(win, data), pars


@pytest.mark.parametrize(
    "method,mode",
    # Mode letters select free sky, FWHM and beta respectively.
    [("g", mode) for mode in ("", "s", "f", "sf")]
    + [("m", mode) for mode in ("", "s", "f", "b", "sf", "sb", "fb", "sfb")],
)
@pytest.mark.parametrize("ndiv", [0, 2])
def test_residual_jacobian_matches_finite_differences(method, mode, ndiv):
    win, pars = make_window(method)
    sigma = np.full(win.data.shape, 2.0)
    # Negative uncertainties mark pixels excluded from the fit.
    sigma[10, 10] = -1.0
    if method == "g":
        model = fitting.Gfit(win, sigma, ndiv, mode, pars[4])
        params = np.array(model.set_par(*pars[:5]))
    else:
        model = fitting.Mfit(win, sigma, ndiv, mode, pars[4], pars[5])
        params = np.array(model.set_par(*pars))

    # Compare the compiled residual Jacobian with central finite differences
    # for every free parameter, including the effect of masking and weights.
    numerical = []
    for i, value in enumerate(params):
        step = 1e-5 * max(1.0, abs(value))
        delta = np.zeros_like(params)
        delta[i] = step
        numerical.append(
            (model.fun(params + delta) - model.fun(params - delta)) / (2 * step)
        )
    np.testing.assert_allclose(
        model.jac(params), np.column_stack(numerical), rtol=1e-6, atol=1e-7
    )
    assert model.fun(params).size == np.count_nonzero(model.mask & (sigma > 0))


@pytest.mark.parametrize("method", ["g", "m"])
@pytest.mark.parametrize("fit_sky", [False, True])
@pytest.mark.parametrize("fixed_fwhm", [False, True])
@pytest.mark.parametrize("fixed_beta", [False, True])
def test_combfit_recovers_synthetic_profile(method, fit_sky, fixed_fwhm, fixed_beta):
    win, truth = make_window(method, sky=12.0 if fit_sky else 0.0)
    sigma = np.ones(win.data.shape)
    sigma[10, 10] = -1.0
    result, errors, extras = fitting.combFit(
        win, sigma, method,
        10.0 if fit_sky else None, 140.0, 21.0, 31.0,
        truth[4] if fixed_fwhm else 4.0, 1.0, fixed_fwhm,
        truth[5] if fixed_beta else 2.3, 10.0, fixed_beta,
        100.0,
    )
    # combFit returns beta=0 for Gaussian fits. Fixed parameters have an
    # uncertainty of -1; an unfitted sky is zero with uncertainty -1.
    expected = truth if method == "m" else (*truth[:5], 0.0)
    # Height has a larger absolute tolerance because its scale is much larger.
    np.testing.assert_allclose(result[1], expected[1], rtol=0, atol=0.03)
    np.testing.assert_allclose(
        [result[i] for i in (0, 2, 3, 4, 5)],
        [expected[i] for i in (0, 2, 3, 4, 5)],
        rtol=0,
        atol=0.005,
    )
    if fit_sky:
        assert errors[0] >= 0
    else:
        assert errors[0] == -1
    if fixed_fwhm:
        assert errors[4] == -1
    else:
        assert errors[4] >= 0
    if method == "g" or fixed_beta:
        assert errors[5] == -1
    else:
        assert errors[5] >= 0
    fit, x, y, chisq, nok, nrej, npar, nfev, message = extras
    assert isinstance(fit, Window)
    assert x.shape == y.shape == fit.data.shape == win.data.shape
    assert chisq >= 0 and np.isfinite(chisq)
    assert nok == np.count_nonzero(sigma > 0)
    # The pre-excluded pixel is counted as rejected within the circular mask.
    assert nrej == 1
    # Height and position are always free; sky, FWHM and beta are optional.
    assert npar == 3 + fit_sky + (not fixed_fwhm) + (method == "m" and not fixed_beta)
    assert nfev > 0 and message
    assert sigma[10, 10] < 0
    np.testing.assert_allclose(fit.data, win.data, atol=0.03)


@pytest.mark.parametrize("method", ["g", "m"])
def test_combfit_recovers_binned_subpixel_profile(method):
    win = Winhead(11, 21, 21, 21, 2, 3)
    x, y = np.meshgrid(win.x(np.arange(win.nx)), win.y(np.arange(win.ny)))
    truth = (12.0, 150.0, 31.7, 51.3, 5.2, 2.5)
    # For 2x3 binning with ndiv=2, average over a 4x6 subpixel grid.
    offsets_x = (np.arange(4) + 0.5) / 2 - 1
    offsets_y = (np.arange(6) + 0.5) / 2 - 1.5
    data = np.mean(
        [
            analytic_profile(method, x + dx, y + dy, *truth)
            for dx in offsets_x
            for dy in offsets_y
        ],
        axis=0,
    )
    data += 0.01 * np.sin(x * 0.7 + y * 0.4)
    wind = Window(win, data)
    result, _, _ = fitting.combFit(
        wind, np.ones(data.shape), method, 10.0, 140.0, 31.0, 52.0,
        5.0, 1.0, False, 2.3, 10.0, False, 100.0, ndiv=2,
    )
    expected = truth if method == "m" else (*truth[:5], 0.0)
    np.testing.assert_allclose(result, expected, rtol=0, atol=0.05)


@pytest.mark.parametrize("method", ["g", "m"])
def test_combfit_constrains_fwhm_at_minimum(method):
    win, _ = make_window(method)
    # The injected FWHM is 4.2, so the unconstrained fit hits this limit.
    result, errors, extras = fitting.combFit(
        win, np.ones(win.data.shape), method, 12.0, 140.0, 21.0, 31.0,
        5.0, 4.5, False, 1.8, 10.0, False, 100.0,
    )
    assert result[4] == 4.5
    assert errors[4] == -1
    assert extras[6] == (4 if method == "g" else 5)


def test_combfit_constrains_beta_at_maximum():
    win, _ = make_window("m")
    # The injected beta is 2.5; the fit must switch to fixed beta=2.
    result, errors, extras = fitting.combFit(
        win, np.ones(win.data.shape), "m", 12.0, 140.0, 21.0, 31.0,
        5.0, 1.0, False, 1.8, 2.0, False, 100.0,
    )
    assert result[5] == 2.0
    assert errors[5] == -1
    assert extras[6] == 5


@pytest.mark.parametrize("method", ["g", "m"])
@pytest.mark.xfail(
    strict=True,
    reason="rejected pixels remain in cached least-squares indices",
)
def test_combfit_rejects_outlier_and_recovers_profile(method):
    win, truth = make_window(method)
    # A single large contaminant should be rejected before the final fit.
    win.data[9, 12] += 500.0
    sigma = np.ones(win.data.shape)
    result, _, extras = fitting.combFit(
        win, sigma, method, 10.0, 140.0, 21.0, 31.0,
        4.0, 1.0, False, 2.3, 10.0, False, 4.0,
    )
    expected = truth if method == "m" else (*truth[:5], 0.0)
    assert sigma[9, 12] < 0
    assert extras[5] == 1
    np.testing.assert_allclose(result, expected, rtol=0, atol=0.05)


@pytest.mark.parametrize("method", ["g", "m"])
def test_combfit_rejects_initial_fwhm_below_minimum(method):
    win, _ = make_window(method)
    with pytest.raises(HipercamError, match="fwhm"):
        fitting.combFit(win, np.ones(win.data.shape), method, 10, 140, 21, 31,
                        1.0, 2.0, False, 2.5, 10, False, 100)


def test_combfit_rejects_initial_beta_above_maximum():
    win, _ = make_window("m")
    with pytest.raises(HipercamError, match="beta"):
        fitting.combFit(win, np.ones(win.data.shape), "m", 10, 140, 21, 31,
                        4.0, 1.0, False, 11, 10, False, 100)


def test_combfit_rejects_unknown_method():
    win, _ = make_window("g")
    with pytest.raises(NotImplementedError, match="Fitting method"):
        fitting.combFit(win, np.ones(win.data.shape), "unknown", 10, 140, 21, 31,
                        4.0, 1.0, False, 2.5, 10, False, 100)
