import unittest

import numpy as np

from hipercam import Winhead
from hipercam.fitting import combFit, fitFastGaussian, gaussian
from hipercam.window import Window


class TestFastGaussian(unittest.TestCase):
    """
    Compares the fast gaussian fit against the standard gaussian fit on
    noise-free (infinite signal-to-noise) data, where both should recover the
    true parameters.
    """

    def setUp(self):
        nx = ny = 31
        self.truth = dict(sky=100.0, height=2000.0, x=15.3, y=16.4, fwhm=4.0)
        x, y = np.meshgrid(np.arange(nx) + 1.0, np.arange(ny) + 1.0)
        t = self.truth
        data = gaussian(
            x, y, t["sky"], t["height"], t["x"], t["y"], t["fwhm"], 1, 1, 0
        )
        self.wind = Window(Winhead(1, 1, nx, ny, 1, 1, "LL"), data)
        self.sigma = np.full(data.shape, 5.0)

    def fit(self, method, fwhm_fix=False, fwhm=3.5):
        t = self.truth
        # deliberately imperfect starting values
        pars, epars, extras = combFit(
            self.wind,
            self.sigma.copy(),
            method,
            t["sky"] - 5.0,
            0.8 * t["height"],
            t["x"] + 0.7,
            t["y"] - 0.6,
            fwhm,
            1.0,
            fwhm_fix,
            0.0,
            0.0,
            True,
            4.0,
            0,
        )
        return pars

    def test_matches_standard_gaussian(self):
        sky_g, h_g, x_g, y_g, fwhm_g, _ = self.fit("g")
        sky_f, h_f, x_f, y_f, fwhm_f, beta_f = self.fit("f")
        self.assertAlmostEqual(sky_f, sky_g, delta=0.1)
        self.assertAlmostEqual(h_f, h_g, delta=0.005 * h_g)
        self.assertAlmostEqual(x_f, x_g, delta=0.01)
        self.assertAlmostEqual(y_f, y_g, delta=0.01)
        self.assertAlmostEqual(fwhm_f, fwhm_g, delta=0.01)
        self.assertEqual(beta_f, 0.0)

    def test_recovers_truth(self):
        t = self.truth
        sky, h, x, y, fwhm, _ = self.fit("f")
        self.assertAlmostEqual(sky, t["sky"], delta=0.1)
        self.assertAlmostEqual(h, t["height"], delta=0.005 * t["height"])
        self.assertAlmostEqual(x, t["x"], delta=0.01)
        self.assertAlmostEqual(y, t["y"], delta=0.01)
        self.assertAlmostEqual(fwhm, t["fwhm"], delta=0.01)

    def test_two_pass_matches_standard_gaussian(self):
        t = self.truth
        sky_g, h_g, x_g, y_g, fwhm_g, _ = self.fit("g")
        (sky, h, x, y, fwhm), _, _ = fitFastGaussian(
            self.wind, self.sigma.copy(), t["sky"] - 5.0, 0.8 * t["height"],
            t["x"] + 0.7, t["y"] - 0.6, 3.5, 1.0, False, 4.0, 0,
            snr_min=3.0, direct_snr=False,
        )
        self.assertAlmostEqual(sky, sky_g, delta=0.1)
        self.assertAlmostEqual(h, h_g, delta=0.005 * h_g)
        self.assertAlmostEqual(x, x_g, delta=0.01)
        self.assertAlmostEqual(y, y_g, delta=0.01)
        self.assertAlmostEqual(fwhm, fwhm_g, delta=0.01)

    def test_fixed_fwhm(self):
        t = self.truth
        # fix the FWHM at its true value
        sky_f, h_f, x_f, y_f, fwhm_f, _ = self.fit(
            "f", fwhm_fix=True, fwhm=t["fwhm"]
        )
        sky_g, h_g, x_g, y_g, fwhm_g, _ = self.fit(
            "g", fwhm_fix=True, fwhm=t["fwhm"]
        )
        self.assertEqual(fwhm_f, t["fwhm"])
        self.assertAlmostEqual(sky_f, sky_g, delta=0.1)
        self.assertAlmostEqual(h_f, h_g, delta=0.005 * h_g)
        self.assertAlmostEqual(x_f, x_g, delta=0.01)
        self.assertAlmostEqual(y_f, y_g, delta=0.01)


class TestFastGaussianOutliers(unittest.TestCase):
    """Bright off-source pixels (cosmic rays) must not wreck the fast fit."""

    def test_cosmic_rays(self):
        rng = np.random.default_rng(1)
        nx = ny = 31
        t = dict(sky=100.0, height=500.0, x=15.3, y=16.4, fwhm=4.0)
        x, y = np.meshgrid(np.arange(nx) + 1.0, np.arange(ny) + 1.0)
        model = gaussian(x, y, t["sky"], t["height"], t["x"], t["y"], t["fwhm"], 1, 1, 0)
        sigma = np.sqrt(model + 25.0)
        for _ in range(30):
            data = model + sigma * rng.standard_normal(model.shape)
            for _ in range(3):
                data[rng.integers(0, ny), rng.integers(0, nx)] += 2000.0
            wind = Window(Winhead(1, 1, nx, ny, 1, 1, "LL"), data)
            pars, _, _ = fitFastGaussian(
                wind, sigma.copy(), 95.0, 400.0, 15.0, 16.0, 3.5, 1.0, False, 4.0, 0
            )
            self.assertAlmostEqual(pars[2], t["x"], delta=0.3)
            self.assertAlmostEqual(pars[3], t["y"], delta=0.3)
            self.assertAlmostEqual(pars[4], t["fwhm"], delta=0.8)


if __name__ == "__main__":
    unittest.main()
