import importlib
import tempfile
import unittest
import warnings
from contextlib import ExitStack
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.wcs import WCS

import hipercam as hcam

shiftadd = importlib.import_module("hipercam.scripts.shiftadd")


@unittest.skipUnless(shiftadd.HAS_REPROJECT, "reproject is not installed")
class TestShiftadd(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.tmp = Path(self.stack.enter_context(tempfile.TemporaryDirectory()))

    def test_tiled_combination_matches_full_stack(self):
        rng = np.random.default_rng(42)
        data = rng.normal(size=(11, 9, 13))
        data[0, 2, 3] = 100
        data[1:4, 4, 5] = np.nan
        data[:, 0, 0] = np.nan
        mapped = np.memmap(
            self.tmp / "stack.dat", mode="w+", dtype="float64", shape=data.shape
        )
        mapped[:] = data

        for method, sigma in (("m", None), ("c", 0), ("c", 1.5)):
            with self.subTest(method=method, sigma=sigma), warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                expected_data = data.copy()
                if method == "m":
                    expected = shiftadd.medianfunc(expected_data, axis=0)
                else:
                    if sigma > 0:
                        mask = sigma_clip(
                            expected_data,
                            sigma_lower=sigma,
                            sigma_upper=sigma,
                            axis=0,
                            copy=False,
                            maxiters=3,
                            cenfunc="mean",
                            stdfunc="std",
                            masked=True,
                        )
                        expected_data[mask.mask] = np.nan
                    expected = shiftadd.meanfunc(expected_data, axis=0)
                actual = shiftadd._combine_stack(
                    mapped, method, sigma, 3, tile_bytes=11 * 8 * 5
                )
                np.testing.assert_allclose(actual, expected, equal_nan=True)
                np.testing.assert_array_equal(mapped, data)

    def test_combination_limits_tile_size_for_many_frames(self):
        shape = (300, 12, 17)
        mapped = np.memmap(
            self.tmp / "stack.dat", mode="w+", dtype="float64", shape=shape
        )
        for frame in range(shape[0]):
            mapped[frame] = frame
        sizes = []
        mean = shiftadd.meanfunc

        def checked_mean(tile, axis):
            sizes.append(tile.nbytes)
            return mean(tile, axis=axis)

        self.stack.enter_context(patch.object(shiftadd, "meanfunc", checked_mean))
        budget = shape[0] * 8 * 7
        result = shiftadd._combine_stack(mapped, "c", 0, 3, tile_bytes=budget)
        np.testing.assert_allclose(result, 149.5)
        self.assertGreater(len(sizes), 1)
        self.assertLessEqual(max(sizes), budget)

    def test_window_mask_matches_original_slicing(self):
        wind = hcam.Winhead(3, 5, 6, 4, 1, 1, "")
        for offset in ((0, 0), (1.2, -2.7), (10, 0), (-20, 3)):
            with self.subTest(offset=offset):
                data = np.arange(12 * 14, dtype=float).reshape(12, 14)
                expected = data.copy()
                mask = np.zeros_like(data, dtype=bool)
                xstart = wind.llx // wind.xbin
                ystart = wind.lly // wind.ybin
                mask[
                    ystart - int(offset[1]) : ystart + wind.ny - int(offset[1]),
                    xstart - int(offset[0]) : xstart + wind.nx - int(offset[0]),
                ] = True
                expected[~mask] = np.nan
                shiftadd._mask_reprojected(data, wind, *offset)
                np.testing.assert_array_equal(data, expected)

    def test_block_reprojection_matches_full_image(self):
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [4, 5]
        wcs.wcs.crval = [40, 30]
        wcs.wcs.cdelt = [-0.01, 0.01]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        input_wcs = wcs.deepcopy()
        input_wcs.wcs.crpix += [0.3, -0.2]
        image = np.random.default_rng(42).random((10, 12)).astype(np.float32)
        shape = (19, 21)
        for name in ("interp", "adaptive", "exact"):
            with self.subTest(method=name):
                function = getattr(shiftadd, "reproject_" + name)
                options = {"boundary_mode": "nearest"} if name == "adaptive" else {}
                expected, _ = function((image, input_wcs), wcs, shape, **options)
                mapped = np.memmap(
                    self.tmp / "reproject.dat",
                    mode="w+",
                    dtype="float64",
                    shape=shape,
                )
                function(
                    (image.astype(np.float64), input_wcs),
                    wcs,
                    shape,
                    output_array=mapped,
                    return_footprint=False,
                    block_size=(7, 8),
                    **options,
                )
                np.testing.assert_allclose(
                    mapped, expected, rtol=1e-8, atol=1e-8, equal_nan=True
                )

    def setup_run(self, method, invalid=False, shifts=None, second_llx=9):
        paths = []
        for frame in range(4):
            ccds = []
            for cnam in ("1", "2"):
                windows = []
                for wnam, llx, level in (("a", 1, 10), ("b", second_llx, 100)):
                    head = hcam.Winhead(llx, 1, 4, 4, 1, 1, "")
                    data = np.full((4, 4), level + frame, dtype=float)
                    if invalid:
                        data[1, 1] = np.nan
                    window = hcam.Window(head, data)
                    window["DSTATUS"] = not (cnam == "2" and frame == 1)
                    windows.append((wnam, window))
                ccds.append((cnam, hcam.CCD(windows, 16, 8)))
            head = fits.Header()
            head["MJDUTC"] = 60000.0 + frame
            path = self.tmp / f"frame{frame}.hcm"
            hcam.MCCD(ccds, head).write(path)
            paths.append(path)

        flist = self.tmp / "frames.lis"
        flist.write_text("".join(f"{path}\n" for path in paths))
        outfile = self.tmp / "out.hcm"
        values = {
            "source": "hf",
            "flist": str(flist),
            "rfilen": "dummy.red",
            "refccd": "1",
            "fthresh": 2.0,
            "reprmethod": "interp",
            "reprorder": 0,
            "trim": False,
            "method": method,
            "sigma": 1.5,
            "maxiters": 3,
            "overwrite": False,
            "output": str(outfile),
        }

        class FakeCline:
            GLOBAL = LOCAL = HIDE = PROMPT = None

            def __init__(self, *args):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

            def register(self, *args):
                pass

            def set_default(self, *args):
                pass

            def get_value(self, name, *args, **kwargs):
                return values[name]

        aperture = SimpleNamespace(x=2, y=2, ref=True)
        rfile = SimpleNamespace(aper={"1": {"ref": aperture}})
        self.stack.enter_context(patch.object(shiftadd, "Cline", FakeCline))
        self.stack.enter_context(
            patch.object(shiftadd.cline, "script_args", lambda args: ("shiftadd", []))
        )
        self.stack.enter_context(
            patch.object(hcam.reduction.Rfile, "read", lambda _: rfile)
        )
        self.stack.enter_context(
            patch.object(
                hcam.reduction,
                "initial_checks",
                lambda *args: ({"1": None}, {"1": None}, True),
            )
        )
        calls = iter([1.0, 1.0, 3.0, 1.0])
        offsets = iter(shifts or [0.0] * 4)

        def move_apertures(*args):
            store = args[-1]
            store["mfwhm"] = next(calls)
            store["ref"] = {"dx": next(offsets), "dy": 0.0}

        self.stack.enter_context(
            patch.object(hcam.reduction, "moveApers", move_apertures)
        )
        self.stack.enter_context(
            patch.object(shiftadd.tempfile, "tempdir", str(self.tmp))
        )
        self.stack.enter_context(
            patch.object(shiftadd.signal, "signal", lambda *args: None)
        )
        return outfile

    def test_shiftadd_windows_and_missing_ccd(self):
        outfile = self.setup_run("c")
        shiftadd.shiftadd()
        result = hcam.MCCD.read(str(outfile))
        for cnam, frame_mean in (("1", 4 / 3), ("2", 1.5)):
            for wnam, level in (("a", 10), ("b", 100)):
                np.testing.assert_allclose(result[cnam][wnam].data, level + frame_mean)
        self.assertEqual(result.head["MJDUTC"], 60001.5)
        self.assertIn("CCD1(3),CCD2(2)", str(result.head["HISTORY"]))
        self.assertFalse(list(self.tmp.glob("tmp*")))

    def test_shifted_windows_match_original_combination(self):
        outfile = self.setup_run("m", shifts=[0, 1, -1, 0], second_llx=5)
        shiftadd.shiftadd()
        result = hcam.MCCD.read(str(outfile))
        orig_wcs = WCS(naxis=2)
        offsets = (0, 1, 0, 0)
        for cnam in ("1", "2"):
            planes = []
            for frame in (0, 1, 3):
                ccd = hcam.MCCD.read(str(self.tmp / f"frame{frame}.hcm"))[cnam]
                if not ccd.is_data():
                    continue
                for wind in ccd.values():
                    input_wcs = orig_wcs.deepcopy()
                    input_wcs.wcs.crpix = [
                        offsets[frame] - wind.llx // wind.xbin,
                        -wind.lly // wind.ybin,
                    ]
                    data, _ = shiftadd.reproject_interp(
                        (wind.data, input_wcs), orig_wcs, (8, 16), order=0
                    )
                    mask = np.zeros(data.shape, dtype=bool)
                    x, y = wind.llx // wind.xbin, wind.lly // wind.ybin
                    mask[
                        y : y + wind.ny,
                        x - offsets[frame] : x + wind.nx - offsets[frame],
                    ] = True
                    data[~mask] = np.nan
                    planes.append(data)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                expected = shiftadd.medianfunc(np.stack(planes), axis=0)
            for wind in result[cnam].values():
                x, y = wind.llx // wind.xbin, wind.lly // wind.ybin
                np.testing.assert_allclose(
                    wind.data, expected[y : y + wind.ny, x : x + wind.nx]
                )

    def test_nan_failure_cleans_temporary_files(self):
        outfile = self.setup_run("m", invalid=True)
        with self.assertRaisesRegex(hcam.HipercamError, "NaN values detected"):
            shiftadd.shiftadd()
        self.assertFalse(outfile.exists())
        self.assertFalse(list(self.tmp.glob("tmp*")))

    def test_reprojection_failure_cleans_temporary_files(self):
        outfile = self.setup_run("m")

        def fail_reproject(*args, **kwargs):
            self.assertIsInstance(kwargs["output_array"], np.memmap)
            self.assertFalse(kwargs["return_footprint"])
            self.assertEqual(kwargs["block_size"], (256, 256))
            raise RuntimeError("reprojection failed")

        self.stack.enter_context(
            patch.object(shiftadd, "reproject_interp", fail_reproject)
        )
        with self.assertRaisesRegex(RuntimeError, "reprojection failed"):
            shiftadd.shiftadd()
        self.assertFalse(outfile.exists())
        self.assertFalse(list(self.tmp.glob("tmp*")))


if __name__ == "__main__":
    unittest.main()
