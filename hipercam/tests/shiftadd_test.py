import importlib
import itertools
import warnings
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from trm.cline import Cline

import hipercam as hcam

shiftadd = importlib.import_module("hipercam.scripts.shiftadd")

pytestmark = pytest.mark.skipif(
    not shiftadd.HAS_REPROJECT, reason="reproject is not installed"
)


# tmp_path is a built-in pytest fixture that provides a
# temporary directory unique to the test invocation.
@pytest.mark.filterwarnings("ignore:Input data contains invalid values")
def test_tiled_combination_matches_full_stack(tmp_path):
    rng = np.random.default_rng(42)
    data = rng.normal(size=(11, 9, 13))
    data[0, 2, 3] = 100
    data[1:4, 4, 5] = np.nan
    data[:, 0, 0] = np.nan
    mapped = np.memmap(
        tmp_path / "stack.dat", mode="w+", dtype="float64", shape=data.shape
    )
    mapped[:] = data

    for method, sigma in (("m", None), ("c", 0), ("c", 1.5)):
        with warnings.catch_warnings():
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
            actual = shiftadd.combine_stack(
                mapped, method, sigma, 3, tile_bytes=11 * 8 * 5
            )
            np.testing.assert_allclose(actual, expected, equal_nan=True)
            np.testing.assert_array_equal(mapped, data)
            in_memory = shiftadd.combine_stack(
                [plane.copy() for plane in data],
                method,
                sigma,
                3,
                tile_bytes=11 * 8 * 5,
            )
            np.testing.assert_allclose(in_memory, expected, equal_nan=True)


def test_combination_limits_tile_size_for_many_frames(tmp_path, monkeypatch):
    shape = (300, 12, 17)
    mapped = np.memmap(tmp_path / "stack.dat", mode="w+", dtype="float64", shape=shape)
    for frame in range(shape[0]):
        mapped[frame] = frame
    sizes = []
    mean = shiftadd.meanfunc

    # monkeypatch the mean function to record the size of each tile processed
    def checked_mean(tile, axis):
        sizes.append(tile.nbytes)
        return mean(tile, axis=axis)

    monkeypatch.setattr(shiftadd, "meanfunc", checked_mean)
    budget = shape[0] * 8 * 7
    result = shiftadd.combine_stack(mapped, "c", 0, 3, tile_bytes=budget)
    np.testing.assert_allclose(result, 149.5)
    assert len(sizes) > 1
    assert max(sizes) <= budget


def test_window_mask_coordinates():
    """
    Tests the mask_reprojected function.

    The mask_reprojected function should correctly identify the coordinates
    of the window within the data array, given various offsets and window
    specifications.

    The window is defined by its lower-left corner (llx, lly), its size (nx, ny),
    and its binning factors (xbin, ybin). The lower-left corrner of the CCD
    is at (1, 1), which is index (0, 0) in the data array.

    (nx, ny) is the size in *binned* pixels.

    We make an 8x8 data array which defines the full-frame reprojection of the
    window. The *unshifted* window is then defined, and we check that the
    masked region of the data array corresponds to the window's coordinates
    after applying the specified offset (in unbinned pixels).
    """
    cases = (
        (hcam.Winhead(1, 1, 3, 2, 1, 1, ""), (0, 0), (0, 2, 0, 3)),
        (hcam.Winhead(3, 5, 3, 2, 1, 1, ""), (1, -1), (5, 7, 1, 4)),
        (hcam.Winhead(1, 1, 3, 2, 1, 1, ""), (1, 1), (0, 1, 0, 2)),
        (hcam.Winhead(1, 1, 3, 2, 1, 1, ""), (10, 0), (0, 0, 0, 0)),
        (hcam.Winhead(5, 7, 3, 2, 2, 3, ""), (2, 3), (1, 3, 1, 4)),
        (hcam.Winhead(5, 7, 3, 2, 2, 3, ""), (-2, -3), (3, 5, 3, 6)),
        (hcam.Winhead(4, 6, 3, 2, 2, 3, ""), (0, 0), (1, 3, 1, 4)),
    )
    for wind, offset, (y1, y2, x1, x2) in cases:
        data = np.arange(8 * 8, dtype=float).reshape(8, 8)
        expected = np.full_like(data, np.nan)
        expected[y1:y2, x1:x2] = data[y1:y2, x1:x2]
        shiftadd.mask_reprojected(data, wind, *offset)
        np.testing.assert_array_equal(data, expected)


def setup_run(
    tmp_path,
    monkeypatch,
    method,
    invalid=False,
    shifts=None,
    second_llx=9,
    full_frame=False,
    binning=1,
    memmap=True,
):
    """
    Mock some MCCD files and inputs for shiftadd.

    This function creates a set of mock MCCD files and input
    parameters for testing the shiftadd functionality. It also
    monkeypatches certain items within the shiftadd module to simulate the behavior of the
    shiftadd script without requiring actual command-line execution or file I/O.

    Parameters
    ----------
    tmp_path : pathlib.Path
        A temporary directory provided by pytest for storing mock files.
    monkeypatch : pytest.MonkeyPatch
        A pytest fixture for monkeypatching attributes and functions.
    method : str
        The combination method to use ('m' for median, 'c' for mean).
    invalid : bool, optional
        If True, introduce NaN values into the mock data to simulate invalid data.
    shifts : list of float, optional
        A list of shifts to apply to the apertures for each frame.
    second_llx : int, optional
        The lower-left x-coordinate for the second window in the mock data.
    full_frame : bool, optional
        If True, create full-frame data instead of windowed data.
    binning : int, optional
        The binning factor to apply to the mock data.
    memmap : bool, optional
        If True, use memory-mapped arrays for processing; otherwise, use in-memory arrays.
    """
    paths = []
    for frame in range(4):
        ccds = []
        for cnam in ("1", "2"):
            windows = []
            window_specs = (
                (("a", 1, 10),)
                if full_frame
                else (("a", 1, 10), ("b", second_llx, 100))
            )
            for wnam, llx, level in window_specs:
                nx, ny = (16 // binning, 8 // binning) if full_frame else (4, 4)
                head = hcam.Winhead(llx, 1, nx, ny, binning, binning, "")
                data = np.full((ny, nx), level + frame, dtype=float)
                if full_frame:
                    data += np.arange(nx * ny).reshape(ny, nx)
                if invalid:
                    data[1, 1] = np.nan
                window = hcam.Window(head, data)
                window["DSTATUS"] = not (cnam == "2" and frame == 1)
                windows.append((wnam, window))
            ccds.append((cnam, hcam.CCD(windows, 16, 8)))
        head = fits.Header()
        head["MJDUTC"] = 60000.0 + frame
        path = tmp_path / f"frame{frame}.hcm"
        hcam.MCCD(ccds, head).write(path)
        paths.append(path)

    flist = tmp_path / "frames.lis"
    flist.write_text("".join(f"{path}\n" for path in paths))
    outfile = tmp_path / "out.hcm"
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
        "memmap": memmap,
        "consflux": True,
        "reprkernel": "Gaussian",
        "kwidth": 1.3,
        "regwidth": 4,
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
    monkeypatch.setattr(shiftadd, "Cline", FakeCline)
    monkeypatch.setattr(shiftadd.cline, "script_args", lambda args: ("shiftadd", []))
    monkeypatch.setattr(hcam.reduction.Rfile, "read", lambda _: rfile)
    monkeypatch.setattr(
        hcam.reduction,
        "initial_checks",
        lambda *args: ({"1": None}, {"1": None}, True),
    )
    calls = itertools.cycle([1.0, 1.0, 3.0, 1.0])
    offsets = itertools.cycle(shifts or [0.0] * 4)

    def move_apertures(*args):
        store = args[-1]
        store["mfwhm"] = next(calls)
        store["ref"] = {"dx": next(offsets), "dy": 0.0}

    monkeypatch.setattr(hcam.reduction, "moveApers", move_apertures)
    monkeypatch.setattr(shiftadd.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(shiftadd.signal, "signal", lambda *args: None)
    return outfile, values


@pytest.mark.filterwarnings("ignore:Input data contains invalid values")
def test_shiftadd_windows_and_missing_ccd(tmp_path, monkeypatch):
    outfile, _ = setup_run(tmp_path, monkeypatch, "c")
    shiftadd.shiftadd()
    result = hcam.MCCD.read(str(outfile))
    for cnam, frame_mean in (("1", 4 / 3), ("2", 1.5)):
        for wnam, level in (("a", 10), ("b", 100)):
            np.testing.assert_allclose(result[cnam][wnam].data, level + frame_mean)
    assert result.head["MJDUTC"] == 60001.5
    assert "CCD1(3),CCD2(2)" in str(result.head["HISTORY"])
    assert not list(tmp_path.glob("tmp*"))


def test_shifted_windows_match_combination(tmp_path, monkeypatch):
    outfile, _ = setup_run(
        tmp_path, monkeypatch, "m", shifts=[0, 1, -1, 0], second_llx=5
    )
    shiftadd.shiftadd()
    result = hcam.MCCD.read(str(outfile))
    orig_wcs = WCS(naxis=2)
    offsets = (0, 1, 0, 0)
    for cnam in ("1", "2"):
        planes = []
        for frame in (0, 1, 3):
            ccd = hcam.MCCD.read(str(tmp_path / f"frame{frame}.hcm"))[cnam]
            if not ccd.is_data():
                continue
            for wind in ccd.values():
                input_wcs = orig_wcs.deepcopy()
                input_wcs.wcs.crpix = [
                    offsets[frame] - (wind.llx - 1) // wind.xbin,
                    -(wind.lly - 1) // wind.ybin,
                ]
                data, _ = shiftadd.reproject_interp(
                    (wind.data, input_wcs), orig_wcs, (8, 16), order=0
                )
                mask = np.zeros(data.shape, dtype=bool)
                x = (wind.llx - 1) // wind.xbin
                y = (wind.lly - 1) // wind.ybin
                mask[
                    y : y + wind.ny,
                    max(0, x - offsets[frame]) : x + wind.nx - offsets[frame],
                ] = True
                data[~mask] = np.nan
                planes.append(data)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            expected = shiftadd.medianfunc(np.stack(planes), axis=0)
        for wind in result[cnam].values():
            x = (wind.llx - 1) // wind.xbin
            y = (wind.lly - 1) // wind.ybin
            np.testing.assert_allclose(
                wind.data, expected[y : y + wind.ny, x : x + wind.nx]
            )


def test_full_frame_preserves_pixels_at_all_edges(tmp_path, monkeypatch):
    outfile, _ = setup_run(tmp_path, monkeypatch, "c", full_frame=True)
    shiftadd.shiftadd()
    result = hcam.MCCD.read(str(outfile))
    expected = np.arange(128).reshape(8, 16)
    for cnam, frame_mean in (("1", 4 / 3), ("2", 1.5)):
        assert result[cnam]["a"].data.shape == (8, 16)
        np.testing.assert_allclose(result[cnam]["a"].data, 10 + frame_mean + expected)


@pytest.mark.filterwarnings("ignore:Input data contains invalid values")
def test_binned_full_frame_shift(tmp_path, monkeypatch):
    outfile, _ = setup_run(
        tmp_path, monkeypatch, "c", full_frame=True, binning=2, shifts=[0, 2, -2, 0]
    )
    shiftadd.shiftadd()
    result = hcam.MCCD.read(str(outfile))
    gradient = np.arange(32).reshape(4, 8)
    expected = gradient + 35 / 3
    expected[:, -1] = gradient[:, -1] + 11.5
    np.testing.assert_allclose(result["1"]["a"].data, expected)
    np.testing.assert_allclose(result["2"]["a"].data, gradient + 11.5)


@pytest.mark.filterwarnings("ignore:Input data contains invalid values")
def test_memory_and_memmap_produce_same_results(tmp_path, monkeypatch):
    outfile, values = setup_run(
        tmp_path, monkeypatch, "m", shifts=[0, 1, -1, 0], second_llx=5
    )
    celestial_wcs = WCS(naxis=2)
    celestial_wcs.wcs.crpix = [4, 5]
    celestial_wcs.wcs.crval = [40, 30]
    celestial_wcs.wcs.cdelt = [-0.01, 0.01]
    celestial_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    monkeypatch.setattr(
        shiftadd, "wcs_from_header", lambda *args, **kwargs: celestial_wcs
    )
    for reprmethod in ("interp", "adaptive", "exact"):
        for method, sigma in (("m", 1.5), ("c", 1.5), ("c", 0)):
            values.update(
                reprmethod=reprmethod, method=method, sigma=sigma, memmap=True
            )
            values["overwrite"] = True
            shiftadd.shiftadd()
            disk = hcam.MCCD.read(str(outfile))
            values["memmap"] = False

            def check_memory(function):
                def run(*args, **kwargs):
                    assert type(kwargs["output_array"]) is np.ndarray
                    assert "block_size" not in kwargs
                    assert not kwargs["return_footprint"]
                    return function(*args, **kwargs)

                return run

            name = "reproject_" + reprmethod
            with monkeypatch.context() as context:
                context.setattr(shiftadd, name, check_memory(getattr(shiftadd, name)))
                shiftadd.shiftadd()
            memory = hcam.MCCD.read(str(outfile))
            for cnam in disk:
                for wnam in disk[cnam]:
                    np.testing.assert_allclose(
                        memory[cnam][wnam].data,
                        disk[cnam][wnam].data,
                        rtol=1e-8,
                        atol=1e-8,
                    )
            assert memory.head["MJDUTC"] == disk.head["MJDUTC"]
            assert str(memory.head["HISTORY"]) == str(disk.head["HISTORY"])
            assert not list(tmp_path.glob("tmp*"))


@pytest.mark.filterwarnings("ignore:Input data contains invalid values")
def test_hidden_memmap_argument_with_real_cline(tmp_path, monkeypatch):
    setup_run(tmp_path, monkeypatch, "m")
    rfile = tmp_path / "dummy.red"
    rfile.touch()
    monkeypatch.setattr(shiftadd, "Cline", Cline)
    function = shiftadd.reproject_interp
    observed = []

    def reproject(*args, **kwargs):
        observed.append(isinstance(kwargs["output_array"], np.memmap))
        return function(*args, **kwargs)

    monkeypatch.setattr(shiftadd, "reproject_interp", reproject)
    for use_memmap in (True, False):
        outfile = tmp_path / f"cline-{use_memmap}.hcm"
        args = [
            "nodefs",
            "source=hf",
            str(tmp_path / "frames.lis"),
            str(rfile),
            "1",
            "interp",
            "no",
            str(outfile),
        ]
        if not use_memmap:
            args.append("memmap=no")
        with monkeypatch.context() as context:
            context.setattr(
                shiftadd.cline, "script_args", lambda *a, **kw: (None, args)
            )
            shiftadd.shiftadd()
        assert outfile.exists()
        assert observed
        assert all(value == use_memmap for value in observed)
        observed.clear()


def test_nan_failure_cleans_temporary_files(tmp_path, monkeypatch):
    outfile, values = setup_run(tmp_path, monkeypatch, "m", invalid=True)
    for memmap in (True, False):
        values["memmap"] = memmap
        with pytest.raises(hcam.HipercamError, match="NaN values detected"):
            shiftadd.shiftadd()
        assert not outfile.exists()
        assert not list(tmp_path.glob("tmp*"))


def test_reprojection_failure_cleans_temporary_files(tmp_path, monkeypatch):
    outfile, _ = setup_run(tmp_path, monkeypatch, "m")

    def fail_reproject(*args, **kwargs):
        assert isinstance(kwargs["output_array"], np.memmap)
        assert not kwargs["return_footprint"]
        raise RuntimeError("reprojection failed")

    monkeypatch.setattr(shiftadd, "reproject_interp", fail_reproject)
    with pytest.raises(RuntimeError, match="reprojection failed"):
        shiftadd.shiftadd()
    assert not outfile.exists()
    assert not list(tmp_path.glob("tmp*"))
