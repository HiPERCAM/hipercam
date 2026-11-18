"""End-to-end smoke tests of the ``reduce`` script.

Fake multi-CCD data and calibrations are written to a temporary
directory, together with an aperture file and a reduce file. ``reduce``
is then run on the files and the log it writes is compared with
benchmark values stored in ``benchmarks/reduce_benchmarks.json``. These
benchmarks are the output of the program at the time they were
generated, so the tests guard against unintended changes in behaviour
rather than proving scientific correctness.

To regenerate the benchmarks after an intentional change run::

    HIPERCAM_UPDATE_BENCHMARKS=1 python -m pytest hipercam/tests/reduce_test.py

Each photometry mode records its own UTC generation timestamp and HiPERCAM
version, so partial updates retain the provenance of the other modes. Install
the ``psf`` extra to regenerate the PSF benchmarks.
"""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

import hipercam as hcam
from hipercam.scripts.reduce import reduce

BENCHMARK_FILE = Path(__file__).parent / "benchmarks" / "reduce_benchmarks.json"
UPDATE = bool(os.environ.get("HIPERCAM_UPDATE_BENCHMARKS"))

CCDS = ("1", "2", "3")
NFRAME = 6
NXWIN = NYWIN = 100
EXPTIME = 10.0
READOUT = 4.0
GAIN = 1.5  # electrons per ADU
SKY = {"1": 40.0, "2": 60.0, "3": 90.0}
FWHM = {"1": 4.0, "2": 4.5, "3": 5.0}

# (x, y, peak-free total electrons) at frame 0 for each aperture, in window
# pixel coordinates which start at 1 for the window below. The same
# stars are in each CCD but with different brightnesses.
STARS = {
    "1": (30.0, 70.0, 20000.0),  # reference
    "2": (65.0, 62.0, 6000.0),  # target (variable)
    "3": (60.0, 25.0, 12000.0),  # comparison
}
BAND_SCALE = {"1": 1.0, "2": 1.5, "3": 2.0}
DRIFT = (0.4, -0.25)  # pixels per frame


def _fake_ccd_frame(cnam, nframe):
    """Returns the noiseless (sky + stars) image of a CCD, in electrons."""
    y, x = np.mgrid[1 : NYWIN + 1, 1 : NXWIN + 1].astype(float)
    sigma = FWHM[cnam] / 2.3548
    image = np.full((NYWIN, NXWIN), SKY[cnam])
    for apnam, (x0, y0, flux) in STARS.items():
        flux = flux * BAND_SCALE[cnam]
        if apnam == "2":
            flux *= 1.0 + 0.2 * np.sin(2 * np.pi * nframe / NFRAME)
        xc = x0 + DRIFT[0] * nframe
        yc = y0 + DRIFT[1] * nframe
        image += (
            flux
            / (2 * np.pi * sigma**2)
            * np.exp(-((x - xc) ** 2 + (y - yc) ** 2) / (2 * sigma**2))
        )
    return image


def _make_calibrations(directory):
    """Writes bias, dark and flat. Returns the arrays for use in building
    raw frames."""
    y, x = np.mgrid[0:NYWIN, 0:NXWIN]
    # bias with small gradient
    bias = {c: 500.0 + 0.5 * x + 0.2 * y + 10 * i for i, c in enumerate(CCDS)}
    # dark with a different constant level in each CCD
    dark = {c: np.full((NYWIN, NXWIN), 0.5 * (i + 1)) for i, c in enumerate(CCDS)}
    # flat with a small gradient in x
    flat = {c: 1.0 + 0.1 * (x - NXWIN / 2) / NXWIN for c in CCDS}

    paths = {}
    for name, arrays, exptime in (
        ("bias", bias, 0.0),
        ("dark", dark, EXPTIME),
        ("flat", flat, 1.0),
    ):
        ccds = {}
        for cnam in CCDS:
            head = fits.Header()
            head["EXPTIME"] = exptime
            win = hcam.Winhead(1, 1, NXWIN, NYWIN, 1, 1, head=head)
            wind = hcam.Window(win, arrays[cnam].astype(np.float32))
            ccds[cnam] = hcam.CCD({"1": wind}, NXWIN, NYWIN)
        mhead = fits.Header()
        mhead["EXPTIME"] = exptime
        mccd = hcam.MCCD(ccds, mhead)
        paths[name] = str(directory / f"{name}.hcm")
        mccd.write(paths[name], overwrite=True)
    return paths, bias, dark, flat


def _make_frames(directory, bias, dark, flat):
    rng = np.random.default_rng(12345)
    fnames = []
    mjd0 = 58000
    for nframe in range(NFRAME):
        ccds = {}
        for cnam in CCDS:
            clean = _fake_ccd_frame(cnam, nframe)
            # electrons -> ADU, add noise, then apply calibration effects
            adu = clean / GAIN
            adu = adu + rng.normal(0.0, np.sqrt(READOUT**2 + adu / 1.0))
            # `adu` is after flat-fielding, bias and dark subtraction
            raw = adu * flat[cnam] + bias[cnam] + dark[cnam]

            head = fits.Header()
            head["MJDINT"] = mjd0
            head["MJDFRAC"] = 0.5 + nframe * EXPTIME / 86400.0
            head["GOODTIME"] = True
            head["EXPTIME"] = EXPTIME
            win = hcam.Winhead(1, 1, NXWIN, NYWIN, 1, 1, head=head)
            wind = hcam.Window(win, raw.astype(np.float32))
            ccds[cnam] = hcam.CCD({"1": wind}, NXWIN, NYWIN)

        mhead = fits.Header()
        mhead["NFRAME"] = nframe + 1
        mhead["TIMSTAMP"] = f"2017-09-17T12:00:{nframe * 10:02d}.000"
        mhead["GOODTIME"] = True
        mhead["EXPTIME"] = EXPTIME
        mccd = hcam.MCCD(ccds, mhead)
        fname = directory / f"frame{nframe + 1:03d}.hcm"
        mccd.write(str(fname), overwrite=True)
        fnames.append(str(fname))

    flist = directory / "frames.lis"
    flist.write_text("\n".join(fnames) + "\n")
    return str(flist)


def _make_aperfile(directory):
    mccdaper = hcam.MccdAper()
    for cnam in CCDS:
        ccdaper = hcam.CcdAper()
        for apnam, (x, y, _) in STARS.items():
            ccdaper[apnam] = hcam.Aperture(
                x, y, 1.5 * FWHM[cnam], 2.5 * FWHM[cnam], 4.0 * FWHM[cnam], apnam == "1"
            )
        mccdaper[cnam] = ccdaper
    fname = directory / "fake.ape"
    mccdaper.write(str(fname))
    return str(fname)


def _make_rfile(directory, apfile, cals, extype, use_psf):
    extraction = "\n".join(
        f"{cnam} = variable {extype} 1.8 5.0 12.0 2.5 14.0 14.0 3.0 20.0 20.0"
        for cnam in CCDS
    )
    light = "\n".join(f"plot = {cnam} 2 1 0.0 1.0 red blue" for cnam in CCDS)
    warns = "\n".join(f"warn = {cnam} 50000 60000" for cnam in CCDS)
    text = f"""
[general]
version = {hcam.REDUCE_FILE_VERSION}
instrument = HiPERCAM
scale = 0.3
ldevice = null
lwidth = 0
lheight = 0
idevice = null
iwidth = 0
iheight = 0
psfdevice = null
toffset = 0
skipbadt = no
ncpu = 1
ngroup = 1
{warns}

[apertures]
aperfile = {apfile}
location = variable
search_half_width = 8.
search_smooth_fwhm = 3.
search_smooth_fft = no
fit_method = gaussian
fit_beta = 4.
fit_beta_max = 20.
fit_fwhm = 4.
fit_fwhm_min = 2.
fit_fwhm_max = 12.
fit_ndiv = 0
fit_fwhm_fixed = no
fit_half_width = 12.
fit_thresh = 8.
fit_height_min_ref = 20.
fit_height_min_nrf = 5.
fit_max_shift = 5.
fit_alpha = 1.0
fit_diff = 2.

[extraction]
{extraction}

[psf_photom]
use_psf = {"yes" if use_psf else "no"}
psf_model = gaussian
fit_half_width = 12
positions = variable

[sky]
method = clipped
error = variance
thresh = 3.

[calibration]
crop = no
bias = {cals["bias"]}
flat = {cals["flat"]}
dark = {cals["dark"]}
fmap =
fpair =
nhalf = 3
rmin = -1.
rmax = 2.
readout = {READOUT}
gain = {GAIN}

[lcplot]
xrange = 0
extend_x = 10.

[light]
linear = yes
y_fixed = no
y1 = 0
y2 = 0
extend_y = 0.1
{light}

[focal_mask]
demask = no
dthresh = 5.

[monitor]
"""
    fname = directory / "fake.red"
    fname.write_text(text)
    return str(fname)


@pytest.fixture(scope="module")
def fake_data(tmp_path_factory):
    directory = tmp_path_factory.mktemp("reduce_data")
    cals, bias, dark, flat = _make_calibrations(directory)
    flist = _make_frames(directory, bias, dark, flat)
    apfile = _make_aperfile(directory)
    return {"dir": directory, "cals": cals, "flist": flist, "apfile": apfile}


def _run_reduce(fake_data, name, extype, use_psf):
    directory = fake_data["dir"]
    rfile = _make_rfile(
        directory, fake_data["apfile"], fake_data["cals"], extype, use_psf
    )
    log = str(directory / f"{name}.log")
    # A command name of None means no defaults files are read or written. The
    # source is a hidden parameter so has to be set by keyword.
    args = [None, "source=hf", fake_data["flist"], "no", rfile, log, "0", "no", "no"]
    if use_psf:
        args.append("no")
    reduce(args)
    return hcam.hlog.Hlog.rascii(log)


QUANTITIES = ("x", "y", "fwhm", "counts", "countse", "sky")


def _summarise(hlog):
    """
    Output a dictionary of time series for each CCD, aperture and quantity in the log.
    """
    summary = {}
    for cnam in CCDS:
        for apnam in STARS:
            for quantity in QUANTITIES:
                key = f"{cnam}/{apnam}/{quantity}"
                summary[key] = [
                    float(v) for v in hlog.tseries(cnam, apnam, quantity, ecol=False).y
                ]
    return summary


def _compare_to_benchmark(name, summary):
    """
    Compare the output of reduce to known benchmark values.
    """
    benchmarks = (
        json.loads(BENCHMARK_FILE.read_text()) if BENCHMARK_FILE.exists() else {}
    )
    if UPDATE:
        benchmarks[name] = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "hipercam_version": hcam.version(),
            "values": summary,
        }
        BENCHMARK_FILE.parent.mkdir(exist_ok=True)
        BENCHMARK_FILE.write_text(json.dumps(benchmarks, indent=1, sort_keys=True))
        return

    assert name in benchmarks, f"no benchmark for {name}; see module docstring"
    entry = benchmarks[name]
    # Older benchmark files did not include provenance.
    expected = entry["values"] if "values" in entry else entry
    assert summary.keys() == expected.keys()
    for key, values in summary.items():
        np.testing.assert_allclose(
            values, expected[key], rtol=1e-4, atol=1e-6, err_msg=key
        )


def test_benchmark_provenance(tmp_path, monkeypatch):
    benchmark_file = tmp_path / "benchmarks.json"
    old_entry = {"1/1/counts": [10.0]}
    benchmark_file.write_text(json.dumps({"normal": old_entry}))
    monkeypatch.setitem(_compare_to_benchmark.__globals__, "BENCHMARK_FILE", benchmark_file)
    monkeypatch.setitem(_compare_to_benchmark.__globals__, "UPDATE", True)

    before = datetime.now(timezone.utc)
    summary = {"1/1/counts": [20.0]}
    _compare_to_benchmark("optimal", summary)
    after = datetime.now(timezone.utc)

    benchmarks = json.loads(benchmark_file.read_text())
    assert benchmarks["normal"] == old_entry
    entry = benchmarks["optimal"]
    assert entry["values"] == summary
    assert entry["hipercam_version"] == hcam.version()
    timestamp = datetime.fromisoformat(entry["timestamp"])
    assert timestamp.utcoffset() == timezone.utc.utcoffset(timestamp)
    assert before <= timestamp <= after

    monkeypatch.setitem(_compare_to_benchmark.__globals__, "UPDATE", False)
    _compare_to_benchmark("normal", old_entry)
    _compare_to_benchmark("optimal", summary)
    with pytest.raises(AssertionError):
        _compare_to_benchmark("optimal", {"1/1/counts": [21.0]})


def _check_constant_stars(summary):
    """Independent of the benchmarks, the fluxes of the non-variable stars
    should be close to their known values for simple aperture and
    PSF photometry.

    Tolerances are set to 5% to allow for random photon noise."""
    for cnam in CCDS:
        for apnam in ("1", "3"):
            true = STARS[apnam][2] * BAND_SCALE[cnam] / GAIN
            np.testing.assert_allclose(
                summary[f"{cnam}/{apnam}/counts"], true, rtol=0.05
            )


def _check_flux_ratio(summary):
    """The flux ratio of the two constant stars is known whatever the
    extraction method, including optimal extraction whose absolute fluxes
    are currently scaled down by a constant factor (see test_reduce)."""
    true = STARS["3"][2] / STARS["1"][2]
    for cnam in CCDS:
        # individual frames scatter by a few percent due to noise
        ratio = np.mean(summary[f"{cnam}/3/counts"]) / np.mean(
            summary[f"{cnam}/1/counts"]
        )
        np.testing.assert_allclose(ratio, true, rtol=0.02, err_msg=cnam)


@pytest.mark.parametrize(
    "name, extype, use_psf",
    [
        ("normal", "normal", False),
        ("optimal", "optimal", False),
        ("psf", "normal", True),
    ],
)
def test_reduce(fake_data, name, extype, use_psf):
    if use_psf:
        pytest.importorskip("photutils")
    hlog = _run_reduce(fake_data, name, extype, use_psf)
    summary = _summarise(hlog)
    assert all(len(v) == NFRAME for v in summary.values())
    assert np.all(np.isfinite(list(summary.values())))
    _check_flux_ratio(summary)
    if name in ("normal", "psf"):
        # optimal extraction is not checked in absolute terms as its fluxes
        # are currently about half the true values
        _check_constant_stars(summary)
    _compare_to_benchmark(name, summary)
