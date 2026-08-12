"""Support routines for image combination and profile evaluation.

The implementation uses a pybind11-backed extension when available and falls
back to a pure-NumPy implementation otherwise.
"""

from __future__ import annotations
import importlib

import numpy as np

try:
    _support_cpp = importlib.import_module("._support_cpp", __package__)
    _avgstd_cpp = _support_cpp.avgstd
    SUPPORT_CPP_AVAILABLE = True
except ImportError:
    SUPPORT_CPP_AVAILABLE = False

__all__ = ["avgstd"]


def avgstd(cube, sigma):
    """Compute clipped mean, standard deviation and accepted-frame counts.

    Parameters
    ----------
    cube : numpy.ndarray
        Three-dimensional array of shape ``(nf, ny, nx)`` with ``float32`` data.
    sigma : float
        Rejection threshold for the iterative clipped-mean routine.

    Returns
    -------
    tuple
        ``(avg, std, num)`` arrays with shapes ``(ny, nx)``.
    """

    cube_arr = np.asarray(cube, dtype=np.float32)
    if cube_arr.ndim != 3:
        raise ValueError("cube must be a 3D array")
    if sigma <= 1.0:
        raise ValueError("sigma must be greater than 1")

    if SUPPORT_CPP_AVAILABLE:
        return _avgstd_cpp(cube_arr, float(sigma))

    return _avgstd_numpy(cube_arr, float(sigma))


def _avgstd_numpy(cube_arr, sigma):
    """Compute clipped mean, standard deviation and accepted-frame counts."""
    nf, ny, nx = cube_arr.shape
    avg = np.empty((ny, nx), dtype=np.float32)
    std = np.empty((ny, nx), dtype=np.float32)
    num = np.empty((ny, nx), dtype=np.int32)

    for iy in range(ny):
        for ix in range(nx):
            vals = cube_arr[:, iy, ix].astype(np.float64)
            ok = np.ones(nf, dtype=bool)
            ncur = nf

            while True:
                if ncur <= 0:
                    avg[iy, ix] = 0.0
                    std[iy, ix] = 0.0
                    num[iy, ix] = 0
                    break

                tavg = vals[ok].mean()
                if ncur > 1:
                    tstd = vals[ok].std(ddof=1)
                else:
                    tstd = 0.0

                thresh = sigma * tstd
                new_ok = ok.copy()
                new_ok[np.abs(vals - tavg) > thresh] = False

                if np.all(new_ok == ok):
                    kept = vals[ok]
                    avg[iy, ix] = float(kept.mean())
                    std[iy, ix] = float(kept.std(ddof=1)) if kept.size > 1 else 0.0
                    num[iy, ix] = int(kept.size)
                    break

                ok = new_ok
                ncur = int(np.count_nonzero(ok))

    return avg, std, num
