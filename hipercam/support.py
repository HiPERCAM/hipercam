"""Support routines for image combination and profile evaluation.

The implementation uses a pybind11-backed extension when available and falls
back to a pure-NumPy implementation otherwise.
"""

import numpy as np

from . import _support_cpp

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

    return _support_cpp.avgstd(cube_arr, float(sigma))
