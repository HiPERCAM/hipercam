#!/usr/bin/env python3
"""Test and benchmark the Numba and C++ fitting functions."""

import argparse
import time

import numpy as np

from hipercam import fitting


def _assert_allclose(name: str, a: np.ndarray, b: np.ndarray, rtol: float, atol: float) -> None:
    """Assert that two arrays are close, raising an error with details if not."""
    if not np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
        diff = np.max(np.abs(a - b))
        rel = np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300))
        raise AssertionError(
            f"{name} mismatch: max_abs={diff:.3e}, max_rel={rel:.3e}, "
            f"rtol={rtol:.1e}, atol={atol:.1e}"
        )


def _make_case(rng: np.random.Generator) -> dict[str, object]:
    """Make a single random test case for fitting functions."""
    # Make a simple 2D grid of x and y coordinates
    ny = int(rng.integers(8, 24))
    nx = int(rng.integers(8, 24))
    x1 = np.linspace(100.25, 100.25 + 0.6 * (nx - 1), nx)
    y1 = np.linspace(200.75, 200.75 + 0.7 * (ny - 1), ny)
    x, y = np.meshgrid(x1, y1)

    return {
        "x": x,
        "y": y,
        "sky": float(rng.uniform(-100.0, 500.0)),
        "height": float(rng.uniform(1.0, 5e4)),
        "xcen": float(rng.uniform(x.min() - 1.0, x.max() + 1.0)),
        "ycen": float(rng.uniform(y.min() - 1.0, y.max() + 1.0)),
        "fwhm": float(rng.uniform(0.5, 8.0)),
        "beta": float(rng.uniform(0.05, 8.0)),
        "xbin": int(rng.integers(1, 5)),
        "ybin": int(rng.integers(1, 5)),
        "ndiv": int(rng.integers(0, 4)),
    }


def _make_cases(ncases: int, seed: int) -> list[dict[str, object]]:
    """Make a list of random test cases for fitting functions."""
    rng = np.random.default_rng(seed)
    return [_make_case(rng) for _ in range(ncases)]


def run_comparison(case_id: int, case: dict[str, object], rtol: float, atol: float) -> None:
    """Run a single test case, comparing C++ and Numba implementations."""
    if not fitting.FITTING_CCP_AVAILABLE:
        raise RuntimeError("hipercam.fitting_cpp is not available. Build/install extension first.")

    # Moffat profile
    results_numba = fitting._moffat_numba(
        case['x'],
        case['y'],
        case['sky'],
        case['height'],
        case['xcen'],
        case['ycen'],
        case['fwhm'],
        case['beta'],
        case['xbin'],
        case['ybin'],
        case['ndiv']
    )
    results_cpp = fitting._moffat_cpp(
        case['x'],
        case['y'],
        case['sky'],
        case['height'],
        case['xcen'],
        case['ycen'],
        case['fwhm'],
        case['beta'],
        case['xbin'],
        case['ybin'],
        case['ndiv']
    )
    _assert_allclose(f"case {case_id} moffat", results_numba, results_cpp, rtol, atol)

    # Moffat derivatives (all flag combinations)
    for comp_dfwhm in (False, True):
        for comp_dbeta in (False, True):
            results_numba = fitting._dmoffat_numba(
                case['x'],
                case['y'],
                case['sky'],
                case['height'],
                case['xcen'],
                case['ycen'],
                case['fwhm'],
                case['beta'],
                case['xbin'],
                case['ybin'],
                case['ndiv'],
                comp_dfwhm,
                comp_dbeta,
            )
            results_cpp = fitting._dmoffat_cpp(
                case['x'],
                case['y'],
                case['sky'],
                case['height'],
                case['xcen'],
                case['ycen'],
                case['fwhm'],
                case['beta'],
                case['xbin'],
                case['ybin'],
                case['ndiv'],
                comp_dfwhm,
                comp_dbeta,
            )
            _assert_allclose(
                f"case {case_id} dmoffat[dfwhm={comp_dfwhm},dbeta={comp_dbeta}]",
                results_numba,
                results_cpp,
                rtol,
                atol,
            )

    # Gaussian profile
    results_numba = fitting._gaussian_numba(
        case['x'],
        case['y'],
        case['sky'],
        case['height'],
        case['xcen'],
        case['ycen'],
        case['fwhm'],
        case['xbin'],
        case['ybin'],
        case['ndiv']
    )
    results_cpp = fitting._gaussian_cpp(
        case['x'],
        case['y'],
        case['sky'],
        case['height'],
        case['xcen'],
        case['ycen'],
        case['fwhm'],
        case['xbin'],
        case['ybin'],
        case['ndiv']
    )
    _assert_allclose(f"case {case_id} gaussian", results_numba, results_cpp, rtol, atol)

    # Gaussian derivatives
    for comp_dfwhm in (False, True):
        results_numba = fitting._dgaussian_numba(
            case['x'],
            case['y'],
            case['sky'],
            case['height'],
            case['xcen'],
            case['ycen'],
            case['fwhm'],
            case['xbin'],
            case['ybin'],
            case['ndiv'],
            comp_dfwhm
        )
        results_cpp = fitting._dgaussian_cpp(
            case['x'],
            case['y'],
            case['sky'],
            case['height'],
            case['xcen'],
            case['ycen'],
            case['fwhm'],
            case['xbin'],
            case['ybin'],
            case['ndiv'],
            comp_dfwhm
        )
        _assert_allclose(
            f"case {case_id} dgaussian[dfwhm={comp_dfwhm}]",
            results_numba,
            results_cpp,
            rtol,
            atol,
            )


def _time_function(
    callable_fn, cases: list[dict[str, object]], nruns: int, *args, **kwargs
) -> np.ndarray:
    """Time a function over multiple runs with given cases and return the elapsed times."""
    times = []
    for _ in range(nruns):
        start = time.perf_counter()
        for case in cases:
            callable_fn(case, *args, **kwargs)
        times.append(time.perf_counter() - start)
    return [time / len(cases) for time in times]  # Return time per case


def _compare_functions(function_numba, function_cpp, cases, nruns, *args, **kwargs):
    """Compare the timing of Numba and C++ functions over multiple runs."""
    times_numba = _time_function(function_numba, cases, nruns, *args, **kwargs)
    times_cpp = _time_function(function_cpp, cases, nruns, *args, **kwargs)
    print(
        f"{np.mean(times_numba) * 1000:5.3f}±{np.std(times_numba) * 1000:5.3f} ms     "
        f"{np.mean(times_cpp) * 1000:5.3f}±{np.std(times_cpp) * 1000:5.3f} ms     "
        f"{np.mean(times_numba) / np.mean(times_cpp):5.2f}x"
    )


def run_benchmarks(cases: list[dict[str, object]], nruns: int) -> None:
    """Run benchmarks comparing C++ and Numba implementations of fitting kernels."""
    print(f"\nBenchmarking {len(cases)} cases over {nruns} runs:")
    print("-" * 85)
    print(f"{'function':32s} {'Numba':18s} {'C++':18s} {'Speedup':8s}")
    print("-" * 85)

    # Run the comparisons
    # Moffat kernel
    def moffat_numba(case):
        return fitting._moffat_numba(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["beta"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
        )
    def moffat_cpp(case):
        return fitting._moffat_cpp(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["beta"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
        )
    print(f"{'moffat':32s}", end=" ")
    _compare_functions(moffat_numba, moffat_cpp, cases, nruns)

    # Gaussian kernel
    def gaussian_numba(case):
        return fitting._gaussian_numba(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
        )
    def gaussian_cpp(case):
        return fitting._gaussian_cpp(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
        )
    print(f"{'gaussian':32s}", end=" ")
    _compare_functions(gaussian_numba, gaussian_cpp, cases, nruns)

    # Moffat derivatives
    def dmoffat_numba(case, comp_dfwhm, comp_dbeta):
        return fitting._dmoffat_numba(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["beta"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
            comp_dfwhm,
            comp_dbeta,
        )
    def dmoffat_cpp(case, comp_dfwhm, comp_dbeta):
        return fitting._dmoffat_cpp(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["beta"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
            comp_dfwhm,
            comp_dbeta,
        )
    for comp_dfwhm in (False, True):
        for comp_dbeta in (False, True):
            print(f"{f'dmoffat[dfwhm={comp_dfwhm},dbeta={comp_dbeta}]':32s}", end=" ")
            _compare_functions(dmoffat_numba, dmoffat_cpp, cases, nruns, comp_dfwhm, comp_dbeta)

    # Gaussian derivatives
    def dgaussian_numba(case, comp_dfwhm):
        return fitting._dgaussian_numba(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
            comp_dfwhm,
        )
    def dgaussian_cpp(case, comp_dfwhm):
        return fitting._dgaussian_cpp(
            case["x"],
            case["y"],
            case["sky"],
            case["height"],
            case["xcen"],
            case["ycen"],
            case["fwhm"],
            case["xbin"],
            case["ybin"],
            case["ndiv"],
            comp_dfwhm,
        )
    for comp_dfwhm in (False, True):
        print(f"{f'dgaussian[dfwhm={comp_dfwhm}]':32s}", end=" ")
        _compare_functions(dgaussian_numba, dgaussian_cpp, cases, nruns, comp_dfwhm)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare fitting timing between the Numpy and C++ implementations"
    )
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--cases", type=int, default=20)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-11)
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    # Make the random test cases
    cases = _make_cases(args.cases, args.seed)

    # Run the equivalence checks for each case
    for case_id, case in enumerate(cases, start=1):
        run_comparison(case_id, case, args.rtol, args.atol)
    print(f"All checks passed: {args.cases} cases, rtol={args.rtol:.1e}, atol={args.atol:.1e}")

    # Run benchmarks (using a new set of cases with a different seed)
    bench_cases = _make_cases(args.cases, args.seed + 1)
    run_benchmarks(bench_cases, args.repeats)


if __name__ == "__main__":
    main()
