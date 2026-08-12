#!/usr/bin/env python3
"""Test and benchmark the Numpy and C++ avgstd implementations."""

import argparse
import time

import numpy as np

from hipercam import support


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare avgstd timing between the Numpy and C++ implementations"
    )
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--frames", type=int, default=32)
    parser.add_argument("--height", type=int, default=64)
    parser.add_argument("--width", type=int, default=64)
    parser.add_argument("--sigma", type=float, default=3.0)
    parser.add_argument("--seed", type=int, default=12345, help="RNG seed")
    args = parser.parse_args()

    # Build a random data cube for testing
    rng = np.random.default_rng(seed=args.seed)
    cube = rng.normal(loc=0.0, scale=1.0, size=(args.frames, args.height, args.width))
    cube = cube.astype(np.float32)
    print(f"Cube shape: {cube.shape}")

    # Call the functions (this also ensures the C++ extension is loaded and compiled for timing)
    numpy_avg, numpy_std, numpy_num = support._avgstd_numpy(cube, args.sigma)
    cpp_avg, cpp_std, cpp_num = support._avgstd_cpp(cube, args.sigma)

    # Raise an error if the results differ
    if not (
        np.allclose(cpp_avg, numpy_avg)
        and np.allclose(cpp_std, numpy_std)
        and np.array_equal(cpp_num, numpy_num)
    ):
        raise RuntimeError("C++ and Numpy avgstd implementations produced different results")
    else:
        print("C++ and Numpy avgstd implementations produced the same results")

    # Now time the functions and print the results
    numpy_times = []
    for _ in range(args.repeats):
        start = time.perf_counter()
        support._avgstd_numpy(cube, args.sigma)
        numpy_times.append(time.perf_counter() - start)

    cpp_times = []
    for _ in range(args.repeats):
        start = time.perf_counter()
        support._avgstd_cpp(cube, args.sigma)
        cpp_times.append(time.perf_counter() - start)

    numpy_time = float(np.mean(numpy_times))
    cpp_time = float(np.mean(cpp_times))
    print(f"Numpy avgstd mean time over {args.repeats} runs: {numpy_time * 1000:.2f} ms")
    print(f"C++ avgstd mean time over {args.repeats} runs: {cpp_time * 1000:.2f} ms")
    print(f"Speedup: {numpy_time / cpp_time:.2f}x")


if __name__ == "__main__":
    main()
