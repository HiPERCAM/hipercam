"""
Minimal setup.py for pybind11 extension support.
All other metadata is in pyproject.toml.
"""

import os
import sys

import numpy as np
import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# OpenMP support
use_openmp = os.environ.get("HIPERCAM_USE_OPENMP", "1") != "0"
openmp_args = []
openmp_link_args = []
if use_openmp and sys.platform.startswith("linux"):
    openmp_args = ["-fopenmp"]
    openmp_link_args = ["-fopenmp"]

# pybind11 extensions
pybind11_extensions = [
    # Profile fitting functions
    Pybind11Extension(
        name="hipercam._fitting_cpp",
        sources=["hipercam/fitting.cpp"],
        include_dirs=[
            np.get_include(),
            pybind11.get_include(),
            pybind11.get_include(user=True),
        ],
        extra_compile_args=["-std=c++11", "-O3", "-ffast-math", "-march=native"] + openmp_args,
        extra_link_args=openmp_link_args,
    ),
    # Support functions
    Pybind11Extension(
        name="hipercam._support_cpp",
        sources=["hipercam/support.cpp"],
        include_dirs=[
            np.get_include(),
            pybind11.get_include(),
            pybind11.get_include(user=True),
        ],
        extra_compile_args=["-std=c++11", "-O3", "-ffast-math", "-march=native"] + openmp_args,
        extra_link_args=openmp_link_args,
    ),
]

setup(
    ext_modules=pybind11_extensions,
    cmdclass={"build_ext": build_ext},
)
