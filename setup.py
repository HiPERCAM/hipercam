"""
Minimal setup.py for Cython and pybind11 extension support.
All other metadata is in pyproject.toml.
"""

import os

# need for Cython and pybind11
import numpy as np
from Cython.Build import cythonize
from setuptools import setup
from setuptools.extension import Extension
import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext


# cython support routine
cython_extensions = [
    Extension(
        "hipercam.support",
        [os.path.join("hipercam", "support.pyx")],
        libraries=["m"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-fno-strict-aliasing"],
    ),
]

# pybind11 extension for profile fitting
pybind11_extensions = [
    Pybind11Extension(
        "hipercam.fitting_cpp",
        ["hipercam/fitting_cpp.cpp"],
        include_dirs=[
            np.get_include(),
            pybind11.get_include(),
            pybind11.get_include(user=True),
        ],
        language="c++",
        extra_compile_args=["-std=c++11", "-O3", "-ffast-math"],
    ),
]

setup(
    ext_modules=cythonize(cython_extensions) + pybind11_extensions,
    cmdclass={"build_ext": build_ext},
)
