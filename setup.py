"""
Minimal setup.py for Cython extension support.
All other metadata is in pyproject.toml.
"""

import os
import numpy as np
from Cython.Build import cythonize
from setuptools import setup
from setuptools.extension import Extension

# cython support routine
extension = [
    Extension(
        "hipercam.support",
        [os.path.join("hipercam", "support.pyx")],
        libraries=["m"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-fno-strict-aliasing"],
    ),
]

setup(
    ext_modules=cythonize(extension),
)
