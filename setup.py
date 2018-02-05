"""
hipercam setup file
"""

import os
from setuptools import setup, find_packages
from setuptools.extension import Extension

# To use a consistent encoding
from codecs import open

# need for Cython
import numpy as np
from Cython.Build import cythonize

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# cython support routine
extension = [
    Extension("hipercam.support",
              [os.path.join('hipercam','support.pyx')],
              libraries=["m"],
              include_dirs=[np.get_include()],
              extra_compile_args=["-fno-strict-aliasing"],),
]

setup(
    name='hipercam',

    # Versions should comply with PEP440. Here
    # we use a version generated automatically
    # from git.
    use_scm_version=True,
    setup_requires=['setuptools_scm'],

    description='hipercam',
    long_description=long_description,

    # The project's main homepage.
    url='http://www.astro.warwick.ac.uk',

    # Author details
    author='Tom Marsh',
    author_email='t.r.marsh@warwick.ac.uk',

    # Choose your license
    license='BSD',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Astronomers',
        'Topic :: Astronomy :: Photometric reduction',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
    ],

    # What does your project relate to?
    keywords='astronomy photometry reduction',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    # extension modules
    ext_modules = cythonize(extension),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['sep','numpy','astropy','matplotlib'],

    # need numpy version 1.12

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'sample': ['package_data.dat'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts' : [
            'add=hipercam.scripts.arith:add',
            'div=hipercam.scripts.arith:div',
            'mul=hipercam.scripts.arith:mul',
            'sub=hipercam.scripts.arith:sub',
            'aligntool=hipercam.scripts.aligntool:aligntool',
            'averun=hipercam.scripts.averun:averun',
            'cadd=hipercam.scripts.carith:cadd',
            'cdiv=hipercam.scripts.carith:cdiv',
            'cmul=hipercam.scripts.carith:cmul',
            'csub=hipercam.scripts.carith:csub',
            'combine=hipercam.scripts.combine:combine',
            'digest=hipercam.scripts.digest:digest',
            'grab=hipercam.scripts.grab:grab',
            'genred=hipercam.scripts.genred:genred',
            'hist=hipercam.scripts.hist:hist',
            'hlogger=hipercam.scripts.hlogger:hlogger',
            'hls=hipercam.scripts.hls:hls',
            'hplot=hipercam.scripts.hplot:hplot',
            'makedata=hipercam.scripts.makestuff:makedata',
            'makebias=hipercam.scripts.makebias:makebias',
            'makefield=hipercam.scripts.makestuff:makefield',
            'makeflat=hipercam.scripts.makeflat:makeflat',
            'plog=hipercam.scripts.plog:plog',
            'reduce=hipercam.scripts.reduce:reduce',
            'rtplot=hipercam.scripts.rtplot:rtplot',
            'setaper=hipercam.scripts.setaper:setaper',
            'stats=hipercam.scripts.stats:stats',
            'times=hipercam.scripts.times:times',
        ],
    },

    # tests
    test_suite = 'nose.collector',
    tests_require = ['nose','numpy','astropy'],

)
