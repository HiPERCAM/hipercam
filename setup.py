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

    # Versions should comply with PEP440. Here we use a version generated
    # automatically from git.
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
        'Programming Language :: Python :: 3.6',
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

    # List run-time dependencies here.  These will be installed by pip
    # when your project is installed. For an analysis of
    # "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # urllib3 requirement to get over a security warning
    install_requires=[
        'sep', 'numpy', 'astropy', 'matplotlib', 'requests',
        'numba', 'websocket-client', 'fitsio', 'pandas', 'Cython',
        'urllib3>=1.26.5', 'keyring', 'lacosmic'
    ],

    # Makes significant use of f-strings which came in python v3.6
    python_requires='>=3.6',

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
            'aligntool=hipercam.scripts.aligntool:aligntool',
            'atanalysis=hipercam.scripts.atanalysis:atanalysis',
            'atbytes=hipercam.scripts.atbytes:atbytes',
            'averun=hipercam.scripts.averun:averun',
            'cadd=hipercam.scripts.carith:cadd',
            'calsearch=hipercam.scripts.calsearch:calsearch',
            'cdiv=hipercam.scripts.carith:cdiv',
            'cmul=hipercam.scripts.carith:cmul',
            'combine=hipercam.scripts.combine:combine',
            'csub=hipercam.scripts.carith:csub',
            'div=hipercam.scripts.arith:div',
            'exploss=hipercam.scripts.exploss:exploss',
            'filtid=hipercam.scripts.filtid:filtid',
            'fits2hcm=hipercam.scripts.fits2hcm:fits2hcm',
            'flagcloud=hipercam.scripts.flagcloud:flagcloud',
            'ftargets=hipercam.scripts.ftargets:ftargets',
            'genred=hipercam.scripts.genred:genred',
            'grab=hipercam.scripts.grab:grab',
            'harchive=hipercam.scripts.harchive:harchive',
            'hfilter=hipercam.scripts.hfilter:hfilter',
            'hinfo=hipercam.scripts.hinfo:hinfo',
            'hist=hipercam.scripts.hist:hist',
            'hlogger=hipercam.scripts.hlogger:hlogger',
            'hlog2col=hipercam.scripts.hlog2col:hlog2col',
            'hlog2fits=hipercam.scripts.hlog2fits:hlog2fits',
            'hls=hipercam.scripts.hls:hls',
            'hmeta=hipercam.scripts.hmeta:hmeta',
            'hpackage=hipercam.scripts.hpackage:hpackage',
            'hplot=hipercam.scripts.hplot:hplot',
            'joinup=hipercam.scripts.joinup:joinup',
            'jtrawl=hipercam.scripts.jtrawl:jtrawl',
            'logsearch=hipercam.scripts.logsearch:logsearch',
            'ltimes=hipercam.scripts.ltimes:ltimes',
            'ltrans=hipercam.scripts.ltrans:ltrans',
            'makebias=hipercam.scripts.makebias:makebias',
            'makedark=hipercam.scripts.makedark:makedark',
            'makefield=hipercam.scripts.makestuff:makefield',
            'makeflat=hipercam.scripts.makeflat:makeflat',
            'makefringe=hipercam.scripts.makefringe:makefringe',
            'makemccd=hipercam.scripts.makestuff:makemccd',
            'makemovie=hipercam.scripts.makemovie:makemovie',
            'mstats=hipercam.scripts.mstats:mstats',
            'mul=hipercam.scripts.arith:mul',
            'ncal=hipercam.scripts.ncal:ncal',
            'pbands=hipercam.scripts.pbands:pbands',
            'pfolder=hipercam.scripts.pfolder:pfolder',
            'plog=hipercam.scripts.plog:plog',
            'psf_reduce=hipercam.scripts.psf_reduce:psf_reduce',
            'redanal=hipercam.scripts.redanal:redanal',
            'redplt=hipercam.scripts.redplt:redplt',
            'reduce=hipercam.scripts.reduce:reduce',
            'rtplot=hipercam.scripts.rtplot:rtplot',
            'nrtplot=hipercam.scripts.nrtplot:nrtplot',
            'rupdate=hipercam.scripts.rupdate:rupdate',
            'setaper=hipercam.scripts.setaper:setaper',
            'setdefect=hipercam.scripts.setdefect:setdefect',
            'setfringe=hipercam.scripts.setfringe:setfringe',
            'splice=hipercam.scripts.splice:splice',
            'stats=hipercam.scripts.stats:stats',
            'sub=hipercam.scripts.arith:sub',
            'tanalysis=hipercam.scripts.tanalysis:tanalysis',
            'tbytes=hipercam.scripts.tbytes:tbytes',
            'uls=hipercam.scripts.uls:uls',
        ],

    },

    # tests
    test_suite = 'nose.collector',
    tests_require = ['nose','numpy','astropy'],

)
