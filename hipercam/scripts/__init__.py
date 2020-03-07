"""
Scripts sub-module of HiPERCAM contains all the commands used from
the terminal. They are all implemented as functions for automatic
inclusion in the documentation and for portability
"""

from .arith import add, div, mul, sub
from .averun import averun
from .carith import cadd, cdiv, cmul, csub
from .combine import combine
from .flagcloud import flagcloud
from .grab import grab
from .genred import genred
from .hist import hist
from .hinfo import hinfo
from .hfilter import hfilter
from .hlogger import hlogger
from .hls import hls
from .hplot import hplot
from .makebias import makebias
from .makedark import makedark
from .makeflat import makeflat
from .makestuff import makemccd, makefield
from .plog import plog
from .reduce import reduce
from .register import register
from .rtplot import rtplot
from .rupdate import rupdate
from .setaper import setaper
from .setdefect import setdefect
from .stats import stats
from .times import times
from .mstats import mstats
from .pfolder import pfolder
from .uls import uls
from .fits2hcm import fits2hcm
from .hlog2fits import hlog2fits
from .redanal import redanal
from .splice import splice

__all__ = [ \
            'add', 'averun',
            'cadd', 'cdiv', 'cmul', 'combine', 'csub',
            'div',
            'genred', 'grab',
            'hist', 'hfilter', 'hlog2fits', 'hlogger', 'hls', 'hplot',
            'makebias', 'makedata', 'makefield', 'makeflat', 'mstats', 'mul',
            'pfolder', 'plog',
            'redanal', 'reduce', 'register', 'rtplot', 'rupdate',
            'setaper', 'setdefect', 'splice', 'stats', 'sub',
            'times',
            'uls',
        ]

try:
    # allow this one to fail
    from .aligntool import aligntool
    __all__.append('aligntool')
except:
    pass

try:
    # optional dependency on photutils, so allow failure
    from .psf_reduce import psf_reduce
    __all__.append('psf_reduce')
except:
    pass

