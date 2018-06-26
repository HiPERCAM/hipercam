"""
Scripts sub-module of HiPERCAM contains all the commands used from
the terminal. They are all implemented as functions for automatic
inclusion in the documentation and for portability
"""

from .arith import add, div, mul, sub
from .averun import averun
from .carith import cadd, cdiv, cmul, csub
from .combine import combine
from .grab import grab
from .genred import genred
from .hist import hist
from .hfilter import hfilter
from .hlogger import hlogger
from .hls import hls
from .hplot import hplot
from .makebias import makebias
from .makeflat import makeflat
from .makestuff import makedata, makefield
from .plog import plog
from .reduce import reduce
from .rtplot import rtplot
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

__all__ = [ \
            'add', 'averun',
            'cadd', 'cdiv', 'cmul', 'combine', 'csub',
            'div',
            'genred', 'grab',
            'hist', 'hfilter', 'hlog2fits', 'hlogger', 'hls', 'hplot',
            'makebias', 'makedata', 'makefield', 'makeflat', 'mstats', 'mul',
            'pfolder', 'plog',
            'redanal', 'reduce', 'rtplot',
            'setaper', 'setdefect', 'stats', 'sub',
            'times',
            'uls',
        ]

try:
    # allow this one to fail
    from .aligntool import aligntool

    __all__.append('aligntool')
except:
    pass

