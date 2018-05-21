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

__all__ = [ \
            'add', 'averun',
            'cadd', 'cdiv', 'cmul', 'combine', 'csub',
            'div',
            'genred', 'grab',
            'hist', 'hfilter', 'hlogger', 'hls', 'hplot',
            'makebias', 'makedata', 'makefield', 'makeflat', 'mstats', 'mul',
            'pfolder', 'plog',
            'reduce', 'rtplot',
            'setaper', 'setdefect', 'stats', 'sub',
            'times',
        ]

try:
    from .aligntool import aligntool

    # insert after add
    __all__.insert(1,'aligntool')
except:
    pass

