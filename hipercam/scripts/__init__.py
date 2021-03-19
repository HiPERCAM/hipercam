"""
Scripts sub-module of HiPERCAM contains all the commands used from
the terminal. They are all implemented as functions for automatic
inclusion in the documentation and for portability
"""

from .arith import add, div, mul, sub
from .atbytes import atbytes
from .atanalysis import atanalysis
from .averun import averun
from .carith import cadd, cdiv, cmul, csub
from .combine import combine
from .filtid import filtid
from .fits2hcm import fits2hcm
from .flagcloud import flagcloud
from .ftargets import ftargets
from .grab import grab
from .genred import genred
from .hist import hist
from .hinfo import hinfo
from .hfilter import hfilter
from .hlogger import hlogger
from .hlog2fits import hlog2fits
from .hlog2col import hlog2col
from .hls import hls
from .hmeta import hmeta
from .hplot import hplot
from .ltimes import ltimes
from .logsearch import logsearch
from .ltrans import ltrans
from .makebias import makebias
from .makedark import makedark
from .makeflat import makeflat
from .makefringe import makefringe
from .makestuff import makemccd, makefield
from .mstats import mstats
from .pfolder import pfolder
from .plog import plog
from .redplt import redplt
from .reduce import reduce
from .rtplot import rtplot
from .rupdate import rupdate
from .setaper import setaper
from .setdefect import setdefect
from .stats import stats
from .redanal import redanal
from .splice import splice
from .tbytes import tbytes
from .tanalysis import tanalysis
from .ulogger import ulogger
from .uls import uls

__all__ = [
    "add",
    "atanalysis",
    "atbytes",
    "averun",
    "cadd",
    "cdiv",
    "cmul",
    "combine",
    "csub",
    "div",
    "genred",
    "grab",
    "hist",
    "hfilter",
    "hlog2fits",
    "hlogger",
    "hls",
    "hplot",
    "ltimes",
    "ltrans",
    "makebias",
    "makedata",
    "makefield",
    "makeflat",
    "mstats",
    "mul",
    "pfolder",
    "plog",
    "redanal",
    "reduce",
    "register",
    "rtplot",
    "rupdate",
    "setaper",
    "setdefect",
    "splice",
    "stats",
    "sub",
    "tanalysis",
    "tbytes",
    "ulogger",
    "uls",
]

try:
    # allow this one to fail
    from .aligntool import aligntool

    __all__.append("aligntool")
except:
    pass

try:
    # optional dependency on photutils, so allow failure
    from .psf_reduce import psf_reduce

    __all__.append("psf_reduce")
except:
    pass
