"""
Scripts sub-module of HiPERCAM contains all the commands used from
the terminal. They are all implemented as functions for automatic
inclusion in the documentation and for portability
"""

from .arith import add, div, mul, sub
from .atanalysis import atanalysis
from .atbytes import atbytes
from .averun import averun
from .calsearch import calsearch
from .carith import cadd, cdiv, cmul, csub
from .combine import combine
from .exploss import exploss
from .filtid import filtid
from .fits2hcm import fits2hcm
from .flagcloud import flagcloud
from .ftargets import ftargets
from .genred import genred
from .grab import grab
from .harchive import harchive
from .hfilter import hfilter
from .hinfo import hinfo
from .hist import hist
from .hlog2col import hlog2col
from .hlog2fits import hlog2fits
from .hlogger import hlogger
from .hls import hls
from .hmeta import hmeta
from .hpackage import hpackage
from .hplot import hplot
from .joinup import joinup
from .jtrawl import jtrawl
from .logsearch import logsearch
from .ltimes import ltimes
from .ltrans import ltrans
from .makebias import makebias
from .makedark import makedark
from .makeflat import makeflat
from .makefringe import makefringe
from .makemovie import makemovie
from .makestuff import makefield, makemccd
from .mstats import mstats
from .ncal import ncal
from .nrtplot import nrtplot
from .pbands import pbands
from .pfolder import pfolder
from .plog import plog
from .redanal import redanal
from .redplt import redplt
from .reduce import reduce
from .rtplot import rtplot
from .rupdate import rupdate
from .setaper import setaper
from .setdefect import setdefect
from .setfringe import setfringe
from .shiftadd import shiftadd
from .splice import splice
from .stats import stats
from .tanalysis import tanalysis
from .tbytes import tbytes
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
    "hpackage",
    "hplot",
    "joinup",
    "ltimes",
    "ltrans",
    "makebias",
    "makedata",
    "makefield",
    "makeflat",
    "makefringe",
    "makemovie",
    "mstats",
    "mul",
    "ncal",
    "pfolder",
    "plog",
    "redanal",
    "reduce",
    "register",
    "rtplot",
    "nrtplot",
    "rupdate",
    "setaper",
    "setdefect",
    "setfringe",
    "shiftadd",
    "splice",
    "stats",
    "sub",
    "tanalysis",
    "tbytes",
    "uls",
]

try:
    # allow this one to fail
    from .aligntool import aligntool

    __all__.append("aligntool")
except:
    pass
