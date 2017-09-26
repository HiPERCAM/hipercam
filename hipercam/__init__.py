# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
reduction package for the high-speed camera HiPERCAM

hipercam is a package for the reduction of data from the 5-band multi-window
high-speed CCD camera HiPERCAM. It can access the raw data, display it, slice
it into individual files containing the 5 CCDs, etc.

When using it for the first time, your first port of entry is probably to look
at the 'scripts' subpackage, e.g. pydoc hipercam.scripts
"""

from .core import *
from . import cline
from .group import *
from .window import *
from .ccd import *
from .aperture import *
from .target import *
from . import mpl
from . import pgp
from . import ucam
from . import hcam
from .spooler import *
