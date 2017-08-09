# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
hipercam is a package for the reduction of data from the 5-band multi-window
high-speed CCD camera HiPERCAM.
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
