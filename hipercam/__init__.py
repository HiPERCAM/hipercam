# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""hipercam - a reduction package for the high-speed camera HiPERCAM

:mod:`hipercam` is a package for the reduction of data from the 5-band
multi-window high-speed CCD camera HiPERCAM. It can access the raw data,
display it, slice it into individual files containing the 5 CCDs, etc. It
provides a combination of an API to access HiPERCAM data and high-level
scripts to enact certain commonly-required operations in reduction. Most users
will be most interested in the latter, and if that is the case for you, your
first port of entry is probably to look at the 'scripts' sub-package,
e.g. ``pydoc hipercam.scripts``.

If you are interested in the API, then read on. The code is compartmentalised
into several sub-modules and sub-packages, e.g. :class:`CCD` and :class:`MCCD`
appear in ``hipercam.ccd``, but they are also imported to the top level, so that
e.g.  ``pydoc hipercam.CCD`` and ``pydoc hipercam.ccd.CCD`` report back the same
thing. This is true for most of the sub-modules, so looking through the
top-level docs covers a large part of the API. However, it is rather long and
if you want to specialise to a restricted part of the documentation, then a
command like 'pydoc hipercam.ccd' can prove useful to find a more focussed
view of related software. There are some exceptions though where you will need
to drill deeper, namely the sub-packages like 'ucam' and 'hcam', the
sub-module :mod:`hipercam.cline`, and the plotting modules :mod:`mpl` and 
:mod:`pgp`.

Perhaps the most important class of all is :class:`MCCD` which is the first
place to go if you have a HiPERCAM '.hcm' file you want to do something with
beyond the current routines (although NB '.hcm' files are just FITS files in
thin disguise of a changed extension and you may be able to do your own stuff
with standard tools for FITS files). Example code::

  >> import hipercam as hcam
  >> mccd = hcam.MCCD.read('myfile.hcm')
  >> print(mccd)
  >> print(mccd['1'])
  >> print(mccd['1']['E1'])
  >> print(mccd['1']['E1'][123,45])
  >> mccd['1']['E1'][123,45] = 0
  >> mccd.write('mynewfile.hcm')

which reads an hcm file in, then prints it, the CCD it contains labelled '1'
(a :class:`CCD` object), then the window data from this CCD labelled 'E1' (a
:class:`Window` object) and finally the pixel value at iy=123, ix=45 (C-like
starting from 0). This shows the dictionary-like interface used for
:class:`MCCD` and :class:`CCD` objects alike (they are both sub-classes of a
container class :class:`Group`). It then sets this pixel to zero and writes
out the result to a new file.
"""

from .core import *
from .group import *
from .window import *
from .ccd import *
from .aperture import *
from .target import *
from . import spooler
from . import cline
from . import mpl
from . import pgp
from . import ucam
from . import utils
from . import hcam
from . import hlog
from . import support
from . import fitting

__all__ = core.__all__ + group.__all__ + window.__all__ + ccd.__all__ + \
    aperture.__all__ + target.__all__
