# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.extern import six

import numpy as np
from hipercam import Window, Windata

win = Window(1,2,300,90,2,3)

assert(win.llx == 1)
assert(win.lly == 2)
assert(win.nx == 300)
assert(win.ny == 90)
assert(win.xbin == 2)
assert(win.ybin == 3)

data = np.ones((90,300))
wind = Windata(win, data)

assert(wind.data[0,0] == 1.)


