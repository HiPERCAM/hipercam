# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines a class to represent a photometric apertures.
"""

import numpy as np
from .core import *

class Aperture(object):

    """Class representing a photometric aperture for measuring the flux of a
    star. At its most basic this consists of 3 circles representing the
    object, radius `r1`, and the sky annulus between radii `r2` and `r3`,
    centered at a specific position, but there are several extra
    characteristics in addition. These are:

    Reference apertures: these indicate (typically) bright stars that should
    be easy to locate. The idea is that when re-positioning apertures, you
    might want to initially do the brightest stars, and then see what the
    overall x,y shift is before attempting fainter ones. One indicates this
    status with a logical flag.

    Linked apertures: some targets are nearly impossible to register, or may
    fade so much that they cannot be detected. Such targets can be "linked" to
    others in the sense that they are offset from them.

    Sky masks: these are fixed circles offset from the aperture in question
    indicating pixels to ignore when estimating the sky background.

    Star masks: these are circles offset from aperture indicating pixels to
    include when summing the object flux.

    """

    def __init__(self, x, y, r1, r2, r3, ref, link=None, mask=None, extra=None):
        """
        Constructor. Arguments::

        x  : (float)
            X position of centre of aperture, or the X offset to apply
            if the aperture is linked from another.

        y  : (float)
            Y position of centre of aperture, or the Y offset to apply
            if the aperture is linked from another.

        r1 : (float)
            Radius (unbinned pixels) of target aperture

        r2 : (float)
            Inner radius (unbinned pixels) of sky annulus

        r3 : (float)
            Outer radius (unbinned pixels) of sky annulus

        ref : (bool)
            True/False to indicate whether this is a reference
            aperture meaning that its position will be re-determined
            before non-reference apertures to provide a shift.

        link : (key value)
            If set, this gives a lookup key to find the aperture this is
            linked from. Set to `None` for no link.

        mask : (list of 3 element tuples)
            Each tuple in the list consists of an x,y offset and a radius in
            unbinned pixels. These are used to mask nearby areas when
            determining the sky value (e.g. to exclude stars)

        extra : (list of 3 element tuples)
            Similar to `mask`, each tuple in the list consists of an x,y
            offset and a radius in unbinned pixels. These however are used
            when summing the total flux and allow one to add in extra regions,
            typically blended stars that cannot be allowed to straddle the
            edge of the inner target aperture.
        """

        self.x = x
        self.y = y
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.ref = ref
        self.link = link
        self.mask = mask
        self.mask = extra

        # A few elementary checks, but mostly it is left up to
        # "duck typing"
        if link is not None and not isinstance(link,list):
            raise ValueError('Aperture.__init__: link must be a list')

        if mask is not None and not isinstance(mask,list):
            raise ValueError('Aperture.__init__: mask must be a list')

        if extra is not None and not isinstance(extra,list):
            raise ValueError('Aperture.__init__: extra must be a list')

    def __repr__(self):
        return 'Aperture(x=:{r}, y={:r}, r1={:r}, r2={:r}, r3={:r}, ref={:r}, link={:r}, mask={:r}, extra={:r})'.format(
            self.x,self.y,self.r1,self.r2,self.r3,
            self.ref,self.link,self.extra)
