# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines a class to represent a photometric apertures.
"""

import numpy as np
import json
from collections import OrderedDict
from .core import *
from .group import *

__all__ = ('Aperture',)

class Aperture(object):

    """Class representing a photometric aperture for measuring the flux of a
    star. At its most basic this consists of 3 circles representing the
    object, radius `rtarg`, and the sky annulus between radii `r2` and `r3`,
    centered at a specific position, but there are several extra
    characteristics in addition. These are:

    Reference apertures: these indicate (typically) bright stars that should
    be easy to locate. The idea is that when re-positioning apertures, you
    might want to initially do the brightest stars, and then see what the
    overall x,y shift is before attempting fainter ones. This status is
    indicated with a logical flag.

    Linked apertures: some targets are nearly impossible to register, or may
    fade so much that they cannot be detected. Such targets can be "linked" to
    others in the sense that they are offset from them.

    Sky masks: these are fixed circles offset from the aperture in question
    indicating pixels to ignore when estimating the sky background.

    Star masks: these are circles offset from aperture indicating pixels to
    include when summing the object flux.

    """

    def __init__(self, x, y, rtarg, rsky1, rsky2, ref, link=None, mask=[], extra=[]):
        """
        Constructor. Arguments::

        x     : (float)
            X position of centre of aperture, or the X offset to apply
            if the aperture is linked from another.

        y     : (float)
            Y position of centre of aperture, or the Y offset to apply
            if the aperture is linked from another.

        rtarg : (float)
            Radius (unbinned pixels) of target aperture

        rsky1 : (float)
            Inner radius (unbinned pixels) of sky annulus

        rsky2 : (float)
            Outer radius (unbinned pixels) of sky annulus

        ref   : (bool)
            True/False to indicate whether this is a reference
            aperture meaning that its position will be re-determined
            before non-reference apertures to provide a shift.

        link  : (Aperture / None)
            If set, this is an Aperture that *this* Aperture is linked from.
            Set to `None` for no link.

        mask  : (list of 3 element tuples)
            Each tuple in the list consists of an x,y offset and a radius in
            unbinned pixels. These are used to mask nearby areas when
            determining the sky value (e.g. to exclude stars)

        extra : (list of 3 element tuples)
            Similar to `mask`, each tuple in the list consists of an x,y
            offset and a radius in unbinned pixels. These however are used
            when summing the total flux and allow one to add in extra regions,
            typically blended stars that cannot be allowed to straddle the
            edge of the inner target aperture.

        Notes: normal practice would be to set link, mask, extra later, having created
        the Aperture
        """

        self.x = x
        self.y = y
        self.rtarg = rtarg
        self.rsky1 = rsky1
        self.rsky2 = rsky2
        self.ref = ref
        self._link = link
        self._mask = mask
        self._extra = extra

    def __repr__(self):
        return 'Aperture(x={!r}, y={!r}, rtarg={!r}, rsky1={!r}, rsky2={!r}, ref={!r}, link={!r}, mask={!r}, extra={!r})'.format(
            self.x,self.y,self.rtarg,self.rsky1,self.rsky2,self.ref,self._link,self._mask,self._extra
            )

    def add_mask(self, xoff, yoff, radius):
        """Adds a mask to the :class:Aperture"""
        self._mask.append((xoff,yoff,radius))

    def add_extra(self, xoff, yoff, radius):
        """Adds a mask to the :class:Aperture"""
        self._extra.append((xoff,yoff,radius))

    def set_link(self, aperture):
        """Links this :class:Aperture to another"""
        self._link = aperture

    def break_link(self):
        """Cancels any link to another :class:Aperture"""
        self._link = None

    def toJson(self, fp):
        """Write in JSON-format to file pointer fp"""
        json.dump(self, fp, cls=_encoder, indent=2)

    @classmethod
    def fromJson(cls, fp):
        """Read from JSON-format file pointer fp"""
        return json.load(fp, cls=_decoder)


# classes to support JSON serialisation of Aperture objects
class _encoder(json.JSONEncoder):
    def default(self, aperture):
        return OrderedDict((
                ('Comment', 'This is an Aperture object'),
                ('x', aperture.x),
                ('y', aperture.y),
                ('rtarg', aperture.rtarg),
                ('rsky1', aperture.rsky1),
                ('rsky2', aperture.rsky2),
                ('link', aperture._link),
                ('mask', aperture._mask),
                ('extra', aperture._extra),
                ))

class _decoder(json.JSONDecoder):

    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, jobj):
        if 'rtarg' in jobj and 'rsky1' in jobj and 'rsky2' in jobj and 'link' in jobj:
            return Aperture(
                jobj['x'],jobj['y'],jobj['rtarg'],jobj['rsky1'],jobj['rsky2'],
                jobj['link'],jobj['mask'],jobj['extra']
            )
        else:
            return jobj


class CcdAper(Group):
    """Class representing all the :class:Apertures for a single CCD.
    Normal usage is to create an empty one and then add
    """

    def __init__(self, aps=Group(None)):
        """Constructs a :class:`CcdAper`.

        Arguments::

          aps : (Group)
              Group of :class:`Aperture` objects
        """
        super().__init__(type(Aperture), aps)
