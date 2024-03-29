# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines classes to represent photometric apertures.

The classes support JSON-style serialisation to allow apertures to be
saved to disk in a fairly easily read and editable format.

"""

import numpy as np
import json
from collections import OrderedDict
from .core import *
from .group import *

__all__ = ("Aperture", "CcdAper", "MccdAper")


class Aperture:

    """Represents an aperture for astronomical photometry

    Essentially this consists of 3 circles representing the object,
    radius `rtarg`, and the sky annulus between radii `r2` and `r3`,
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

    COMPO apertures: comparison stars injected onto the science image by
    COMPO have a difference world coordinate system to the main science image.
    Their motion on the CCD is flipped, so that an (x, y) shift in the stars
    on the science CCD produces a (-x, -y) shift in the COMPO stars. This
    status is indicated with a logical flag.

    Sky masks: ('mask) these are fixed circles offset from the aperture in
    question indicating pixels to ignore when estimating the sky background.

    Star masks: ('extra') these are circles offset from aperture indicating
    pixels to include when summing the object flux. This is to combat problems
    caused by blended objects. These *also* act as sky masks so there should
    be no need to also mask such object.

    Parameters
    ----------
    x : float
        X position of centre of aperture, or the X offset to apply
        if the aperture is linked from another.
    y : float
        Y position of centre of aperture, or the Y offset to apply
        if the aperture is linked from another.
    rtarg : float
        Radius (unbinned pixels) of target aperture
    rsky1 : float
        Inner radius (unbinned pixels) of sky annulus
    rsky2 : float
        Outer radius (unbinned pixels) of sky annulus
    ref : bool
        True/False to indicate whether this is a reference
        aperture meaning that its position will be re-determined
        before non-reference apertures to provide a shift.
    compo : bool
        True/False to indicate whether this is a compo
        aperture meaning that its shifts are inverted w.r.t non-compo
        apertures.
    mask : list of 3 element tuples
        Each tuple in the list consists of an x,y offset and a radius in
        unbinned pixels. These are used to mask nearby areas when
        determining the sky value (e.g. to exclude stars)
    extra : list of 2 element tuples
        Similar to `mask`, but each tuple in the list consists only of an
        x,y offset. These however are used as centres of additional target
        apertures to allow blended stars to be include in the total flux.
        They are given the same radius as the target aperture. They are also
        used to exclude sky pixels.
    link : str
        If != '', this is a string label for another :class:`Aperture` that
        *this* :class:`Aperture` is linked from. The idea is that this label
        can be used to lookup the :class:`Aperture`.

     .. note::

        Normal practice would be to set link, mask, extra later,
        having created the Aperture. Attributes of the same name as
        all the arguments are defined.  We copy the mask and extra
        apertures to avoid propagating references.

    """

    def __init__(
        self, x, y, rtarg, rsky1, rsky2, ref, compo=False, mask=[], extra=[], link=""
    ):
        self.x = x
        self.y = y
        self.rtarg = rtarg
        self.rsky1 = rsky1
        self.rsky2 = rsky2
        self.ref = ref
        self.compo = compo
        self.mask = mask.copy()
        self.extra = extra.copy()
        self.link = link

    def __repr__(self):
        return "Aperture(x={!r}, y={!r}, rtarg={!r}, rsky1={!r}, rsky2={!r}, ref={!r}, compo={!r}, mask={!r}, extra={!r}, link={!r})".format(
            self.x,
            self.y,
            self.rtarg,
            self.rsky1,
            self.rsky2,
            self.ref,
            self.compo,
            self.mask,
            self.extra,
            self.link,
        )

    def copy(self, memo=None):
        """Returns with a copy of the Aperture"""
        return Aperture(
            self.x,
            self.y,
            self.rtarg,
            self.rsky1,
            self.rsky2,
            self.ref,
            self.compo,
            self.mask.copy(),
            self.extra.copy(),
            self.link,
        )

    def add_mask(self, xoff, yoff, radius):
        """Adds a mask to the :class:Aperture"""
        self.mask.append((xoff, yoff, radius))

    def add_extra(self, xoff, yoff):
        """Adds a mask to the :class:Aperture"""
        self.extra.append((xoff, yoff))

    def set_link(self, aplabel):
        """Links this :class:Aperture to a lookup label for another"""
        self.link = aplabel

    def break_link(self):
        """Cancels any link to another :class:Aperture"""
        self.link = ""

    @property
    def linked(self):
        """Returns True if the :class:Aperture is linked to another"""
        return self.link != ""

    def check(self):
        """Run a few checks on an :class:Aperture. Raises a ValueError if there are
        problems.

        """
        if self.rtarg <= 0:
            raise ValueError(
                "Aperture = {!r}\nTarget aperture radius = "
                "{:.2f} <= 0".format(self, self.rtarg)
            )

        elif self.rsky1 > self.rsky2:
            raise ValueError(
                "Aperture = {!r}\nInner sky aperture radius "
                "(={:.2f}) > outer radius (={:.2f})".format(
                    self, self.rsky1, self.rsky2
                )
            )

        elif (
            not isinstance(self.link, str)
            or not isinstance(self.mask, list)
            or not isinstance(self.ref, bool)
            or not isinstance(self.compo, bool)
        ):
            raise ValueError(
                "Aperture = {!r}\nOne or more of link, mask, ref or compo "
                "has the wrong type".format(self)
            )

    def write(self, fname):
        """Dumps Aperture in JSON format to a file called fname"""
        with open(fname, "w") as fp:
            json.dump(self, cls=_Encoder, indent=2)

    def toString(self):
        """Returns Aperture as a JSON-type string"""
        return json.dumps(self, cls=_Encoder, indent=2)

    @classmethod
    def read(cls, fname):
        """Read from JSON-format file fname"""
        with open(fname) as fp:
            aper = json.load(fp, cls=_Decoder)
        aper.check()
        return aper


class CcdAper(Group):
    """Class representing all the :class:Apertures for a single CCD.
    Normal usage is to create an empty one and then add apertures via
    the usual mechanism for updating dictionaries, i.e. ccdap[label] =
    aperture.
    """

    def __init__(self, aps=Group(Aperture)):
        """Constructs a :class:`CcdAper`.

        Arguments::

          aps : (Group)
              Group of :class:`Aperture` objects
        """
        super().__init__(Aperture, aps)

    def __repr__(self):
        return "{:s}(aps={:s})".format(self.__class__.__name__, super().__repr__())

    def check(self):
        """Checks for problems with links"""
        for apnam, aper in self.items():
            if aper.linked:
                if aper.link not in self:
                    raise ValueError(
                        "Aperture = {!r} links to anon-existent aperture".format(self)
                    )
                elif self[aper.link].linked:
                    raise ValueError(
                        "Aperture = {!r} is linked to an aperture which is itself linked".format(
                            self
                        )
                    )

    def write(self, fname):
        """Dumps ccdAper in JSON format to a file called fname"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ["hipercam.CcdAper"] + list(self.items)
        with open(fname, "w") as fp:
            json.dump(listify, fp, cls=_Encoder, indent=2)

    def copy(self, memo=None):
        return CcdAper(super().copy(memo))


class MccdAper(Group):
    """Class representing all the :class:Apertures for multiple CCDs.
    Normal usage is to create an empty one and then add apertures via
    the usual mechanism for updating dictionaries, e.g.

      >> mccdap = MccdAper()
      >> mccdap['ccd1'] = CcdAper()
      >> mccdap['ccd2'] = CcdAper()
      >> mccdap['ccd1']['ap1'] = Aperture(100,200,10,15,125,False)

    etc.
    """

    def __init__(self, aps=Group(CcdAper)):
        """Constructs a :class:`CcdAper`.

        Arguments::

          aps : (Group)
              Group of :class:`CcdAper` objects
        """
        super().__init__(CcdAper, aps)

    def __repr__(self):
        return "{:s}(aps={:s})".format(self.__class__.__name__, super().__repr__())

    def write(self, fname):
        """Dumps MccdAper in JSON format to a file called fname"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ["hipercam.MccdAper"] + list(
            (
                (key, ["hipercam.CcdAper"] + list(val.items()))
                for key, val in self.items()
            )
        )
        with open(fname, "w") as fp:
            json.dump(listify, fp, cls=_Encoder, indent=2)

    def toString(self):
        """Returns MccdAper in JSON format as a string"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ["hipercam.MccdAper"] + list(
            (
                (key, ["hipercam.CcdAper"] + list(val.items()))
                for key, val in self.items()
            )
        )
        return json.dumps(listify, cls=_Encoder, indent=2)

    @classmethod
    def read(cls, fname):
        """Read from JSON-format file fname. Since such files can fairly easily be
        corrupted by injudicious editing, some consistency checks are
        run. File loading does not happen often, so this should not be a
        serious overhead

          fp : a file-like object opened for reading of text

        Returns an MccdAper object.

        """
        with open(fname) as fp:
            obj = json.load(fp, cls=_Decoder)
        listify = [(v1, CcdAper(v2[1:])) for v1, v2 in obj[1:]]
        mccdaper = MccdAper(listify)
        for cnam, ccdaper in mccdaper.items():
            ccdaper.check()
        return mccdaper


# classes to support JSON serialisation of Aperture objects
class _Encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Aperture):
            return OrderedDict(
                (
                    ("Comment", "hipercam.Aperture"),
                    ("x", obj.x),
                    ("y", obj.y),
                    ("rtarg", obj.rtarg),
                    ("rsky1", obj.rsky1),
                    ("rsky2", obj.rsky2),
                    ("ref", obj.ref),
                    ("compo", obj.compo),
                    ("mask", obj.mask),
                    ("extra", obj.extra),
                    ("link", obj.link),
                )
            )

        return super().default(obj)


class _Decoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        # looks out for Aperture objects. Everything else done by default
        if "rtarg" in obj and "rsky1" in obj and "rsky2" in obj and "link" in obj:
            if "compo" not in obj:
                obj["compo"] = False
            return Aperture(
                obj["x"],
                obj["y"],
                obj["rtarg"],
                obj["rsky1"],
                obj["rsky2"],
                obj["ref"],
                obj["compo"],
                obj["mask"],
                obj["extra"],
                obj["link"],
            )

        return obj
