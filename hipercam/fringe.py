# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to specify pairs of points for fringe
measurements. See 'setfringe' for usage.
"""

import numpy as np
import json
from collections import OrderedDict

from .core import *
from .group import *
from . import utils

class FringePair:

    """Class representing a pair of positions for fringe
    measurement. The idea is that ratio of the difference
    in intensity at the two positions compared between two
    frames is used to scale a reference fringe image before
    subtracting from data.
    """

    def __init__(self, x1, y1, x2, y2):
        """Constructor. Arguments::

           severity     : Severity
              indicator of the severity of a defect
        """
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2

    def copy(self, memo=None):
        """Returns with a copy of the Aperture"""
        return FringePair(
            self.x1,
            self.y1,
            self.x2,
            self.y2
        )

    def dist(self, x, y):
        ex = self.x2-self.x1
        ey = self.y2-self.y1
        ax = x-self.x1
        ay = y-self.y1
        if ex == 0. or ey == 0.:
            d = np.sqrt(ax**2+ay**2)
        else:
            ell = np.sqrt(ex**2+ey**2)
            ex /= ell
            ey /= ell
            mu = ex*ax+ey*ay
            mu = max(0,min(mu,1))
            xm = self.x1+mu*ex
            ym = self.y1+mu*ey
            d = np.sqrt((x-xm)**2+(y-ym)**2)
        return d

    def __repr__(self):
        return \
            f"FringePair(x1={self.x1:.2f}, y1={self.y1:.2f}," \
            f" x2={self.x2:.2f}, y2={self.y2:.2f)"

class CcdFringePair(Group):
    """Class representing all the :class:`FringePair`s for a single CCD.
    Normal usage is to create an empty one and then add FringePairs
    via the usual mechanism for updating dictionaries,
    i.e. ccddef[label] = fpair

    """

    def __init__(self, fpairs=Group(FringePair)):
        """Constructs a :class:`CcdFringePair`.

        Arguments::

          fpairs : Group
              Group of :class:`FringePair` objects
        """
        super().__init__(FringePair, defs)

    def __repr__(self):
        return "{:s}(fpairs={:s})".format(self.__class__.__name__, super().__repr__())

    def write(self, fname):
        """Dumps ccdAper in JSON format to a file called fname"""

        # dumps as list to retain order through default iterator
        # encoding that buggers things otherwise
        listify = ["hipercam.CcdDefect"] + list(self.items)
        with open(fname, "w") as fp:
            json.dump(listify, fp, cls=_Encoder, indent=2)

    def copy(self, memo=None):
        return CcdDefect(super().copy(memo))


class MccdFringePair(Group):
    """Class representing all the :class:FringePairs for multiple CCDs.
    Normal usage is to create an empty one and then add apertures via
    the usual mechanism for updating dictionaries, e.g.

      >> mccddef = MccdFringePair()
      >> mccddef['ccd1'] = CcdFringePair()
      >> mccddef['ccd2'] = CcdFringePair()
      >> mccddef['ccd1']['ap1'] = FringePair(100,200,10,15,125,False)

    etc.
    """

    def __init__(self, fpairs=Group(CcdFringePair)):
        """Constructs a :class:`MccdFringePair`.

        Arguments::

          fpairs : (Group)
              Group of :class:`CcdFringePair` objects
        """
        super().__init__(CcdFringePair, fpairs)

    def __repr__(self):
        return "{:s}(defs={:s})".format(self.__class__.__name__, super().__repr__())

    def write(self, fname):
        """Dumps a MccdFringePair in JSON format to a file called fname"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ["hipercam.MccdFringePair"] + \
            [
                (key, ["hipercam.fringe.CcdFringePair"] + list(val.items()))
                for key, val in self.items()
            ]

        with open(fname, "w") as fp:
            json.dump(listify, fp, cls=_Encoder, indent=2)

    def toString(self):
        """Returns MccdFringePair in JSON format as a string"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ["hipercam.fringe.MccdFringePair"] + list(
            (
                (key, ["hipercam.fringe.CcdFringePair"] + list(val.items()))
                for key, val in self.items()
            )
        )
        return json.dumps(listify, cls=_Encoder, indent=2)

    @classmethod
    def read(cls, fname):
        """Read from JSON-format file fname.

          fp : a file-like object opened for reading of text

        Returns an MccdFringePair object.

        """
        with open(fname) as fp:
            obj = json.load(fp, cls=_Decoder)
        listify = [(v1, CcdFringePair(v2[1:])) for v1, v2 in obj[1:]]
        mccd_def = MccdFringePair(listify)
        return mccd_def


# classes to support JSON serialisation of FringePair objects
class _Encoder(json.JSONEncoder):
    def default(self, obj):

        if isinstance(obj, FringePair):
            return OrderedDict(
                (
                    ("Comment", "hipercam.fringe.FringePair"),
                    ("x1", obj.x1),
                    ("y1", obj.y1),
                    ("x2", obj.x2),
                    ("y2", obj.y2),
                )
            )

        return super().default(obj)


class _Decoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        if "Comment" in obj and obj["Comment"] == "hipercam.fringe.FringePair":
            return FringePair(obj["x1"], obj["y1"], obj["x2"], obj["y2"])

        return obj
