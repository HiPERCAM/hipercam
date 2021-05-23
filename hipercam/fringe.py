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
        # record which window encloses this
        # for speed in operations
        self.wnam = None

    def copy(self, memo=None):
        """Returns with a copy of the Aperture"""
        return FringePair(
            self.x1,
            self.y1,
            self.x2,
            self.y2
        )

    def dist(self, x, y):
        """
        Defines a distance of the FringePair from a point
        """
        ex = self.x2-self.x1
        ey = self.y2-self.y1
        ax = x-self.x1
        ay = y-self.y1
        if ex == 0. and ey == 0.:
            d = np.sqrt(ax**2+ay**2)
        else:
            ell = np.sqrt(ex**2+ey**2)
            ex /= ell
            ey /= ell
            lam = ex*ax+ey*ay
            lam = max(0,min(lam,ell))
            xm = self.x1 + lam*ex
            ym = self.y1 + lam*ey
            d = np.sqrt((x-xm)**2+(y-ym)**2)
        return d

    def inside(self, ccd, nhalf):
        # the nhalf doesn't really test properly in the next two lines
        # because it is measured in unbinned pixels
        wnam1 = ccd.inside(self.x1, self.y1, nhalf)
        wnam2 = ccd.inside(self.x2, self.y2, nhalf)
        if wnam1 is not None and wnam2 is not None and wnam1 == wnam2:
            self.wnam = wnam1

            # Final check that we are far enough away from the edge to
            # measure full averages
            wind = ccd[self.wnam]
            xlo, xhi, ylo, yhi = wind.extent()

            if self.x1 <= xlo + nhalf*wind.xbin or \
               self.x2 <= xlo + nhalf*wind.xbin or \
               self.x1 >= xhi - nhalf*wind.xbin or \
               self.x2 >= xhi - nhalf*wind.xbin or \
               self.y1 <= ylo + nhalf*wind.ybin or \
               self.y2 <= ylo + nhalf*wind.ybin or \
               self.y1 >= yhi - nhalf*wind.ybin or \
               self.y2 >= yhi - nhalf*wind.ybin:
                self.wnam = None
        else:
            self.wnam = None
        return self.wnam

    def diff(self, ccd, nhalf):
        """Returns difference in intensity. Only run after having run
        "inside" with the same value of nhalf because it assumes the
        window name has been set

        nhalf : int
           will measure sum over region +/-nhalf binned pixels around
           nearest pixel of either end of FringePair
        """

        wind = ccd[self.wnam]

        # pixel centres
        ix1 = int(round(wind.x_pixel(self.x1)))
        iy1 = int(round(wind.y_pixel(self.y1)))
        ix2 = int(round(wind.x_pixel(self.x2)))
        iy2 = int(round(wind.y_pixel(self.y2)))

        # regions
        xlo1, xhi1 = ix1-nhalf, ix1+nhalf+1
        ylo1, yhi1 = iy1-nhalf, iy1+nhalf+1
        xlo2, xhi2 = ix2-nhalf, ix2+nhalf+1
        ylo2, yhi2 = iy2-nhalf, iy2+nhalf+1

        # extract sums over the regions
        mean1 = wind.data[ ylo1:yhi1, xlo1:xhi1].mean()
        mean2 = wind.data[ ylo2:yhi2, xlo2:xhi2].mean()

        return mean2 - mean1

    def __repr__(self):
        return \
            f"FringePair(x1={self.x1:.2f}, y1={self.y1:.2f}," \
            f" x2={self.x2:.2f}, y2={self.y2:.2f})"

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
        super().__init__(FringePair, fpairs)

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

    def crop(self, ccd, nhalf):
        """
        Returns version of CcdFringePair with any
        *not* in ccd removed.
        """
        tccdfpair = CcdFringePair()
        for fpnam, fpair in self.items():
            if fpair.inside(ccd, nhalf) is not None:
                tccdfpair[fpnam] = fpair
        return tccdfpair

    def diff(self, ccd, nhalf):
        """Returns array of differences for all FringePairs in the ccd "crop"
        should have been run using ccd before running this.

        """
        diffs = np.array([fpair.diff(ccd, nhalf) for fpair in self.values()])
        return diffs

    def scale(self, ccd, ccdref, nhalf):
        """
        Measures scale factor needed to subtract the fringes in
        a reference CCD `ccdref` from `ccd`. Measures median of
        the FringePair amplitude ratios.
        """
        dccd = self.diff(ccd, nhalf)
        dccdref = self.diff(ccdref, nhalf)
        ratios = dccd / dccdref
        return np.nanmedian(ratios)

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

    def crop(self, mccd, nhalf):
        tmccdfpair = MccdFringePair()
        for cnam in mccd:
            if cnam in self:
                tccdfpair = self[cnam].crop(mccd[cnam], nhalf)
                if len(tccdfpair):
                    tmccdfpair[cnam] = tccdfpair
        return tmccdfpair

    def __repr__(self):
        return "{:s}(defs={:s})".format(self.__class__.__name__, super().__repr__())

    def write(self, fname):
        """Dumps a MccdFringePair in JSON format to a file called fname"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ["hipercam.fringe.MccdFringePair"] + \
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
