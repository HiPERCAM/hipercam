# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to represent CCD defects. Base class Defect from
which others are inherited.
"""

import numpy as np
import json
from collections import OrderedDict
from abc import ABC, abstractmethod
from enum import Enum
from .core import *
from .group import *

__all__ = ('Severity', 'Defect',)

class Severity(Enum):
    """
    Enum class to define defect severity levels. There are
    two: MODERATE and SEVERE
    """
    MODERATE = 1
    SEVERE = 2

class Defect(ABC):

    """Abstract class representing a CCD defect.

    Attributes::

      severity   : Severity
         severity of the Defect
    """

    def __init__(self, severity):
        """Constructor. Arguments::

           severity     : Severity
              indicator of the severity of a defect
        """

        self.severity = severity

    @abstractmethod
    def copy(self, memo=None):
        """Returns with a copy of the Aperture"""
        return Aperture(
            self.x, self.y, self.rtarg, self.rsky1, self.rsky2,
            self.ref, self.mask.copy(), self.extra.copy(),
            self.link
        )

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def dist(self, x, y):
        """Returns the distance of a defect from a position x,y. This defines
        what is meant by the distance from a defect
        """
        pass

    def write(self, fname):
        """Dumps Aperture in JSON format to a file called fname"""
        with open(fname,'w') as fp:
            json.dump(self, cls=_Encoder, indent=2)

    def toString(self):
        """Returns Aperture as a JSON-type string"""
        return json.dumps(self, fp, cls=_Encoder, indent=2)

    @classmethod
    def read(cls, fname):
        """Read from JSON-format file fname"""
        with open(fname) as fp:
            aper = json.load(fp, cls=_Decoder)
        aper.check()
        return aper

class Pixel(ABC):

    """Single pixel defect.

    Attributes::

      severity   : Severity
         severity of the Defect

      x, y       : int, int
         coordinates of the pixel, starting (1,1) at lower-left
    """

    def __init__(self, severity, x, y):
        """Constructor. Arguments::

           severity     : Severity
              indicator of the severity of a defect

           x, y         : int, int
              coordinates of the pixel, starting (1,1) at lower-left
        """
        super().__init__(severity)
        self.x = x
        self.y = y

    def copy(self, memo=None):
        """Returns with a copy of the Pixel"""
        return Pixel(self.severity, self.x, self.y)

    def __repr__(self):
        return 'Pixel(severity={!r}, x={!r}, y={!r})'.format(
            self.severity, self.x, self.y
        )

    def dist(self, x, y):
        return np.sqrt((x-self.x)**2+(y-self.y)**2)

class CcdDefect(Group):
    """Class representing all the :class:Defects for a single CCD.
    Normal usage is to create an empty one and then add Defects via
    the usual mechanism for updating dictionaries, i.e. ccddef[label] =
    defect
    """

    def __init__(self, defs=Group(Defect)):
        """Constructs a :class:`CcdDefect`.

        Arguments::

          defs : Group
              Group of :class:`Defect` objects
        """
        super().__init__(Defect, defs)

    def __repr__(self):
        return '{:s}(defs={:s})'.format(
            self.__class__.__name__, super().__repr__()
            )

    def write(self, fname):
        """Dumps ccdAper in JSON format to a file called fname"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ['hipercam.CcdDefect'] + list(self.items)
        with open(fname,'w') as fp:
            json.dump(listify, fp, cls=_Encoder, indent=2)

    def copy(self, memo=None):
        return CcdDefect(
            super().copy(memo)
        )

class MccdDefect(Group):
    """Class representing all the :class:Defects for multiple CCDs.
    Normal usage is to create an empty one and then add apertures via
    the usual mechanism for updating dictionaries, e.g.

      >> mccddef = MccdDefect()
      >> mccddef['ccd1'] = CcdDefect()
      >> mccddef['ccd2'] = CcdDefect()
      >> mccddef['ccd1']['ap1'] = Defect(100,200,10,15,125,False)

    etc.
    """

    def __init__(self, defs=Group(CcdDefect)):
        """Constructs a :class:`MccdDefect`.

        Arguments::

          defs : (Group)
              Group of :class:`CcdDefect` objects
        """
        super().__init__(CcdDefect, defs)

    def __repr__(self):
        return '{:s}(defs={:s})'.format(
            self.__class__.__name__, super().__repr__()
            )

    def write(self, fname):
        """Dumps a MccdDefect in JSON format to a file called fname"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ['hipercam.MccdDefect'] + list(
            ((key,['hipercam.CcdDefect']+list(val.items())) \
             for key, val in self.items())
        )
        with open(fname,'w') as fp:
            json.dump(listify, fp, cls=_Encoder, indent=2)

    def toString(self):
        """Returns MccdDefect in JSON format as a string"""

        # dumps as list to retain order through default iterator encoding
        # that buggers things otherwise
        listify = ['hipercam.MccdDefect'] + list(
            ((key,['hipercam.CcdDefect']+list(val.items())) \
             for key, val in self.items())
        )
        return json.dumps(listify, cls=_Encoder, indent=2)

    @classmethod
    def read(cls, fname):
        """Read from JSON-format file fname.

          fp : a file-like object opened for reading of text

        Returns an MccdDefect object.

        """
        with open(fname) as fp:
            obj = json.load(fp, cls=_Decoder)
        listify = [(v1,CcdDefect(v2[1:])) for v1,v2 in obj[1:]]
        mccd_def = MccdDefect(listify)
        return mccd_def

# classes to support JSON serialisation of Defect objects
class _Encoder(json.JSONEncoder):

    def default(self, obj):

        if isinstance(obj, Pixel):
            return OrderedDict(
                (
                    ('Comment', 'hipercam.defect.Pixel'),
                    ('severity', obj.severity),
                    ('x', obj.x),
                    ('y', obj.y),
                    )
                )

        return super().default(obj)

class _Decoder(json.JSONDecoder):

    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        # look out for Defect objects. Everything else done by default
        if 'Comment' in obj and obj['Comment'] == 'hipercam.defect.Pixel':
            return Pixel(
                obj['severity'], obj['x'], obj['y']
            )

        return obj
