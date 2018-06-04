# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes to represent CCD defects. Base class Defect from
which others are inherited. The derived classes are :class:`Point`,
:class:`Line` and :class:`Hot` representing point, line and hot pixel
defects.
"""

import numpy as np
import json
from collections import OrderedDict
from abc import ABC, abstractmethod
from enum import Enum
from .core import *
from .group import *
from . import utils

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

class Point(Defect):

    """Single Point defect.

    Attributes::

      severity   : Severity
         severity of the Defect

      x, y       : float, float
         coordinates of the pixel, starting (1,1) at lower-left
    """

    def __init__(self, severity, x, y):
        """Constructor. Arguments::

           severity     : Severity
              indicator of the severity of a defect

           x, y         : float, float
              coordinates of the pixel, starting (1,1) at lower-left
        """
        super().__init__(severity)
        self.x = x
        self.y = y

    def copy(self, memo=None):
        """Returns with a copy of the Point"""
        return Point(self.severity, self.x, self.y)

    def __repr__(self):
        return 'Point(severity={!r}, x={!r}, y={!r})'.format(
            self.severity, self.x, self.y
        )

    def dist(self, x, y):
        return np.sqrt((x-self.x)**2+(y-self.y)**2)

class Line(Defect):

    """Line defect class.

    Attributes::

      severity   : Severity
         severity of the Defect

      x1, y1     : float, float
         coordinates of one end of the line, starting (1,1) at lower-left

      x2, y2     : float, float
         coordinates of the other end of the line, starting (1,1) at lower-left
    """

    def __init__(self, severity, x1, y1, x2, y2):
        """Constructor. Arguments::

           severity     : Severity
              indicator of the severity of a defect

           x1, y1         : float, float
              coordinates of one end of the line, starting (1,1) at lower-left

           x2, y2         : float, float
              coordinates of the other end of the line
        """
        super().__init__(severity)
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2

    def copy(self, memo=None):
        """Returns with a copy of the Line"""
        return Line(self.severity, self.x1, self.y1, self.x2, self.y2)

    def __repr__(self):
        return 'Line(severity={!r}, x1={!r}, y1={!r}, x2={!r}, y2={!r})'.format(
            self.severity, self.x1, self.y1, self.x2, self.y2
        )

    def dist(self, x, y):
        """Defines 'distance' from a line defect as the minimum distance
        of a point (x,y) from any portion of the line between its two
        end-points
        """
        p = utils.Vec2D(x,y)
        l1 = utils.Vec2D(self.x1, self.y1)
        l2 = utils.Vec2D(self.x2, self.y2)
        pl1 = p-l1
        pl2 = p-l2

        # first the minimum distance from either end-point
        d = min(pl1.length(), pl2.length())

        l12 = l2-l1
        ll12 = l12.length()
        if ll12 > 0:
            lam = utils.dot(pl1,l12)/ll12**2
            if lam > 0 and lam < 1:
                d = min(d, (p-l1-lam*l12).length())
        return d

class Hot(Defect):

    """Hot pixel defect.

    Attributes::

      severity   : Severity
         severity of the Defect

      x, y       : float, float
         coordinates of the pixel, starting (1,1) at lower-left
    """

    def __init__(self, severity, x, y):
        """Constructor. Arguments::

           severity     : Severity
              indicator of the severity of a defect

           x, y         : float, float
              coordinates of the pixel, starting (1,1) at lower-left
        """
        super().__init__(severity)
        self.x = x
        self.y = y

    def copy(self, memo=None):
        """Returns with a copy of the Hot"""
        return Hot(self.severity, self.x, self.y)

    def __repr__(self):
        return 'Hot(severity={!r}, x={!r}, y={!r})'.format(
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

        if isinstance(obj, Point):
            return OrderedDict(
                (
                    ('Comment', 'hipercam.defect.Point'),
                    ('severity', obj.severity.name),
                    ('x', obj.x),
                    ('y', obj.y),
                    )
                )

        elif isinstance(obj, Line):
            return OrderedDict(
                (
                    ('Comment', 'hipercam.defect.Line'),
                    ('severity', obj.severity.name),
                    ('x1', obj.x1),
                    ('y1', obj.y1),
                    ('x2', obj.x2),
                    ('y2', obj.y2),
                    )
                )

        elif isinstance(obj, Hot):
            return OrderedDict(
                (
                    ('Comment', 'hipercam.defect.Hot'),
                    ('severity', obj.severity.name),
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
        if 'Comment' in obj:
            if obj['Comment'] == 'hipercam.defect.Point':
                return Point(
                    getattr(Severity,obj['severity']), obj['x'], obj['y']
                )
            elif obj['Comment'] == 'hipercam.defect.Line':
                return Line(
                    getattr(Severity,obj['severity']), obj['x1'], obj['y1'],
                    obj['x2'], obj['y2']
                )
            elif obj['Comment'] == 'hipercam.defect.Hot':
                return Hot(
                    getattr(Severity,obj['severity']), obj['x'], obj['y']
                )

        return obj
