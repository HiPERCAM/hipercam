# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes which contain groups of objects of identical type. Each
object must support a method called "clash" where clash(self, other)
returns True if self and other conflict is some way. For instance,
two CCD sub-windows might be said to clash if they contain any pixels
in common.
"""

# Imports for 2 / 3 compatibility
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

from collections import OrderedDict

from .core import *


class Group(OrderedDict):
    """Base class for a container of multiple objects of identical
    type. Subclassed from OrderedDict, this class adds checks on consistency
    of the objects and that they don't "clash". Dictionaries are used to allow
    flexible indexing; the use of OrderedDict is to preserve the ordering
    which allows some short cuts when comparing related Group objects. The
    meaning of "clash" is defined by the objects which must support a method
    "clash".

    """

    def __init__(self, dct):
        super(Group,self).__init__(dct)

        # check that all objects are of the same type
        vals = iter(self.values())
        first_type = type(next(vals))
        for obj in vals:
            if not isinstance(obj, first_type):
                raise HipercamError(
                    'Group.__init__: more than one object type')

        # check for any conflicts
        vals = list(self.values())
        for i,val1 in enumerate(vals):
            for val2 in vals[i+1:]:
                if val1.clash(val2):
                    raise HipercamError(
                        'Group.__init__: object clash encountered')

    def __repr__(self):
        return 'Group(dct=' + super(Group,self).__repr__() + ')'

    def __setitem__(self, key, item):
        """Adds an item `item` keyed by `key`
        checking that its type matches and that it does
        does clash with any current member of the :class:
        `Group`.
        """
        if len(self):
            # check that new item has same type as current ones
            type_item = type(item)
            for k,obj in self.items():
                if not isinstance(obj,type_item):
                    raise HipercamError(
                        'Group.__setitem__: object type differs from existing Group members')
                if item.clash(obj):
                    raise HipercamError(
                        'Group.__setitem__: object clash encountered')

        # OK now set the new item
        super(Group,self).__setitem__(key, item)
