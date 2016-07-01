# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines classes which contain groups of identical objects. Each
object must support a method called "clash" where clash(self, other)
returns True if self and other conflict is some way. For instance,
two CCD sub-windows might be said to clash if they contain any pixels
in common.
"""

# Standard pre-amble from astropy
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.extern import six

import numpy as np
from .core import *

class Group(dict):
    """Base class for a container of multiple
    objects of identical type. Subclassed from dictionaries,
    this class adds checks on consistency of the objects
    and that they don't clash. Dictionaries are used to allow
    flexible indexing.
    """

    def __init__(self, dct):
        super(Group,self).__init__(dct)
        if len(self) > 1:
            # check that all objects are of the same type
            tone = type(self.values()[0])
            for obj in self.values()[1:]:
                if not isinstance(obj,tone):
                    raise HipercamError(
                        'Group.__init__: more than one object type')

            # check for any conflicts
            keys = self.keys()
            for i in xrange(len(keys)-1):
                for j in xrange(i+1,len(keys)):
                    if self[keys[i]].clash(self[keys[j]]):
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
            tone = type(item)
            for ob in self:
                if not isinstance(ob,tone):
                    raise HipercamError(
                        'Group.__setitem__: object type differs from existing Group members')
                if obj.clash(ob):
                    raise HipercamError(
                        'Group.__setitem__: object clash encountered')

        # OK now set the new item
        super(Group,self).__setitem__(key, item)

