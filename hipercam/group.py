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
    type. Subclassed from :class:`OrderedDict`, this class checks the
    consistent types of the objects, that they don't `clash`, and that
    the keys are either integers or string. Dictionaries are
    used to allow flexible indexing; the use of :class:`OrderedDict` is to
    preserve the ordering which allows some short cuts when comparing related
    Group objects. The meaning of `clash` is defined by the objects which must
    support a method `clash` with signature `clash(self, other)` that returns
    True if the object clashes with `other`.

    The restriction to integer and string keys is to allow easy writing to
    FITS files.

    """

    def __init__(self, odct=None):
        """Group constructor from an OrderedDict.

        Arguments::

          odct : (:class:`OrderedDict`)
               An :class:`OrderedDict` containing objects of identical type, keyed on
               integers and/or strings.

        """
        if isinstance(odct, OrderedDict):
            super(Group,self).__init__(odct)
        elif odct is None:
            super(Group,self).__init__()
        else:
            raise HipercamError(
                    'Group.__init__: odct must be None or an OrderedDict')

        # rather un-"pythonic" level of checking here, but better IMO
        # in this case to fail during construction than at some later
        # random time which could be horribly obscure.

        # keys must be either integers or strings
        if any(type(key) != type(0) and type(key) != type('')
               for key in self.keys()):
            raise HipercamError(
                    'Group.__init__: odct keys must be integers or strings')

        if len(self) > 1:

            # iterator over stored objects
            oiter = iter(self.values())

            # type of first object
            tfobj = type(next(oiter))

            # check the rest have the same type
            if any(type(obj) != tfobj for obj in oiter):
                raise HipercamError(
                    'Group.__init__: more than one object type')

            # check for clashes, an N*(N-1)/2 problem
            objs = list(self.values())
            for i, obj1 in enumerate(objs):
                if any(obj.clash(obj1) for obj in objs[i+1:]):
                    raise HipercamError(
                        'Group.__init__: object clash encountered')

    def __setitem__(self, key, item):
        """Adds an item `item` keyed by `key`
        checking that its type matches and that it does
        does clash with any current member of the :class:
        `Group`.
        """

        if type(key) != type(0) and type(key) != type(''):
            raise HipercamError(
                'Group.__setitem__: key must be an integer or a string')

        # check that the new item has same type as current ones
        if any(type(obj) != type(item) for obj in self.values()):
            raise HipercamError(
                'Group.__setitem__: type of item differs from existing Group members')

        # check that item does not clash with any current one
        if any(obj.clash(item) for obj in self.values()):
            raise HipercamError(
                'Group.__setitem__: object clash encountered')

        # checks passed, set the new item
        super(Group,self).__setitem__(key, item)


    def __repr__(self):
        """Returns string representation of object"""
        return 'Group(odct=' + super(Group,self).__repr__() + ')'

