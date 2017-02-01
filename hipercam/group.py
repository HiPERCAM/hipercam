# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines classes which contain groups of objects of identical type. Each
object must support a method called "clash" where clash(self, other) raises an
exception if `self` and `other` conflict in some way. For instance, two CCD
sub-windows might be said to clash if they contain any pixels in common.

"""

from collections import OrderedDict

from .core import *

class Group(OrderedDict):
    """Base class for a container of objects of identical type. Subclassed from
    :class:`OrderedDict`, this class checks that the objects are of the same
    type, that they don't `clash`, and that the keys are
    integers. Dictionaries are used to allow flexible indexing; the use of
    :class:`OrderedDict` is to preserve the ordering which allows some short
    cuts when comparing related Group objects. The meaning of `clash` is
    defined by the objects which must support a method `clash` with signature
    `clash(self, other)` that raises an exception if the object clashes with
    `other`.

    The restriction to integer keys is to allow easy writing to FITS files
    and referencing within scripts.

    """

    def __init__(self, *args, **kwargs):
        """Group constructor from an :class:`OrderedDict`; see that for
        possible arguments. The keys are checked to ensure that they
        are integers, and the object types are checked to see that they
        are identical.

        Sets attribute `otype`, the type of the stored objects (used for
        simple input checks)

        """
        self.otype = None
        super().__init__(*args, **kwargs)

        # rather un-"pythonic" level of checking here, but better IMO
        # in this case to fail during construction than at some later
        # random time which could be horribly obscure.

        # keys must be integers
        if any(type(key) != type(0) for key in self.keys()):
            raise HipercamError('Group.__init__: keys must be integers')

        if len(self):
            # iterator over stored objects
            oiter = iter(self.values())

            # type of first object
            self.otype = type(next(oiter))

            # check the rest have the same type
            if any(type(obj) != self.otype for obj in oiter):
                raise HipercamError('Group.__init__: more than one object type')

            # check for clashes, an N*(N-1)/2 problem. 'clash'
            # should raise an exception if there is a problem.
            objs = list(self.values())
            for i, ob in enumerate(objs):
                for obj in objs[i+1:]:
                    ob.clash(obj)

    def __setitem__(self, key, item):
        """Adds an item `item` keyed by `key`
        checking that its type matches and that it does
        does clash with any current member of the :class:
        `Group`.
        """
        if type(key) != type(0):
            raise HipercamError(
                'Group.__setitem__: key must be an integer')

        # check that the new item has same type as current ones
        if self.otype is not None and type(item) != self.otype:
            raise HipercamError(
                'Group.__setitem__: key = ' + str(key) + ', item type (=' + str(type(item)) +
                ') differs from existing Group data type (=' + str(self.otype) + ')')

        # check that the new item does not clash with any current one
        # clash should raise an exception if there is a problem
        for obj in self.values():
            item.clash(obj)

        # checks passed, set the new item
        super(Group,self).__setitem__(key, item)

class Agroup(Group):
    """A :class:`Group` which defines arithmetic methods +=, +, etc which must be
    supported by whatever objects the :class:`Group` contains. This allows the
    same approach to be taken when operating on :class:`CCD` and :class:`MCCD`
    objects which are based on this class.

    """

    def __iadd__(self, other):
        """Adds `other` to a :class:`Agroup` as 'self += other'. If `other` is another
        :class:`Agroup` with the same object type (`otype`) as self, then the
        operation will be applied to each pair of objects with matching
        keys. Otherwise `other` will be regarded as a constant object to
        add to each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.otype == self.otype:
            for key, obj in self.items():
                if key in other:
                    obj += other[key]
        else:
            for obj in self.values():
                obj += other

        return self

    def __isub__(self, other):
        """Subtracts `other` from a :class:`Agroup` as 'self -= other'. If `other` is another
        :class:`Agroup` with the same object type (`otype`) as self, then the
        operation will be applied to each pair of objects with matching
        keys. Otherwise `other` will be regarded as a constant object to
        subtract from each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.otype == self.otype:
            for key, obj in self.items():
                if key in other:
                    obj -= other[key]
        else:
            for obj in self.values():
                obj -= other

        return self

    def __imul__(self, other):
        """Multiplies a :class:`Agroup` by `other` as 'self *= other'. If `other` is
        another :class:`Agroup` with the same object type (`otype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to multiply each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.otype == self.otype:
            for key, obj in self.items():
                if key in other:
                    obj *= other[nwin]
        else:
            for obj in self.values():
                obj *= other

        return self

    def __itruediv__(self, other):
        """Divides a :class:`Agroup` by `other` as 'self /= other'. If `other` is
        another :class:`Agroup` with the same object type (`otype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to divide into each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.otype == self.otype:
            for key, obj in self.items():
                if key in other:
                    obj /= other[nwin]
        else:
            for obj in self.values():
                obj /= other

        return self


