# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines classes which contain groups of objects of identical type. Can
check for conflicts between the objects if they support a method called
`clash` with signature `clash(self, other)` which raises an exception if
`self` and `other` conflict in some way. For instance, two CCD sub-windows
might be said to clash if they contain any pixels in common.

"""

import copy
from collections import OrderedDict

from .core import *

__all__ = ('Group', 'Agroup')

class Group(OrderedDict):
    """A specialized OrderedDict for storing objects which all match in terms of
    "instance" and are indexed by strings only. It will check for conflicts
    between the stored objects if the objects have a method `clash` with
    signature `clash(self, other)` which raises an exception if `self` and
    `other` conflict in some way. The objects should support a `copy` method
    to return a deepcopy. This is used in the :class:`Group`s copy operation.

    """

    def __init__(self, *args, **kwargs):
        """Group constructor; see `dict` for possible arguments. The keys are checked
        to ensure that they are integers, and the object types are checked to
        see that they are identical. The order of the keys is tracked to help
        with file output.

        Sets attribute `otype`, the type of the stored objects (used for
        simple input checks)

        """
        super().__init__(*args, **kwargs)

        # rather un-"pythonic" level of checking here, but better IMO
        # in this case to fail during construction than at some later
        # random time which could be horribly obscure.

        # keys must be strings
        if any(not isinstance(key,str) for key in self.keys()):
            raise HipercamError('Group.__init__: keys must be strings')

        if len(self):
            # iterator over stored objects
            oiter = iter(self.values())

            # type of first object
            self.otype = type(next(oiter))

            # check the rest have the same type
            if any(not isinstance(obj, self.otype) for obj in oiter):
                raise HipercamError('Group.__init__: more than one object type')

            try:
                # check for clashes, an N*(N-1)/2 problem. 'clash'
                # should raise an exception if there is a problem.
                objs = list(self.values())
                for i, ob in enumerate(objs):
                    for obj in objs[i+1:]:
                        ob.clash(obj)
            except AttributeError:
                # we are OK if no 'clash' is defined
                pass

    def __setitem__(self, key, item):
        """Adds an item `item` keyed by `key`
        checking that its type matches and that it does
        does clash with any current member of the :class:
        `Group`.
        """
        if not instance(key,str):
            raise KeyError(
                'Group.__setitem__: key must be an string')

        # store or check that the new item matches in type
        if not hasattr(self, 'otype'):
            if len(self):
                raise ValueError(
                    'Group.__setitem__: non-zero elements but otype not defined')
            else:
                self.otype = type(item)
        else:
            if not isinstance(item, self.otype):
                raise HipercamError(
                    'Group.__setitem__: key = ' + str(key) + ', item type (=' + str(type(item)) +
                    ') differs from existing Group data type (=' + str(self.otype) + ')')

        try:
            # check that the new item does not clash with any current one
            # clash should raise an exception if there is a problem
            for obj in self.values():
                item.clash(obj)

        except AttributeError:
            # ok if no 'clash' defined
            pass

        # checks passed, set the new item and add the key
        super().__setitem__(key, item)

    def copy(self, memo=None):
        """Copy operation. The stored objects must have a `copy(self, memo)` method.

        """
        group = Group()
        group.otype = self.otype
        for key, val in self.items():
            group[key] = val.copy(memo)
        return group

    def __copy__(self):
        """Copy operation for copy.copy

        """
        return self.copy()

    def __deepcopy__(self, memo):
        """Copy operation for copy.deepcopy

        """
        return self.copy(memo)

    def __repr__(self):
        return 'Group([{}])'.format(
            ', '.join('({!r}, {!r})'.format(key,val)
                      for key, val in self.items())
        )

class Agroup(Group):
    """A :class:`Group` which defines arithmetic methods +=, +, etc which must be
    supported by whatever objects the :class:`Group` contains. This allows the
    same approach to be taken when operating on :class:`CCD` and :class:`MCCD`
    objects which are based on this class. The objects stored in the
    :class:`Agroup` must define a method `copy` which returns a deep copy;
    this is required by several operators.

    """

    def copy(self, memo=None):
        """Copy operation.

        """
        return Agroup(super().copy(memo))

    def __iadd__(self, other):
        """Adds `other` to the :class:`Agroup` as 'self += other'. If `other` is
        another :class:`Agroup` with the same object type (`otype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to add to each object in the :class:`Agroup`.

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
        """Subtracts `other` from the :class:`Agroup` as 'self -= other'. If `other`
        is another :class:`Agroup` with the same object type (`otype`) as
        self, then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to subtract from each object in the :class:`Agroup`.

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
        """Multiplies the :class:`Agroup` by `other` as 'self *= other'. If `other` is
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
                    obj *= other[key]
        else:
            for obj in self.values():
                obj *= other

        return self

    def __itruediv__(self, other):
        """Divides the :class:`Agroup` by `other` as 'self /= other'. If `other` is
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
                    obj /= other[key]
        else:
            for obj in self.values():
                obj /= other

        return self

    def __add__(self, other):
        """Adds `other` to the :class:`Agroup` as '= self + other'. If `other` is
        another :class:`Agroup` with the same object type (`otype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to add to each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """
        cself = self.copy()
        cself += other
        return cself

    def __sub__(self, other):
        """Subtracts `other` from the :class:`Agroup` as '= self - other'. If `other` is
        another :class:`Agroup` with the same object type (`otype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to subtract from each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """
        cself = self.copy()
        cself -= other
        return cself

    def __mul__(self, other):
        """Multiplies the :class:`Agroup` by `other` as '= self * other'. If `other`
        is another :class:`Agroup` with the same object type (`otype`) as
        self, then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to multiply each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """
        cself = self.copy()
        cself *= other
        return cself

    def __truediv__(self, other):
        """Divides the :class:`Agroup` by `other` as '= self / other'. If `other` is
        another :class:`Agroup` with the same object type (`otype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to divide into each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """
        cself = self.copy()
        cself /= other
        return cself

    def __radd__(self, other):
        """Adds the :class:`Agroup` to other as '= other + self'.
        """
        return self.__add__(other)

    def __rsub__(self, other):
        """Subtracts the :class:`Agroup` from `other` as '= other - self'.
        """
        cself = self.copy()
        cself *= -1. # in-place faster
        cself += other
        return cself

    def __rmul__(self, other):
        """Multiplies the :class:`Agroup` by `other` as '= other * self'.
        """
        return self.__mul__(other)

    def __truediv__(self, other):
        """Divides a :class:`Agroup` by `other` as '= other / self'.
        """
        cself = self.copy()
        cself /= other
        return cself

    def __repr__(self):
        return 'Agroup([{}])'.format(
            ', '.join('({!r}, {!r})'.format(key,val)
                      for key, val in self.items())
        )
