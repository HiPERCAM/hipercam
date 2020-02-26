# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines classes which contain groups of objects of identical type
keyed with strings only.

"""

import copy
from collections import OrderedDict

from .core import *

__all__ = ('Group', 'Agroup')

class Group(OrderedDict):
    """A specialized OrderedDict that enforces string keys only and a
    single class for the values. The objects stored as values should
    support a `copy` method to return a deepcopy for use in the
   :class:`Group`s copy operation.

    :class:`Group` is designed to amalgamate several of the classes used for
    objects associated with CCDs such as :class:`Window` (see :class:`CCD`),
    and :class:`Aperture` (see :class:`CcdAper`)

    """

    def __init__(self, ftype, *args, **kwargs):
        """Group constructor; see :class:`OrderedDict` for possible values of
        args and kwargs. The first argument, `ftype`, is the fixed
        type that will be enforced. e.g. it could be Window, Winhead,
        CCD etc. An attribute of the same name will be stored.

        Arguments::

          ftype : type
            the type of the objects that will be allowed.

        """
        # set the example type
        self.ftype = ftype

        # set the key / value pairs
        super().__init__(*args, **kwargs)

        # rather un-"pythonic" level of checking here, but better IMO
        # in this case to fail during construction than at some later
        # random time which could be horribly obscure.

        # keys must be strings
        if any(not isinstance(key,str) for key in self.keys()):
            raise HipercamError('keys must be strings')

        # values must have specific type
        if any(not isinstance(obj,ftype) for obj in self.values()):
            raise HipercamError(
                'an input object is not a {!s} object'.format(ftype)
            )

    def __setitem__(self, key, item):
        """Adds an item `item` keyed by `key`
        checking that its type matches 
        """
        if not isinstance(key,str):
            raise KeyError('key must be an string')

        # store or check that the new item matches in type
        if not isinstance(item, self.ftype):
            raise HipercamError(
                'key = {:s}: item = {!s} does not have'
                ' the expected type = {!s}'.format(key,item,self.ftype)
            )

        # checks passed, set the new item and add the key
        super().__setitem__(key, item)

    def copy(self, memo=None):
        """Copy operation. The stored objects must have a `copy(self, memo)` method.

        """
        group = Group(self.ftype)
        for key, val in self.items():
            group[key] = val.copy(memo)
        return group

    def check(self):
        """Checks validity of a Group if its objects have a method
        'clash' defined to check that they are compatible. This should
        have signature 'clash(self, other)', but does not have to be
        defined for all objects storable in Groups. See :class:`Winhead` for
        an example of an object for which this makes sense.

        If the Group fails this check, a ValueError is raised.
        """
        keys = list(self.keys())
        for key1 in keys[:-1]:
            for key2 in keys[1:]:
                try:
                    self[key1].clash(self[key2])
                except AttributeError:
                    # if object doesn't support clash, break out early
                    break

    def get_num(self, num):
        """Returns the num element of the Group (starting from 0)"""
        keys = list(self.keys())
        return self[keys[num]]

    def __copy__(self):
        """Copy operation for copy.copy

        """
        return self.copy()

    def __deepcopy__(self, memo):
        """Copy operation for copy.deepcopy

        """
        return self.copy(memo)

    def __repr__(self):
        return 'Group(ftype={!r}, [{}])'.format(
            self.ftype,
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
        agroup = Agroup(self.ftype)
        for key, val in self.items():
            agroup[key] = val.copy(memo)
        return agroup

    def __iadd__(self, other):
        """Adds `other` to the :class:`Agroup` as 'self += other'. If `other` is
        another :class:`Agroup` with the same object type (`ftype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to add to each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.ftype == self.ftype:
            for key, obj in self.items():
                if key in other:
                    obj += other[key]
        else:
            for obj in self.values():
                obj += other

        return self

    def __isub__(self, other):
        """Subtracts `other` from the :class:`Agroup` as 'self -= other'. If `other`
        is another :class:`Agroup` with the same object type (`ftype`) as
        self, then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to subtract from each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.ftype == self.ftype:
            for key, obj in self.items():
                if key in other:
                    obj -= other[key]
        else:
            for obj in self.values():
                obj -= other

        return self

    def __imul__(self, other):
        """Multiplies the :class:`Agroup` by `other` as 'self *= other'. If `other` is
        another :class:`Agroup` with the same object type (`ftype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to multiply each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.ftype == self.ftype:
            for key, obj in self.items():
                if key in other:
                    obj *= other[key]
        else:
            for obj in self.values():
                obj *= other

        return self

    def __itruediv__(self, other):
        """Divides the :class:`Agroup` by `other` as 'self /= other'. If `other` is
        another :class:`Agroup` with the same object type (`ftype`) as self,
        then the operation will be applied to each pair of objects with
        matching keys. Otherwise `other` will be regarded as a constant object
        to divide into each object in the :class:`Agroup`.

        In the first case, if self has keys that are not in `other`, then they
        will be left untouched. Any keys in `other` not in self are ignored.

        """

        if isinstance(other, Agroup) and other.ftype == self.ftype:
            for key, obj in self.items():
                if key in other:
                    obj /= other[key]
        else:
            for obj in self.values():
                obj /= other

        return self

    def __add__(self, other):
        """Adds `other` to the :class:`Agroup` as '= self + other'. If `other` is
        another :class:`Agroup` with the same object type (`ftype`) as self,
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
        another :class:`Agroup` with the same object type (`ftype`) as self,
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
        is another :class:`Agroup` with the same object type (`ftype`) as
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
        another :class:`Agroup` with the same object type (`ftype`) as self,
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

    def __rtruediv__(self, other):
        """Divides a :class:`Agroup` by `other` as '= other / self'.
        """

        cself = self.copy()
        for obj in cself.values():
            obj = other / cself

        return self

    def __repr__(self):
        return 'Agroup(ftype={!r}, [{}])'.format(
            self.ftype,
            ', '.join('({!r}, {!r})'.format(key,val)
                      for key, val in self.items())
        )
