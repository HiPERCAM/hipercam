# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines a class to contain header information compatible with FITS headers.
The aim is to provide a faster replacement for astropy.io.fits.Header
"""

from collections import OrderedDict as odict
from astropy.io import fits

__all__ = ('Header',)

class Header(odict):

    """Simulates some basic functionality of astropy.io.fits.Header
    objects while trying to be more efficient. Values are stored as
    two element tuples consisting of the actual value and a comment.
    Extra comments and history are allowed for as separate entities.
    The idea is to maintain headers that are compatible with writing
    to and reading from FITS, so for instance, they are case
    insensitive and 'HIERARCH ' at the start will be stripped.

    """

    def __init__(self, od=odict()):
        self.comments = []
        self.history = []
        super().__init__(od)


    @classmethod
    def from_fits(cls, fhead):
        """Sets a Header from a FITS header

        Arguments::

          fhead : astropy.io.fits.Header
            the FITS header to copy.

        """
        hd = cls()
        for key,item in fhead.items():
            hd[key] = item
        return hd

    @property
    def to_fits(self):
        """Returns the Header as an astropy.io.fits.Header
        which is useful when saving to disk
        """
        fh = fits.Header()
        for key,item in super().items():
            fh[key] = item

        for comment in self.comments:
            fh.add_comment(comment)

        for history in self.history:
            fh.add_history(history)

        return fh

    @staticmethod
    def _process_key(key):
        key = key.upper()
        if key.startswith('HIERARCH '):
            key = key[9:]
        return key

    def __setitem__(self, key, item):

        key = Header._process_key(key)

        if key == 'COMMENT' or key == 'HISTORY':
            if not isinstance(item,str):
                raise ValueError(
                    'Can only use string values for key = COMMENT or HISTORY'
                )
            if key == 'COMMENT':
                self.comments.append(item)
            elif key == 'HISTORY':
                self.history.append(item)

        elif isinstance(item,tuple):
            if len(item) != 2:
                raise ValueError(
                    'Can only set Header items with '
                    'values or 2-element (value,comment) pairs'
                )
            super().__setitem__(key, item)

        else:

            super().__setitem__(key,(item,''))

    def __getitem__(self, key):
        key = Header._process_key(key)
        return super().__getitem__(key)[0]

    def get(self, key, default):
        key = Header._process_key(key)
        return super().get(key,(default,''))[0]

    def get_full(self, key):
        """Returns with 2-element (value,comment) tuple for the given keyword"""
        key = Header._process_key(key)
        return super().__getitem__(key)

    def get_comment(self, key):
        key = Header._process_key(key)
        return super().__getitem__(key)[1]

    def copy(self):
        hd = Header(super().copy())
        hd.comments = self.comments.copy()
        hd.history = self.history.copy()
        return hd

if __name__ == '__main__':

    fh = fits.Header()
    fh['a'] = 1.23
    fh['A'] = 1.234
    fh['b'] = (2.34,'test comment')
    fh['this is a long keyword'] = (4.5,'test')
    fh['COMMENT'] = ' a comment'
    fh.add_comment('hello 1')
    fh.add_comment('hello 2')
    fh.add_history('did something')
    fh.add_history('did something else')
    print(fh)
    hd = Header.from_fits(fh)
    print(hd)

    fhc = hd.to_fits
    print(fhc)

    print(fhc['a'])
