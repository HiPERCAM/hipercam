# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines a class to contain header information compatible with FITS headers.
The aim is to provide a faster replacement for astropy.io.fits.Header
"""

from collections import OrderedDict as odict
from astropy.io.fits import Header as FITS_Header

__all__ = ('Header',)

class Header:

    """Simulates the basic functionality of astropy.io.fits.Header objects while
    being lighter weight and thus more efficient. The idea is to allow headers
    that are compatible with writing to and reading from FITS files, i.e.
    that have (uppercase) keywords, values, comments and allows for 'HIERARCH'
    long keywords. History and comments are added at the end of the
    header. There is no attempt to replicate all the methods of
    astropy.io.fits.Header, although many are similar in nature to reduce the
    number of changes to the code required to implement this. This class has
    does not have all the checks of astropy.io.fits.Header objects so it is
    possible to create invalid headers, but this should become obvious when
    attempting conversion to FITS. Some example code::

      >> import Header from hipercam
      >> head = Header()
      >> head['a'] = (1.23, 'a comment')
      >> head['object'] = 'IP Peg'

    :class:`Header` objects are used as the base class of :class:`Winhead` and
    :class:`Window` objects which build into :class:`CCD` and :class:`MCCD`
    objects.
    """

    SPECIAL_KEYWORDS = ('COMMENT', 'HISTORY', '')

    def __init__(self, head=None, copy=False):
        """Initialiser. 'head' can be (i) another Header, (ii) an ordered
        dictionary with values set to 2-element tuples containing
        (value,comment) pairs, (iii) a list of (key,value,comment) tuples
        similar to astropy.io.fits.Header 'cards' or (iv) an
        astropy.io.fits.Header object.

        No comments, history or blank lines can be passed via the ordered
        dictionary option, i.e. it is for genuine header only.

        Arguments::

          head : list of (k,v,c) tuples | Header | OrderedDict : fits.Header
            initial data to set the header.

          copy : bool
            if 'head' is a Header or a list of (key,value,comment) tuples,
            this controls whether it the data are copied by value (copy=True)
            or reference.

        """

        # Data is held in a list of three-element tuples (key,value,comment)
        # in cards. Actual header items can be looked up via the _lookup
        # dictionary. cards can also contain blank lines (all elements = '')
        # comments: key='COMMENTS', value=the comment, comment='', and
        # history:  key='HISTORY', value=the history, comment=''.
        self.cards = []
        self._lookup = {}

        if isinstance(head, Header):
            # Another Header
            if copy:
                self.cards = head.cards.copy()
                self._lookup = head._lookup.copy()
            else:
                self.cards = head.cards
                self._lookup = head._lookup

        elif isinstance(head, odict):
            # OrderedDict
            self.cards = []
            self.keys = {}
            for key, (value, comment) in head.items():
                key = Header._process_key(key)
                if key not in Header.SPECIAL_KEYWORDS:
                    self.cards.append((key, value, comment))
                    ukey = key.upper()
                    if ukey not in self._lookup:
                        self._lookup[ukey] = len(self.cards)-1
                    else:
                        raise ValueError(
                            ('key = {:s} appears more than once'
                             ' (case-insensitive)').format(ukey)
                        )

        elif isinstance(head, FITS_Header):
            # Build from an astropy.io.fits.HEADER
            for n, card in enumerate(head.cards):
                key = card.keyword
                self.cards.append((key,card.value,card.comment))

                ukey = key.upper()
                if ukey not in Header.SPECIAL_KEYWORDS:
                    if ukey not in self._lookup:
                        self._lookup[ukey] = n
                    else:
                        raise ValueError(
                            ('key = {:s} appears more than once'
                             ' (case-insensitive)').format(ukey)
                        )
        else:
            # A list of (key,value,comment) tuples.
            if head is None:
                head = []

            if copy:
                self.cards = head.copy()
            else:
                self.cards = head

            for n, (key, value, comment) in enumerate(head):
                ukey = key.upper()
                if ukey not in Header.SPECIAL_KEYWORDS:
                    if ukey not in self._lookup:
                        self._lookup[ukey] = n
                    else:
                        raise ValueError(
                            ('key = {:s} appears more than once'
                             ' (case-insensitive)').format(ukey)
                        )

        # Calculate pointers to where COMMENTS and HISTORY start
        # and end so we can add in at the right place.
        self._hstart = self._hstop = self._cstart = self._cstop = -1

        for n, (key, value, comment) in enumerate(self.cards):
            if key == 'COMMENT':
                if self._cstart == -1:
                    self._cstart = n
                self._cstop = n+1

            elif key == 'HISTORY':
                if self._hstart == -1:
                    self._hstart = n
                self._hstop = n+1

        if self._cstart == -1:
            self._cstart = self._cstop = len(self.cards)

        if self._hstart == -1:
            self._hstart = self._hstop = len(self.cards)

    @property
    def to_fits(self):
        """Returns the Header as an astropy.io.fits.Header
        which is useful when saving to disk
        """
        # Translate the cards to avoid warnings about long keys
        cards = []
        for key,value,comment in self.cards:
            if len(key) > 8:
                key = 'HIERARCH ' + key
            cards.append((key,value,comment))

        return FITS_Header(cards)

    @staticmethod
    def _process_key(key):
        """This converts keys of 8 characters or less to uppercase.
        Any leading 'HIERARCH ' is stripped"""
        if len(key) <= 8:
            key = key.upper()
        elif key.upper().startswith('HIERARCH '):
            key = key[9:]
        return key

    def __setitem__(self, key, item):
        """Sets the value of a header item using the [] form, e.g. `head['NREC'] =
        2345`.  Reserved keywords 'COMMENT', 'HISTORY' and '' (blank) are
        trapped.

        Arguments::

          key : string
             the name of the header item.

          item : tuple | string | int | float
             the (value,comment) pair or just the value (in which
             case the comment will be blank.

        """
        key = Header._process_key(key)

        if key in Header.SPECIAL_KEYWORDS:
            raise ValueError(
                "Keywords 'COMMENT', 'HISTORY' and '' are reserved"
            )

        elif isinstance(item,tuple):
            value, comment = item

        else:
            value, comment = item, None

        ukey = key.upper()
        if ukey in self._lookup:
            # Overwrite pre-existing value, re-covering the old comment if
            # no new one supplied
            comment = self.cards[self._lookup[ukey]][2] if \
                      comment is None else comment
            self.cards[self._lookup[ukey]] = (key,value,comment)
        else:
            # Insert new item before the history or comments start
            index = min(self._hstart,self._cstart)
            self.cards.insert(index, (key,value,comment))
            self._lookup[key] = index
            self._hstart += 1
            self._hstop += 1
            self._cstart += 1
            self._cstop += 1

    def __getitem__(self, key):
        """Returns the value associated with header item 'key' using the
        [] form, e.g. `nrec = head['NRECORD']`."""
        key = Header._process_key(key)
        return self.cards[self._lookup[key]][1]

    def get(self, key, default):
        """Returns the value associated with header item 'key', returning
        a default value if the header item does not exist."""
        key = Header._process_key(key)
        if key in self._lookup:
            return self.cards[self._lookup[key]][1]
        else:
            return default

    def get_full(self, key):
        """Returns with 2-element (value,comment) tuple for the given keyword"""
        key = Header._process_key(key)
        return self.cards[self._lookup[key]][1:]

    def get_comment(self, key):
        """Returns with the comment for the given keyword"""
        key = Header._process_key(key)
        return self.cards[self._lookup[key]][2]

    def add_comment(self, comment):
        """Adds a comment to the end of the 'COMMENT' section"""
        if isinstance(comment,str):
            if self._hstart >= self._cstop:
                self._hstart += 1
            if self._hstop >= self._cstop:
                self._hstop += 1
            self.cards.insert(self._cstop,('COMMENT',comment,''))
            self._cstop += 1
        else:
            raise ValueError('only string comments allowed')

    def add_history(self, history):
        """Adds a line of history to the end of the 'HISTORY' section"""
        if isinstance(history,str):
            if self._cstart >= self._hstop:
                self._cstart += 1
            if self._cstop >= self._hstop:
                self._cstop += 1
            self.cards.insert(self._hstop,('HISTORY',history,''))
            self._hstop += 1
        else:
            raise ValueError('only string history lines allowed')

    def copy(self):
        return Header(self, True)

    def __repr__(self):
        return 'Header(head={!r})'.format(self.cards)

    def __delitem__(self, key):
        """Deletes a particular header item when called e.g. as
        `del head[key]`."""

        key = Header._process_key(key)
        index = self._lookup[key]
        del self.cards[index]
        del self._lookup[key]

        # Now need to correct comment and history pointers
        # along with any later index pointers.
        self._hstart -= 1
        self._hstop -= 1
        self._cstart -= 1
        self._cstop -= 1
        for key,ind in self._lookup.items():
            if ind > index:
                self._lookup[key] = ind-1

    def __contains__(self, key):
        """Defines the operation 'in' for a Header. Look up whether
        a given key is in the header"""
        key = Header._process_key(key)
        return key in self._lookup

    def update(self, head):
        """Update the Header with the contents of another header (excluding comments
and history)"""
        nadded = 0
        for key,value,comment in head.cards:
            if key.upper() not in Header.SPECIAL_KEYWORDS:
                self[key] = (value,comment)
