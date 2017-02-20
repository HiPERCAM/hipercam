"""
Defines a class for feeding out image via a context manager-come-iterator to provide
a uniform data access
"""

from .ccd import rfits
from . import ucam

__all__ = ('Feeder',)

class Feeder:

    """Uniform interface context manager-come-iterator to allow a list of CCD /
    MCCD FITS files to be read in a for loop from either a file list or raw
    data files from any of ULTRACAM / ULTRASPEC / HiPERCAM.

    """

    # used to indicate which instrument is involved.
    ULTRA = 1
    HIPER = 2

    def __init__(self, run=None, server=False, inst=HIPER, flist=None): 
        """Sets up info for establishing a source of data frames. Possibilities are
        (1) via a list of CCD or MCCD files, (2) via a server and (3) from a
        single raw data file on a local disk. It copes with data from ULTRACAM/SPEC or
        HiPERCAM.

        Arguments::

           run : (string | None)
              If nor None, this must contain the run number as in 'run003' to
              access, whether via a server or the local disk.

           server : (bool | None)
              If run is not None, then this must be a bool, True for file
              access via a server, False for file accessed from local disk.

           inst : (string)
              either Feeder.ULTRA or Feeder.HIPER for ULTRACAM/SPEC vs
              HiPERCAM as the instrument source.

           flist : (string | None)
              If run is None, this must be set to the name of a list of files.

        """

        self.run = run
        self.server = server
        self.inst = inst
        self.flist = flist

        if self.run is not None:
            if not self.server:
                if self.inst == Feeder.ULTRA:
                    self._iter = ucam.Rdata(self.run)
                else:
                    raise ValueError('Feeder.__init__: unrecognised argument combination')
            else:
                raise ValueError('Feeder.__init__: unrecognised argument combination')
        else:
            raise ValueError('Feeder.__init__: unrecognised argument combination')

    # next two define the actions associated with the context manager
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._iter.__exit__(args)

    # next two are to make this an iterator
    def __iter__(self):
        return self

    def __next__(self):
        return self._iter.__next__()
