"""Defines a class for feeding out images via a context manager-come-iterator
to provide a uniform data access.

"""

from .ccd import rfits
from . import ucam

__all__ = ('Spooler',)

class Spooler:

    """A common requirement is the need to loop through a stack of images. The aim
    of this class is to provide uniform access via an iterable context
    manager. This allows it to cover 5 cases in the same way: access to raw
    data via a server or a local disk file, whether from ULTRACAM/SPEC or
    HiPERCAM, and access to lists of files.

    """

    # used to indicate which instrument is involved.
    ULTRA = 1
    HIPER = 2

    def __init__(self, run=None, server=False, inst=Spooler.HIPER, flt=True,
                 flist=None):
        """Sets up info for establishing a source of data frames. Possibilities are
        (1) via a list of CCD or MCCD files, (2) via a server and (3) from a
        single raw data file on a local disk. It copes with data from
        ULTRACAM/SPEC or HiPERCAM.

        Arguments::

           run : (string | None)
              If not None, this must contain the run number as in 'run003' to
              access, whether via a server or the local disk.

           server : (bool | None)
              If run is not None, then this must be a bool, True for file
              access via a server, False for file accessed from local disk.

           inst : (string)
              If run is not None, this must be set to either Spooler.ULTRA or
              Spooler.HIPER for ULTRACAM/SPEC vs HiPERCAM as the data source.

           flt : (bool)
              If run is not None, flt=True/False determines whether the raw
              data are converted to 4-bytes floats on input or not.

           flist : (string | None)
              If run is None, this must be set to the name of a list of files.

        """

        self.run = run
        self.server = server
        self.inst = inst
        self.flt = flt
        self.flist = flist

        if self.run is not None:
            if self.inst == Spooler.ULTRA:
                self._iter = ucam.Rdata(self.run, flt=self.flt,
                                        server=self.server)
            else:
                raise ValueError('Spooler.__init__: unrecognised argument combination')
        else:
            raise ValueError('Spooler.__init__: unrecognised argument combination')

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
