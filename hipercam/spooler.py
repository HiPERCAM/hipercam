"""Defines classes and methods for reading data of varying formats and
for iterating through multiple images from a variety of sources.


"""

from astropy.io import fits
from .ccd import CCD, MCCD
from . import ucam

__all__ = ('Spooler', 'data_source', 'rhcam')

def rhcam(fname):
    """Reads a HiPERCAM file containing either CCD or MCCD data and returns one or
    the other.

    Argument::

       fname : (string)

          Name of file. This should be a FITS file with a particular
          arrangement of HDUs. The keyword HIPERCAM should be present in the
          primary HDU's header and set to either CCD for single 'CCD' data or
          'MCCD' for multiple CCDs. A ValueError will be raised if this is not
          the case. See the docs on CCD and MCCD for more.

    """

    # Read HDU list
    with fits.open(fname) as hdul:
        htype = hdul[0].header['HIPERCAM']
        if htype == 'CCD':
            return CCD.rhdul(hdul)
        elif htype == 'MCCD':
            return MCCD.rhdul(hdul)
        else:
            raise ValueError(
                'Could not find keyword "HIPERCAM" in primary header of file = {:s}'.format(
            fname)
            )

def data_source(inst, server=None, flist=True):
    """Helper routine to return the data source needed by the :class:`Spooler`
    class given an instrument name and whether access is via a server (else
    disk) or a file list.

    Arguments::

       inst : (string)
          Instrument name. Current choices: 'ULTRA' for ULTRACAM/SPEC, 'HIPER'
          for 'HiPERCAM'.

       server : (bool | None)
          If server == True, access via a server is expected. If False, a local
          disk file is assumed. If None, then flist should be set to True.

       flist : (bool | None)
          If server is None, then flist should be set to True, otherwise it
          should be set to None or False. 

    Returns::

       An integer corresponding to one of the :class:`Spooler` class attributes
       representing the supported data sources.

    Exceptions::

       A ValueError will be raised if the inputs are not recognised or
       conflict.
    """

    if inst == 'ULTRA':
        if server is None:
            raise ValueError(
                'spooler.data_source: for inst = "ULTRA", server must be True or False')
        elif server:
            return Spooler.ULTRA_SERV
        else:
            return Spooler.ULTRA_DISK

    elif inst == 'HIPER':
        if server is None:
            if flist is None or not flist:
                raise ValueError(
                    'spooler.data_source: for inst = "HIPER" and server is None, flist must be True')
            else:
                return Spooler.HIPER_LIST

        elif server:
            return Spooler.HIPER_SERV
        else:
            return Spooler.HIPER_DISK
    else:
        raise ValueError(
            'spooler.data_source: inst = {!s} not recognised'.format(inst)
        )

class Spooler:

    """A common requirement is the need to loop through a stack of images. With a
    variety of possible data sources, one requires handling of multiple
    possible ways of accessing the data. The aim of this class is to provide
    uniform access via an iterable context manager.

    """

    # supported types of data access
    HIPER_DISK = 1 # HiPERCAM raw data from a local disk file
    HIPER_SERV = 2 # HiPERCAM raw data from a server
    HIPER_LIST = 3 # HiPERCAM from list of hcm files
    ULTRA_DISK = 4 # ULTRACAM/SPEC raw data from a local disk file
    ULTRA_SERV = 5 # ULTRACAM/SPEC raw data from a server

    def __init__(self, ident, source=HIPER_SERV, first=1, flt=False):
        """Sets up info for establishing a source of data frames.

        Arguments::

           ident : (string)

              An identifier of the data source whose nature depends on the
              value of `source`. For instance if source == HIPER_DISK, this
              should be a run number e.g. 'run003' or 'data/run004'. If source
              == HIPER_LIST then it should be the name of a file list.

           source : (int)
              Data source. The possibilities are defined by a set of class
              attributes such as `HIPER_DISK` and `ULTRA_SERV`. See below
              for a full listing of their meanings. You will need to specify
              the class as in `Spooler.ULTRA_DISK`.

           first : (int)
              The first frame to access (ignored for the file list sources on
              the basis that lists can be tuned to suit usage)

           flt : (bool)
              If True, a conversion to 4-byte floats on input will be attempted
              if the data are of integer form.

        Current source options and their implications for the meaning of
        `ident`::

           HIPER_DISK : HiPERCAM raw data from a local disk file. `ident`
                        should specify the run in this case.

           HIPER_SERV : HiPERCAM raw data from a server. `ident` should
                        specify the run; the server name will be grabbed from
                        an environment variable `HIPERCAM_DEFAULT_URL`

           HIPER_LIST : HiPERCAM from list of hcm files. `ident` should
                        specify the name of a list of .hcm files.

           ULTRA_DISK : ULTRACAM/SPEC raw data from a local disk file. `ident`
                        should specify the run.

           ULTRA_SERV : ULTRACAM/SPEC raw data from a server. `ident` should
                        specify the run; the server name will be grabbed from
                        an environment variable `ULTRACAM_DEFAULT_URL`

        """

        if source == Spooler.HIPER_DISK:
            raise ValueError(
                'Spooler.__init__: source == HIPER_DISK not yet supported')

        elif source == Spooler.HIPER_SERV:
            raise ValueError(
                'Spooler.__init__: source == HIPER_SERV not yet supported')

        elif source == Spooler.ULTRA_DISK:
            self._iter = ucam.Rdata(run, first, flt, False)

        elif source == Spooler.ULTRA_SERV:
            self._iter = ucam.Rdata(run, first, flt, True)

        elif source == Spooler.HIPER_LIST:
            self._iter = open(ident)

        else:
            raise ValueError(
                'Spooler.__init__: source == {:d} not recognised'.format(source)
            )

        self.source = source

    # next two define the actions associated with the context manager
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._iter.__exit__(args)

    # next two are to make this an iterator
    def __iter__(self):
        return self

    def __next__(self):
        if self.source == Spooler.HIPER_LIST:
            for fname in self._iter:
                if not fname.startswith('#'):
                    break
            else:
                raise StopIteration

            return rhcam(fname.strip())
        else:
            return self._iter.__next__()
