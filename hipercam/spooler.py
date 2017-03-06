"""Defines classes and methods for reading data of varying formats and
for iterating through multiple images from a variety of sources.


"""

from abc import ABC, abstractmethod
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

class SpoolerBase(ABC):

    """A common requirement is the need to loop through a stack of images. With a
    variety of possible data sources, one requires handling of multiple
    possible ways of accessing the data. The aim of this class is to provide
    uniform access via an iterable context manager. It is written as an abstract
    class 

    """

    def __enter__(self):
        return self

    @abstractmethod
    def __exit__(self, *args):
        pass

    def __iter__(self):
        return self

    @abstractmethod
    def __next__(self):
        pass

class UcamDiskSpool(SpoolerBase):

    """Provides an iterable context manager to loop through frames within
    a raw ULTRACAM or ULTRASPEC disk file.
    """

    def __init__(self, run, first=1, flt=False):
        """Attaches the UcamDiskSpool to a run.

        Arguments::

           run : (string)

              The run number, e.g. 'run003' or 'data/run004'. 

           first : (int)
              The first frame to access.

           flt : (bool)
              If True, a conversion to 4-byte floats on input will be attempted
              if the data are of integer form.

        """
        self._iter = ucam.Rdata(run, first, flt, False)

    def __exit__(self, *args):
        self._iter.__exit__(args)

    def __next__(self):
        return self._iter.__next__()

class UcamServSpool(SpoolerBase):

    """Provides an iterable context manager to loop through frames within
    a raw ULTRACAM or ULTRASPEC raw file served from the ATC FileServer
    """

    def __init__(self, run, first=1, flt=False):
        """Attaches the UcamDiskSpool to a run.

        Arguments::

           run : (string)

              The run number, e.g. 'run003' or 'data/run004'. 

           first : (int)
              The first frame to access.

           flt : (bool)
              If True, a conversion to 4-byte floats on input will be attempted
              if the data are of integer form.

        """
        self._iter = ucam.Rdata(run, first, flt, True)

    def __exit__(self, *args):
        self._iter.__exit__(args)

    def __next__(self):
        return self._iter.__next__()


class HcamListSpool(SpoolerBase):

    """Provides an iterable context manager to loop through frames within
    a list of HiPERCAM .hcm files.
    """

    def __init__(self, lname):
        """Attaches the :class:`HcamListSpool` to a list of files

        Arguments::

           lname : (string)
              Name of a list of HiPERCAM files. # at the start of a line is recognized
              as a comment flag.

        """
        self._iter = open(lname)

    def __exit__(self, *args):
        self._iter.close()

    def __next__(self):
        for fname in self._iter:
            if not fname.startswith('#'):
                break
        else:
            raise StopIteration


def data_source(inst, resource, server=None, flist=none, first=1, flt=True):
    """Returns a context manager needed to run through a set of exposures.
    This is basically a wrapper around the various context managers that
    hook off the SpoolerBase class.

    Arguments::

       inst : (string)
          Instrument name. Current choices: 'ULTRA' for ULTRACAM/SPEC, 'HIPER'
          for 'HiPERCAM'.

       resource : (string)
          File name. Either a run number for ULTRA or HIPER and server==True,
          of a file with a list of files for HIPER and server==False.

       server : (bool | None)
          If server == True, access via a server is expected. If False, a local
          disk file is assumed. If None, then flist should be set to True.

       flist : (bool | None)
          If server is None, then flist should be set to True, otherwise it
          should be set to None or False.

       first : (int)
          If a raw disk file is being read, either locally or via a server,
          this parameter sets where to start in the file.

       flt : (bool)
          If True, convert to 32-bit floats on input.

    Returns::

       A :class:`SpoolerBase` object that can appear inside a "with" 
       statement, e.g.

         with data_source('ULTRA', 'run003') as dsource:
            for frame in dsource:

    Exceptions::

       A ValueError will be raised if the inputs are not recognised or
       conflict.
    """

    if inst == 'ULTRA':
        if server is None:
            raise ValueError(
                'spooler.data_source: for inst = "ULTRA", server must be True or False')
        elif server:
            return UcamServSpool(resource,first,flt)
        else:
            return UcamDiskSpool(resource,first,flt)

    elif inst == 'HIPER':
        if server is None:
            if flist is None or not flist:
                raise ValueError(
                    'spooler.data_source: for inst = "HIPER" and server is None, flist must be True')
            else:
                return HcamListSpool(resource)
        elif server:
                raise ValueError(
                    'spooler.data_source: HiPERCAM server not implemented yet'
                )
        else:
                raise ValueError(
                    'spooler.data_source: HiPERCAM raw data nor defined yet'
                )
    else:
        raise ValueError(
            'spooler.data_source: inst = {!s} not recognised'.format(inst)
        )

