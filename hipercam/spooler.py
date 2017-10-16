"""Defines classes and methods for reading data of varying formats and
for iterating through multiple images from a variety of sources in a
standardised manner.


"""

import time
from abc import ABC, abstractmethod
from collections import OrderedDict

from astropy.io import fits
from . import ccd
from . import ucam
from . import hcam
from . import core

__all__ = (
    'SpoolerBase', 'data_source', 'rhcam', 'UcamServSpool',
    'UcamDiskSpool', 'HcamListSpool', 'get_ccd_pars',
    'hang_about'
)

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
            return ccd.CCD.rhdul(hdul)
        elif htype == 'MCCD':
            return ccd.MCCD.rhdul(hdul)
        else:
            raise ValueError(
                'Could not find keyword "HIPERCAM" in primary header of file = {:s}'.format(
            fname)
            )

class SpoolerBase(ABC):

    """A common requirement is the need to loop through a stack of
    images. With a variety of possible data sources, one requires handling of
    multiple possible ways of accessing the data. The aim of this class is to
    provide uniform access via an iterable context manager. It is written as
    an abstract class.
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

    def __init__(self, run, first=1):
        """Attaches the UcamDiskSpool to a run.

        Arguments::

           run : (string)

              The run number, e.g. 'run003' or 'data/run004'. 

           first : (int)
              The first frame to access.

        """
        self._iter = ucam.Rdata(run, first, False)

    def __exit__(self, *args):
        self._iter.__exit__(args)

    def __next__(self):
        return self._iter.__next__()

class UcamServSpool(SpoolerBase):

    """Provides an iterable context manager to loop through frames within
    a raw ULTRACAM or ULTRASPEC raw file served from the ATC FileServer
    """

    def __init__(self, run, first=1):
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
        self._iter = ucam.Rdata(run, first, True)

    def __exit__(self, *args):
        self._iter.__exit__(args)

    def __next__(self):
        return self._iter.__next__()

class HcamDiskSpool(SpoolerBase):

    """Provides an iterable context manager to loop through frames within
    a raw HiPERCAM file.
    """

    def __init__(self, run, first=1):
        """Attaches the HcamDiskSpool to a run.

        Arguments::

           run : (string)

              The run number, e.g. 'run003' or 'data/run004'. 

           first : (int)
              The first frame to access.

        """
        self._iter = hcam.Rdata(run, first, False)

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
        return ccd.MCCD.rfits(core.add_extension(fname.strip(), core.HCAM))

class HcamServSpool(SpoolerBase):

    """Provides an iterable context manager to loop through frames within
    a raw HiPERCAM file served from Stu's FileServer
    """

    def __init__(self, run, first=1):
        """Attaches the HcamServSpool to a run.

        Arguments::

           run : (string)

              The run number, e.g. 'run003' or 'data/run004'.

           first : (int)
              The first frame to access, 0 to always try to return the last complete one.

        """
        self._iter = hcam.Rdata(run, first, True)

    def __exit__(self, *args):
        self._iter.__exit__(args)

    def __next__(self):
        return self._iter.__next__()


def data_source(inst, resource, flist, server, first):
    """Returns a context manager needed to run through a set of exposures.
    This is basically a wrapper around the various context managers that
    hook off the SpoolerBase class.

    Arguments::

       inst     : (string)
          Instrument name. Current choices: 'ULTRA' for ULTRACAM/SPEC, 'HIPER'
          for 'HiPERCAM'. 'HIPER' is the one to use if you are looking at a list
          of HiPERCAM-format 'hcm' files.

       resource : (string)
          File name. Either a run number for ULTRA or HIPER or a file list
          if flist == True.

       flist    : (bool)
          True if 'resource' is a list of files.

       server   : (bool)
          If flist == False, then 'server' controls whether access will be attempted via
          a server. If not, a local disk file is assumed. Note that this parameter is
          ignored if flist == True.

       first    : (int)
          If a raw disk file is being read, either locally or via a server,
          this parameter sets where to start in the file. 0 for latest.

    Returns::

       A :class:`SpoolerBase` object that can appear inside a "with" 
       statement, e.g.

         with data_source('ULTRA', 'run003') as dsource:
            for frame in dsource:
    """

    if inst == 'ULTRA':
        if flist:
            raise NotImplementedError(
                'spooler.data_source: file lists have not been implemented for inst = "ULTRA"'
                )

        if server:
            return UcamServSpool(resource,first)
        else:
            return UcamDiskSpool(resource,first)

    elif inst == 'HIPER':
        if flist:
            return HcamListSpool(resource)
        elif server:
            return HcamServSpool(resource,first)
        else:
            return HcamDiskSpool(resource,first)

    else:
        raise ValueError(
            'spooler.data_source: inst = {!s} not recognised'.format(inst)
        )

def get_ccd_pars(inst, resource, flist):
    """Returns a list of tuples of CCD labels, and maximum dimensions, i.e. (label, nxmax, nymax)
    for each CCD pointed at by the input parameters.

    Arguments::

       inst     : (string)
          Instrument name. Current choices: 'ULTRA' for ULTRACAM/SPEC, 'HIPER'
          for 'HiPERCAM'. 'HIPER' is the one to use if you are looking at a list
          of HiPERCAM-format 'hcm' files.

       resource : (string)
          File name. Either a run number for ULTRA or HIPER or a file list
          if flist == True.

       flist    : (bool)
          True if 'resource' is a list of files.

    Returns with a list of tuples with the information outlined above. In the
    case of file lists these are extracted from the first file of the list only.
    """

    if inst == 'ULTRA':
        if flist:
            raise NotImplementedError(
                'spooler.get_ccd_pars: file lists have not been implemented for inst = "ULTRA"'
                )
        rhead = ucam.Rhead(resource)
        if rhead.instrument == 'ULTRACAM':
            # ULTRACAM raw data file: fixed data
            return OrderedDict((('1',(1080,1032)), ('2',(1080,1032)), ('3',(1080,1032))))

        elif rhead.instrument == 'ULTRASPEC':
            # ULTRASPEC raw data file: fixed data
            return OrderedDict(('1',(1056,1072)))

        else:
            raise ValueError(
                'spooler.get_ccd_pars: instrument {:s} not supported'.format(rhead.instrument)
                )

    elif inst == 'HIPER':
        if flist:
            # file list: we access the first file of the list to read the key and dimension
            # info on the assumption that it is the same, for all files.
            with open(resource) as flp:
                for fname in flp:
                    if not fname.startswith('#'):
                        break
                else:
                    raise ValueError(
                        'spooler.get_ccd_pars: failed to find any file names in {:s}'.format(resource)
                        )
            return ccd.get_ccd_info(core.add_extension(fname.strip(), core.HCAM))

        else:
            # HiPERCAM raw data file: fixed data
            return OrderedDict(
                (
                    ('1',(hcam.HCM_NXTOT, hcam.HCM_NYTOT)),
                    ('2',(hcam.HCM_NXTOT, hcam.HCM_NYTOT)),
                    ('3',(hcam.HCM_NXTOT, hcam.HCM_NYTOT)),
                    ('4',(hcam.HCM_NXTOT, hcam.HCM_NYTOT)),
                    ('5',(hcam.HCM_NXTOT, hcam.HCM_NYTOT))
                    )
                )

def hang_about(mccd, twait, tmax, total_time):
    """
    Carries out some standard actions when we loop through
    frames which are common to rtplot, reduce and grab. This
    is a case of seeing whether we want to try again for a
    frame that may have arrived while we wait.

    Arguments::

         mccd   : (MCCD | none)
            the frame just read. If None, that is a trigger to wait
            a bit before trying again

         twait  : (float)
            how long to wait each time, seconds

         tmax   : (float)
            maximum time to wait, seconds

         total_time : (float)
            total time waited. Should be initialised to zero

    Returns: (give_up, try_again, total_time)

       give_up == True   ==> stop trying
       try_again == True ==> have another go

    If both are False, that indicates success and we carry on.
    """

    if mccd is None:

        if tmax < total_time + twait:
            print(
                ' ** last frame unchanged for {:.1f} sec.'
                ' cf tmax = {:.1f}; will wait no more'.format(
                    total_time, tmax)
            )
            give_up, try_again = True, False

        else:

            if total_time == 0:
                # add a blank line the first time round
                print()

            print(
                ' ** last frame unchanged for {:.1f} sec.'
                ' cf tmax = {:.1f}; will wait another'
                ' twait = {:.1f} sec.'.format(
                    total_time, tmax, twait)
            )

            # pause
            time.sleep(twait)

            # increment the waiting time
            total_time += twait

            # have another go
            give_up, try_gain = False, True

    else:
        # we have a frame, reset the total time waited
        # set flags to show we don't want to give up or
        # try again
        total_time = 0
        give_up, try_again = False, False

    return (give_up, try_again, total_time)

