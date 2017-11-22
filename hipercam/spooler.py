"""Classes and methods to support reading data in multiple formats and for iterating
through multiple images from a variety of sources in a standardised manner.

For instance this piece of code::

   with spooler.data_source(source, resource, first) as spool:
       for mccd in spool:
           ...

creates an object 'spool' that then acts a source of :class:`MCCD` objects. The
:meth:`data_source` method returns this for a variety of possible sources of data
such as raw HiPERCAM and ULTRACAM files and a list of hcm files. Important methods
are :meth:`data_source` and :meth:`get_ccd_pars`. You may also want to look at
classes such as :class:`HcamListSpool` if you are interested in accessing data in a
specific manner.
"""

import time
from abc import ABC, abstractmethod
from collections import OrderedDict

from astropy.io import fits
from . import ccd
from . import ucam
from . import hcam
from . import core
from . import utils

__all__ = (
    'SpoolerBase', 'data_source', 'rhcam', 'UcamServSpool',
    'UcamDiskSpool', 'HcamListSpool', 'get_ccd_pars',
    'hang_about', 'HcamServSpool', 'HcamDiskSpool',
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
    a list of HiPERCAM .hcm files, either a file listing the files, or a Python
    list. Use as::

        with spooler.HcamListSpool(flist) as spool:
            for mccd in spool:
                ...

    to get :class:`MCCD`  objects or::


        with spooler.HcamListSpool(flist,cnam) as spool:
            for ccd in spool:
                ...

    to get :class:`CCD` objects.
    """

    def __init__(self, lname, cnam=None):
        """Attaches the :class:`HcamListSpool` to a list of files

        Arguments::

           lname : string or list
              Name of a file list of HiPERCAM files, or a Python list of file
              names (more generally any iterable array of names). Blank
              entries or any starting with '#' are ignored.

           cnam  : string or None
              CCD label if you want to return individual CCDs rather than
              MCCDs. This is used in 'combine' to save memory.
        """
        if isinstance(lname, string):
            self._iter = open(lname)
            self._close_at_end = True
        else:
            self._iter = lname
            self._close_at_end = False

        self.cnam = cnam

    def __exit__(self, *args):
        if self._close_at_end:
            # close in the case of file lists
            self._iter.close()

    def __next__(self):
        # returns next image from a list. adds in the file name
        # as a header parameter
        for fname in self._iter:
            if not fname.startswith('#') and not fname.isspace():
                break
        else:
            raise StopIteration

        if self.cnam is None:
            mccd = ccd.MCCD.read(utils.add_extension(fname.strip(), core.HCAM))
            mccd.head['FILENAME'] = fname
            return mccd
        else:
            ccd = ccd.CCD.read(
                utils.add_extension(fname.strip(), core.HCAM), self.cnam
            )
            ccd.head['FILENAME'] = fname
            return ccd

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


def data_source(source, resource, first=1):
    """Returns a context manager needed to run through a set of exposures.
    This is basically a wrapper around the various context managers that
    hook off the SpoolerBase class.

    Arguments::

       source    : (string)

          Data source. Options are 'hl', 'hs', 'ul', 'us', 'hf'. The leading
          ('h' | 'u') indicates ULTRA(CAM|SPEC), or HiPERCAM. The trailing
          ('l' | 's' | 'f') refers to access through a local file for 'l'
          (.fits for HiPERCAM or .dat/.xml for ULTRA(CAM|SPEC), a server for
          's' (either the ATC fileserver for ULTRA(CAM|SPEC) or Stu
          Littlefair's HiPERCAM server), or from a list of fils in HiPERCAM's
          FITS-based hcm format.

       resource : (string)
          File name. A run number if source=('?l'|'?s') or a file list
          for source='hl'.

       first    : (int)
          If a raw disk file is being read, either locally or via a server,
          this parameter sets where to start in the file. 0 to always try
          to get the last. See also 'hang_about' in this case. This parameter
          is ignored if source=='hf'.

    Returns::

       A :class:`SpoolerBase` object that can appear inside a "with"
       statement, e.g.

         >> with data_source('us', 'run003') as dsource:
         >>    for frame in dsource:
         >>       ... do something with 'frame'  
    """

    if source == 'us':
        return UcamServSpool(resource, first)
    elif source == 'ul':
        return UcamDiskSpool(resource, first)
    elif source == 'hs':
        return HcamServSpool(resource, first)
    elif source == 'hl':
        return HcamDiskSpool(resource, first)
    elif source == 'hf':
        return HcamListSpool(resource)
    else:
        raise ValueError(
            '{!s} is not a recognised data source'.format(source)
        )

def get_ccd_pars(source, resource):
    """Returns a list of tuples of CCD labels, and maximum dimensions,
    i.e. (label, nxmax, nymax) for each CCD pointed at by the input
    parameters.

    Arguments::

        source   : (string)
           Data source. See 'data_source' for details.

       resource : (string)
          File name. Either a run number or a file list. Again, see
          'data_source' for details.

    Returns with a list of tuples with the information outlined above. In the
    case of file lists these are extracted from the first file of the list only.

    """

    if source.startswith('u'):

        # ULTRA(CAM|SPEC)
        rhead = ucam.Rhead(resource)
        if rhead.instrument == 'ULTRACAM':
            # ULTRACAM raw data file: fixed data
            return OrderedDict(
                (('1',(1080,1032)), ('2',(1080,1032)), ('3',(1080,1032)))
            )

        elif rhead.instrument == 'ULTRASPEC':
            # ULTRASPEC raw data file: fixed data
            return OrderedDict(('1',(1056,1072)))

        else:
            raise ValueError(
                'instrument = {:s} not supported'.format(rhead.instrument)
            )

    elif source.startswith('h'):
        if source.endswith('f'):
            # file list: we access the first file of the list to read the key
            # and dimension info on the assumption that it is the same, for
            # all files.
            with open(resource) as flp:
                for fname in flp:
                    if not fname.startswith('#') and not fname.isspace():
                        break
                else:
                    raise ValueError(
                        'failed to find any file names in {:s}'.format(resource)
                        )
            return ccd.get_ccd_info(
                utils.add_extension(fname.strip(), core.HCAM))

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

def hang_about(obj, twait, tmax, total_time):
    """
    Carries out some standard actions when we loop through frames which are
    common to rtplot, reduce and grab. This is a case of seeing whether we
    want to try again for a frame or a time that may have arrived while we wait.

    Arguments::

         obj   : object or None
            something, e.g. an MCCD or a Time, it does not matter what,
            because the only thing that is tested is whether it is 'None'. If
            it is, that is a trigger to wait a bit before trying again

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

    if obj is None:

        if tmax < total_time + twait:
            if tmax > 0.:
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
            give_up, try_again = False, True

    else:
        # we have a frame, reset the total time waited
        # set flags to show we don't want to give up or
        # try again
        total_time = 0
        give_up, try_again = False, False

    return (give_up, try_again, total_time)

