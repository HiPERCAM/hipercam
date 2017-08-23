"""
Section for handling the raw HiPERCAM data files. These are FITS files
with a single HDU containing extensive header information and the data in
a "FITS cube" of unusual nature.
"""

import struct
import warnings
import numpy as np
import urllib.request

from numpy.lib.stride_tricks import as_strided
import fitsio
from astropy.io import fits

from hipercam import (CCD, Group, MCCD, Windat, Window, HCM_NXTOT, HCM_NYTOT)
from hipercam import (HRAW, add_extension) 
from .herrors import *

__all__ = ['Rhead', 'Rdata']

class Rhead:
    """Reads an interprets header information from a HiPERCAM run and
    generates the window formats needed to interpret the data. The file is
    kept open while the Rhead exists.

    The following attributes are set::

      cheads  : (list of astropy.io.fits.Header)
         the headers for each CCD which contain any info needed per CCD, above
         all the 'skip' numbers.

      ffile   : (fitsio.FITS)
         the opened file

      fname   : (string)
         the file name

      header  : (fitsio??)
         the header of the FITS file

      thead   : (astropy.io.fits.Header)
         the top-level header

      windows : (list of Windows)
         the windows are the same for each CCD, so these are just one set
         in EFGH order for window 1 and then the same for window 2 if there
         is one, so either 4 or 8 Windows.

      mode    : (string)
         the readout mode used.

    """

    def __init__(self, fname, server=False, full=True):
        """Opens a fits file containing HiPERCAM data, reads and interprets
        the header of the primary HDU.

        Arguments::

           fname  : (string)
              file name

           server : (bool)
              True to access the data from a server. It uses a URL set in an
              environment variable "HIPERCAM_DEFAULT_URL" in this instance.

           full   : (bool)
              Flag controlling the amount of header detail to gather (as a
              time saver) True for detail, False for the bare minimum, the
              latter might be useful for large numbers of very small format
              images where you don't to waste effort.
        """

        # store the file name
        self.fname = fname

        # open the file with fitsio
        self.ffile = fitsio.FITS(add_extension(fname,HRAW))

        # read the header
        self.header = self.ffile[0].read_header()

        # create the top-level header
        self.thead = fits.Header()

        # set the mode, one or two windows per quadrant, drift etc.
        # This is essential.
        self.mode = self.header['ESO DET READ CURNAME']
        self.thead['MODE'] = (self.mode,'HiPERCAM readout mode')

        # tuples of fixed data for each quadrant
        LLX = (1, 1025, 1025, 1)
        LLY = (1, 1, 521, 521)
        READOUT_Y = (1, 1, 1040, 1040)
        OFFSET = (1, 1, -1, -1)
        ADD_YSIZE = (0, 0, 1, 1)
        QUAD = ('E', 'F', 'G', 'H')

        # extract the binning factors
        xbin = self.header['ESO DET BINX1']
        ybin = self.header['ESO DET BINY1']

        # extract the data to build the first 4 windows
        if self.mode.startswith('FullFrame'):
            nx = 1024 // xbin
            ny = 520 // ybin
            llxs = LLX
            llys = LLY

        elif self.mode.startswith('OneWindow') or \
             self.mode.startswith('TwoWindows'):
            nx = self.header['ESO DET WIN1 NX'] // xbin
            ny = self.header['ESO DET WIN1 NY'] // ybin
            llxs = [llx + self.header['ESO DET WIN1 XS{}'.format(quad)]
                    for llx, quad in zip(LLX, QUAD)]
            llys = [
                readout_y + offset*self.header['ESO DET WIN1 YS'] -
                add_ysize*self.header['ESO DET WIN1 NY']
                for readout_y, offset, add_ysize in zip(READOUT_Y, OFFSET, ADD_YSIZE)
                ]

        elif self.mode.startswith('Drift'):
            nx = self.header['ESO DET DRWIN NX'] // xbin
            ny = self.header['ESO DET DRWIN NY'] // ybin
            llxs = [llx + self.header['ESO DET DRWIN XS{}'.format(quad)]
                    for llx, quad in zip(LLX, QUAD)]
            llys = [
                readout_y + offset*self.header['ESO DET DRWIN1 YS'] -
                add_ysize*self.header['ESO DET DRWIN NY']
                for readout_y, offset, add_ysize in zip(READOUT_Y, OFFSET, ADD_YSIZE)
                ]

        else:
            msg = 'mode {} not currently supported'.format(self.mode)
            raise ValueError(msg)

        # add in the first four Windows
        self.windows = []
        for llx, lly in zip(llxs, llys):
            self.windows.append(Window(llx, lly, nx, ny, xbin, ybin))

        if self.mode.startswith('TwoWindows'):
            # add extra four for this mode
            nx = self.header['ESO DET WIN2 NX'] // xbin
            ny = self.header['ESO DET WIN2 NY'] // ybin
            llxs = [llx + self.header['ESO DET WIN2 XS{}'.format(quad)]
                    for llx, quad in zip(LLX, QUAD)]
            llys = [
                readout_y + offset*self.header['ESO DET WIN2 YS'] -
                add_ysize*self.header['ESO DET WIN2 NY']
                for readout_y, offset, add_ysize in zip(READOUT_Y, OFFSET, ADD_YSIZE)
                ]

            for llx, lly in zip(llxs, llys):
                self.windows.append(Window(llx, lly, nx, ny, xbin, ybin))

        # Build (more) header info
        if 'DATE' in self.header:
            self.thead['DATE'] = (self.header.get('DATE'), self.header.get_comment('DATE'))
        if full and 'ESO DET GPS' in self.header:
            self.thead['GPS'] = (self.header.get('ESO DET GPS'), self.header.get_comment('ESO DET GPS'))

        # Header per CCD
        self.cheads = []
        for n in range(1,6):
            chead = fits.Header()

            # Essential items
            hnam = 'ESO DET NSKIP{:d}'.format(n)
            if hnam in self.header:
                chead['NSKIP'] = (self.header.get(hnam), self.header.get_comment(hnam))

            # Nice-if-you-can-get-them items
            if full:
                # whether this CCD has gone through a reflection
                hnam = 'ESO DET REFLECT{:d}'.format(n)
                if hnam in self.header:
                    chead['REFLECT'] = (self.header.get(hnam), 'is image reflected')

                # readout noise
                hnam = 'ESO DET CHIP{:d} RON'.format(n)
                if hnam in self.header:
                    chead['RONOISE'] = (self.header.get(hnam), self.header.get_comment(hnam))

                # gain
                hnam = 'ESO DET CHIP{:d} GAIN'.format(n)
                if hnam in self.header:
                    chead['GAIN'] = (self.header.get(hnam), self.header.get_comment(hnam))

            self.cheads.append(chead)

        # end of constructor / initialiser

    def __del__(self):
        """Destructor closes the file"""
        self.ffile.close()

    def npix(self):
        """
        Returns number of (binned) pixels per CCD
        """
        np = 0
        for win in self.windows:
            np += win.nx*win.ny
        return np

class Rdata (Rhead):
    """Callable, iterable object to represent HiPERCAM raw data files.

    The idea is to instantiate an Rdata from a HiPERCAM FITS file and then the
    object generated can be used to deliver MCCD objects by specifying a frame
    number e.g.::

      rdat = Rdata('run045.fits')
      fr10 = rdat(10)
      fr11 = rdat()

    reads frame number 10 and then 11 in from 'run045.fits', or sequentially::

      for mccd in Rdata('run045.fits'):
         print('nccd = ',len(mccd))

    which just prints out the number of CCDs from every exposure in the file
    (=5 in all cases),

    """

    def __init__(self, fname, nframe=1, server=False):
        """Connects to a raw HiPERCAM FITS file for reading. The file is kept
        open.  The Rdata object can then generate MCCD objects through being
        called as a function or iterator.

        Arguments::

           fname : (string)
              run name, e.g. 'run036'.

           nframe : (int)
              the frame number to read first [1 is the first].

           server : (bool)
              True/False for server vs local disk access
        """

        # read the header
        Rhead.__init__(self, fname, server)
        self.nframe = nframe
        self.server = server

    # Want to run this as a context manager
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.__del__()

    # and as an iterator.
    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self.__call__()
        except (HendError, urllib.error.HTTPError):
            raise StopIteration

    def ntotal(self):
        """
        Returns the total number of complete frames
        """
        if self.server:
            raise NotImplementedError('needs HiPERCAM server to be implemented')
        else:
            ntot = self.header['NAXIS3']
        return ntot

    def time(self, nframe=None):
        """Returns timing information of frame nframe (starts from 1). This saves
        effort reading the data in some cases. Once done, it moves to the next
        frame.

        Arguments::

           nframe : (int | None)
              frame number to get, starting at 1. 0 for the last (complete)
              frame. If nframe == None, the current frame will be used.

        See utimer for what gets returned by this. See also Rtime for a class
        dedicated to reading times only.
        """

        if self.server:
            raise NotImplementedError('needs HiPERCAM server to be implemented')
        else:
            raise NotImplementedError('needs HiPERCAM timing info implementing')

    def __call__(self, nframe=None):
        """Reads one exposure from the run the :class:`Rdata` is attached
        to. It works on the assumption that the internal file pointer in the
        :class:`Rdata` is positioned at the start of a frame. If `nframe` is
        None, then it will read the frame it is positioned at. If nframe is an
        integer > 0, it will try to read that particular frame; if nframe ==
        0, it reads the last complete frame.  nframe == 1 gives the first
        frame. This returns an MCCD object. It raises an exception if it fails
        to read data.  The data are stored internally as either 4-byte floats
        or 2-byte unsigned ints.

        Arguments::

           nframe : (int)
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted.

        Returns an MCCD for ULTRACAM, a CCD for ULTRASPEC.

        Apart from reading the raw bytes, the main job of this routine is to
        divide up and re-package the bytes read into Windats suitable for
        constructing CCD objects.

        """

        if self.server:
            raise NotImplementedError('needs HiPERCAM server access to be implemented')
        else:
            # timing bytes will need to be implemented

            # read data
            if nframe == 0:
                self.nframe = self.ntotal()
            elif nframe is not None:
                self.nframe = nframe

            if self.nframe > self.ntotal():
                self.nframe = 1
                raise HendError('Rdata.__call__: tried to access a frame exceeding the maximum')

            # read in frame
            frame = self.ffile[0][self.nframe-1,:,:]

            # Now build up Windats-->CCDs-->MCCD
            cnams = ('1', '2', '3', '4', '5')
            wnams1 = ('E1', 'F1', 'G1', 'H1')
            wnams2 = ('E2', 'F2', 'G2', 'H2')

            # Below, the intermediate arrays data1 and data2 are 4D arrays
            # indexed by (ccd,window,ny,nx) containing the image data. The crucial
            # step is contained in the as_strided routine

            ccds = Group()

            if self.mode.startswith('TwoWindows'):
                # get an example of the first set of Windows
                win1 = self.windows[0]

                # get view covering first 4 windows of data
                data_size1 = 20*win1.nx*win1.ny
                frame1 = frame[:,:,:data_size1]

                # re-view the data as a 4D array indexed by (ccd,window,y,x)
                data1 = as_strided(
                    frame1, strides=(8, 2, 40*win1.nx, 40),
                    shape=(5, 4, win1.ny, win1.nx)
                )

                # get an example of the second set of Windows
                win2 = self.windows[4]

                # get view covering second 4 windows of data
                # skip 32 timing bytes / 16 words at the end]
                frame2 = frame[:,:,data_size1:-16]

                # re-view as 4D array indexed by (ccd,window,y,x)
                data2 = as_strided(
                    frame2, strides=(8, 2, 40*win2.nx, 40),
                    shape=(5, 4, win2.ny, win2.nx)
                )

                # load into CCDs
                for nccd, (cnam, chead) in enumerate(zip(cnams,self.cheads)):
                    # now the Windats
                    windats = Group()
                    for nwin, (wnam1, wnam2) in enumerate(zip(wnams1, wnams2)):
                        windats[wnam1] = Windat(self.windows[nwin], data1[nccd, nwin])
                        windats[wnam2] = Windat(self.windows[nwin+4], data2[nccd, nwin])

                    # compile the CCD
                    ccds[cnam] = CCD(windats, HCM_NXTOT, HCM_NYTOT, chead)
            else:
                # get an example of the first set of Windows
                win1 = self.windows[0]

                # get view covering first 4 windows of data
                data_size1 = 20*win1.nx*win1.ny
                frame1 = frame[:,:,:data_size1]

                # re-view as 4D array indexed by (ccd,window,y,x)
                data1 = as_strided(
                    frame1, strides=(8, 2, 40*win1.nx, 40),
                    shape=(5, 4, win1.ny, win1.nx)
                )

                # load into CCDs
                for nccd, (cnam, chead) in enumerate(zip(cnams,self.cheads)):
                    # now the Windats
                    windats = Group()
                    for nwin, wnam1 in enumerate(wnams1):
                        windats[wnam1] = Windat(self.windows[nwin], data1[nccd, nwin])

                    # compile the CCD
                    ccds[cnam] = CCD(windats, HCM_NXTOT, HCM_NYTOT, chead)

            # finally create the MCCD
            mccd = MCCD(ccds, self.thead)

        # update the frame counter
        self.nframe += 1

        return mccd
