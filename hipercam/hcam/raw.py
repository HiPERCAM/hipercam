"""
Section for handling the raw HiPERCAM data files. These are FITS files with a
single HDU containing header information on the window format amonst other
things with the data in a "FITS cube" of unusual nature in that it has
dimensions of nframe x 1 x lots, where "lots" contains all pixels of all
windows of all CCDs plus timing bytes at the end. The pixels come in an odd
order (pixel 1 of each window and then the same for each CCD etc), and the main
work of the code is in disentangling them to produce :class:MCCD objects.
"""

import os
import struct
import warnings
import json
import websocket

import numpy as np
from numpy.lib.stride_tricks import as_strided
from astropy.io import fits
from astropy.time import Time

from hipercam import (CCD, Group, MCCD, Windat, Window, HCM_NXTOT, HCM_NYTOT)
from hipercam import (HRAW, add_extension) 
from .herrors import *

__all__ = ['Rhead', 'Rdata', 'Rtime']

# Set the URL of the FileServer from the environment or default to localhost.
URL = os.environ['HIPERCAM_DEFAULT_URL'] if 'HIPERCAM_DEFAULT_URL' in os.environ else \
      'ws://localhost:8007/'

# FITS offset to convert from signed to unsigned 16-bit ints and vice versa
BZERO = (1 << 15)

class Rhead:
    """Reads an interprets header information from a HiPERCAM run (3D FITS file)
    and generates the window formats needed to interpret the data. The file is
    kept open while the Rhead exists. The file can be on the local disk or obtained
    via the fileserver.

    The following attributes are set::

      cheads    : (list of astropy.io.fits.Header)
         the headers for each CCD which contain any info needed per CCD, above
         all the 'skip' numbers, which are used when creating MCCD frames.

      fname     : (string)
         the file name

      header    : (astropy.io.fits.Header)
         the header of the FITS file.

      mode      : (string)
         the readout mode used, read from 'ESO DET READ CURNAME'

      ntbytes   : (int)
         number of timing bytes, deduced from the dimensions of the FITS cube
         (NAXIS1 contains all pixels and timing bytes) compared to the window
         dimensions.

      thead     : (astropy.io.fits.Header)
         the top-level header that will be used to create MCCDs.

      windows   : (list of Windows)
         the windows are the same for each CCD, so these are just one set
         in EFGH order for window 1 and then the same for window 2 if there
         is one, so either 4 or 8 Windows.

    """

    def __init__(self, fname, server=False, full=True):
        """Opens a fits file containing HiPERCAM data, reads and interprets
        the header of the primary HDU.

        Arguments::

           fname  : (string)
              file name. '.fits' will be added to it. NB If server==True, this
              must have the form 'run010' (possibly with a path) as the server
              uses this to distinguish runs from other stuff.

           server : (bool)
              True/False for server vs local disk access. Server access goes
              through a websocket.  It uses a base URL taken from the
              environment variable "HIPERCAM_DEFAULT_URL", or, if that is not
              set, "ws://localhost:8007/".  The server here is Stu Littlefair's
              Python-based server that defaults to port 8007.

           server : (bool)
              True to access the data from a server. It uses a URL set in an
              environment variable "HIPERCAM_DEFAULT_URL" in this instance.

           full   : (bool)
              Flag controlling the amount of header detail to gather (as a
              time saver) True for detail, False for the bare minimum, the
              latter might be useful for large numbers of very small format
              images where you don't to waste effort.
        """

        # store the file name and whether server being used
        self.fname = fname
        self.server = server

        if server:
            # open socket connection to server
            self._ws = websocket.create_connection(URL + fname)

            # read the header from the server
            data = json.dumps(dict(action='get_hdr'))
            self._ws.send(data)
            self.header = fits.Header.fromstring(self._ws.recv())

        else:
            # open the file
            self._ffile = open(add_extension(fname, HRAW),'rb')

            # read the header
            self.header = fits.Header.fromfile(self._ffile)

            # store number of bytes read from header. used to offset later
            # requests for data
            self._hbytes = self._ffile.tell()

        # calculate the framesize in bytes
        bitpix = abs(self.header['BITPIX'])
        self._framesize = (bitpix*self.header['NAXIS1']) // 8

        # store the scaling and offset
        self._bscale = self.header['BSCALE']
        self._bzero = self.header['BZERO']

        # create the top-level header
        self.thead = fits.Header()

        # set the mode, one or two windows per quadrant, drift etc.
        # This is essential.
        self.mode = self.header['ESO DET READ CURNAME']
        self.thead['MODE'] = (self.mode,'HiPERCAM readout mode')

        # now we start on decoding the header parameters into Windows
        # required to construct MCCD objects

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

        # check for overscan / prescan
        oscan = self.header['ESO DET INCOVSCY']
        pscan = self.header['ESO DET INCPRSCX']
        yframe = 520 if oscan else 512
        xframe = 1074 if pscan else 1024

        # extract the data to build the first 4 windows
        if self.mode.startswith('FullFrame'):
            nx = xframe // xbin
            ny = yframe // ybin
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

        # track total number of pixels in order to compute the number of timing bytes
        npixels = 20*nx*ny

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

            # add in the pixels
            npixels += 20*nx*ny

        # any extra bytes above those needed for the data are timing
        self.ntbytes = self._framesize - (bitpix*npixels) // 8

        # Build (more) header info
        if 'DATE' in self.header:
            self.thead['DATE'] = (self.header['DATE'], self.header.comments['DATE'])
        if full and 'ESO DET GPS' in self.header:
            self.thead['GPS'] = (self.header['ESO DET GPS'], self.header.comments['ESO DET GPS'])

        # Header per CCD
        self.cheads = []
        for n in range(1,6):
            chead = fits.Header()

            # Essential items
            hnam = 'ESO DET NSKIP{:d}'.format(n)
            if hnam in self.header:
                chead['NSKIP'] = (self.header[hnam], self.header.comments[hnam])

            # Nice-if-you-can-get-them items
            if full:
                # whether this CCD has gone through a reflection
                hnam = 'ESO DET REFLECT{:d}'.format(n)
                if hnam in self.header:
                    chead['REFLECT'] = (self.header[hnam], 'is image reflected')

                # readout noise
                hnam = 'ESO DET CHIP{:d} RON'.format(n)
                if hnam in self.header:
                    chead['RONOISE'] = (self.header[hnam], self.header.comments[hnam])

                # gain
                hnam = 'ESO DET CHIP{:d} GAIN'.format(n)
                if hnam in self.header:
                    chead['GAIN'] = (self.header[hnam], self.header.comments[hnam])

            self.cheads.append(chead)

        # end of constructor / initialiser

    def __del__(self):
        """Destructor closes the file or web socket"""
        if self.server:
            if hasattr(self, '_ws'): self._ws.close()
        else:
            if hasattr(self, '_ffile'): self._ffile.close()

    def npix(self):
        """
        Returns the number of (binned) pixels per CCD
        """
        np = 0
        for win in self.windows:
            np += win.nx*win.ny
        return np

    def ntotal(self):
        """
        Returns the total number of complete frames
        """
        if self.server:
            raise NotImplementedError('needs HiPERCAM server to be implemented')
        else:
            ntot = self.header['NAXIS3']
        return ntot

    # Want to run this as a context manager
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.__del__()


class Rdata (Rhead):
    """Callable, iterable object to represent HiPERCAM raw data files.

    The idea is to instantiate an Rdata from a HiPERCAM FITS file and then the
    object generated can be used to deliver MCCD objects by specifying a frame
    number e.g.::

      >> rdat = Rdata('run045.fits')
      >> fr10 = rdat(10)
      >> fr11 = rdat()

    reads frame number 10 and then 11 in from 'run045.fits', or sequentially::

      >> for mccd in Rdata('run045.fits'):
      >>    print('nccd = ',len(mccd))

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
              the frame number to read first [1 is the first]. This initialises an attribute
              of the same name that is used when reading frames sequentially.

           server : (bool)
              True/False for server vs local disk access. Server access goes
              through a websocket.  It uses a base URL taken from the
              environment variable "HIPERCAM_DEFAULT_URL", or, if that is not
              set, "ws://localhost:8007/".  The server here is Stu Littlefair's
              Python-based server that defaults to port 8007.
        """

        # read the header
        Rhead.__init__(self, fname, server)

        # move to the start of the frame
        self.nframe = nframe
        if not server:
            self._ffile.seek(self._hbytes+self._framesize*(self.nframe-1))

    # Rdata objects functions as iterators.
    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self.__call__()
        except (HendError):
            raise StopIteration

    def __call__(self, nframe=None):
        """Reads one exposure from the run the :class:`Rdata` is attached
        to. If `nframe` is None, then it will read the frame it is positioned
        at. If nframe is an integer > 0, it will try to read that particular
        frame; if nframe == 0, it reads the last complete frame.  nframe == 1
        gives the first frame. If all works, an MCCD object is returned. It
        raises an HendError if it reaches the end of a local disk file. If it
        fails to read from the server it returns None, but does not raise an
        Exception in order to allow continued attempts as reading from what
        might be a growing file.

        Arguments::

           nframe : (int)
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted.

        Apart from reading the raw bytes, the main job of this routine is to
        divide up and re-package the bytes read into Windats suitable for
        constructing CCD objects.

        If access via a server is requested, it is assumed that the file being
        accessed could be being added to. In
        """

        # set the frame to be read and whether the file pointer needs to be
        # reset
        if nframe is None:
            # just read whatever frame we are on
            reset = False
        else:
            if nframe == 0:
                # go for the last one
                nframe = self.ntotal()

            # update frame counter
            reset = self.nframe != nframe
            self.nframe = nframe

        # don't attempt to read off end of disk file
        if not self.server and self.nframe > self.ntotal():
            self.nframe = 1
            raise HendError('Rdata.__call__: tried to access a frame exceeding the maximum')

        if self.server:
            # define command to send to server
            if reset:
                # requires a seek as well as a read
                data = json.dumps(dict(action='get_frame', frame_number=self.nframe))
            else:
                # just read next set of bytes
                data = json.dumps(dict(action='get_next'))

            # get data
            self._ws.send(data)
            raw_bytes = self._ws.recv()

            if len(raw_bytes) == 0:
                # if we are trying to access a file that has not yet been written
                # 0 bytes will be returned. In this case we return with None. NB we
                # do not upfate self.nframe in this case.
                return None

            # separate into the frame and timing data, correcting the frame
            # bytes for the FITS BZERO offset
            frame  = np.fromstring(raw_bytes[:-self.ntbytes], '>u2')
            frame += BZERO
            tbytes = raw_bytes[-self.ntbytes:]

        else:
            # find the start of the frame if necessary
            if reset:
                self._ffile.seek(self._hbytes+self._framesize*(self.nframe-1))

            # read in frame and then the timing data, correcting the frame for
            # the standard FITS BZERO offset
            frame  = np.fromfile(self._ffile, '>u2', (self._framesize - self.ntbytes) // 2)
            frame += BZERO
            tbytes = self._ffile.read(self.ntbytes)
            if len(tbytes) != self.ntbytes:
                raise IOError('failed to read frame from disk file')

        # we now have all the pixels and timing bytes of the frame of interest
        # read into a 1D ndarray of unsigned 2-byte integers. Now interpret them.

        # first the time
        frameCount, timeStampCount, years, day_of_year, hours, mins, seconds, nanoseconds, nsats, synced = \
            decode_timing_bytes(tbytes)

        if frameCount+1 != self.nframe:
            # check that frameCount is as expected
            warnings.warn(
                'unexpected clash timestamp frame number = {:d} with that expected = {:d}'.format(
                    frameCount+1, self.nframe
                    ))

        time_string = '{}:{}:{}:{}:{:.7f}'.format(
            years, day_of_year, hours, mins, seconds+nanoseconds/1e9
            )
        tstamp = Time(time_string, format='yday')
        self.thead['TIMSTAMP'] = (tstamp.isot, 'Raw frame timestamp, UTC')
        self.thead['NFRAME'] = (frameCount+1, 'Frame number')

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
            frame1 = frame[:data_size1]

            # re-view the data as a 4D array indexed by (ccd,window,y,x)
            data1 = as_strided(
                frame1, strides=(8, 2, 40*win1.nx, 40),
                shape=(5, 4, win1.ny, win1.nx)
                )

            # get an example of the second set of Windows
            win2 = self.windows[4]

            # get view covering second 4 windows of data
            frame2 = frame[data_size1:]

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

            # re-view as 4D array indexed by (ccd,window,y,x)
            data1 = as_strided(
                frame, strides=(8, 2, 40*win1.nx, 40),
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

class Rtime (Rhead):
    """Callable, iterable object to represent the timestamps only of HiPERCAM raw data files.

    This is very nearly the same as Rdata except for speed it only reads the
    timing data of each frame in the case of local disk files and it does no
    processing of raw data other than the timing in any circumstance.  Each
    call to it returns an astropy.time.Time object, thus

       >> for time in Rtime('run045'):
       >>     print(time.isot)

    lists all the timestamps of the run.
    """

    def __init__(self, fname, nframe=1, server=False):
        """Connects to a raw HiPERCAM FITS file for reading. The file is kept
        open.  The Rdata object can then generate MCCD objects through being
        called as a function or iterator.

        Arguments::

           fname : (string)
              run name, e.g. 'run036'.

           nframe : (int)
              the frame number to read first [1 is the first]. This initialises an attribute
              of the same name that is used when reading frames sequentially.

           server : (bool)
              True/False for server vs local disk access. Server access goes
              through a websocket.  It uses a base URL taken from the
              environment variable "HIPERCAM_DEFAULT_URL", or, if that is not
              set, "ws://localhost:8007/".  The server here is Stu Littlefair's
              Python-based server that defaults to port 8007.
        """

        # read the header
        Rhead.__init__(self, fname, server)

        # move to the start of the timing bytes
        self.nframe = nframe
        if not server:
            self._ffile.seek(self._hbytes+self._framesize*self.nframe-self.ntbytes)

    # Rtime objects functions as iterators.
    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self.__call__()
        except (HendError):
            raise StopIteration

    def __call__(self, nframe=None):
        """Reads the timing data of one frame from the run the :class:`Rtime`
        is attached to. If `nframe` is None, then it will read the frame it is
        positioned at. If nframe is an integer > 0, it will try to read that
        particular frame; if nframe == 0, it reads the last complete frame.
        nframe == 1 gives the first frame. If all works, an astropy.time.Time
        object is returned. 

        Arguments::

           nframe : (int)
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted.

        """

        # set the frame to be read and whether the file pointer needs to be
        # reset
        if nframe is None:
            # just read whatever frane we are on
            reset = False
        else:
            if nframe == 0:
                # go for the last one
                nframe = self.ntotal()

            # update frame counter
            reset = self.nframe != nframe
            self.nframe = nframe

        # ?? need a maximum frame number routine in the server 
        # at the moment just don't run ntotal if server
        if not self.server and self.nframe > self.ntotal():
            self.nframe = 1
            raise HendError('Rdata.__call__: tried to access a frame exceeding the maximum')

        if self.server:
            # define command to send to server
            if reset:
                # requires a seek as well as a read
                data = json.dumps(dict(action='get_frame', frame_number=self.nframe))
            else:
                # just read next set of bytes
                data = json.dumps(dict(action='get_next'))

            # get data
            self._ws.send(data)
            raw_bytes = self._ws.recv()

            # extract the timing data
            tbytes = raw_bytes[-self.ntbytes:]

        else:
            # find the start of the timing data if necessary
            if reset:
                self._ffile.seek(self._hbytes+self._framesize*self.nframe-self.ntbytes)

            # read the timing data, skip forward to next set
            tbytes = self._ffile.read(self.ntbytes)
            self._ffile.seek(self._framesize-self.ntbytes, 1)

        # Now interpret them.
        frameCount, timeStampCount, years, day_of_year, hours, mins, seconds, nanoseconds, nsats, synced = \
            decode_timing_bytes(tbytes)

        time_string = '{}:{}:{}:{}:{:.7f}'.format(
            years, day_of_year, hours, mins, seconds+nanoseconds/1e9
            )
        tstamp = Time(time_string, format='yday')

        # update the frame counter
        self.nframe += 1

        # Return
        return tstamp

def decode_timing_bytes(tbytes):
    """Decode the timing bytes tacked onto the end of every HiPERCAM frame in the 3D FITS file.

    The timing bytes are encoded as 8x4-byte little-endian unsigned integers followed
    by 2x2-byte little-endian . On saving to FITS they are mangled
    because:

      1) they are split into 2-byte uints
      2) 2**15 is subtracted from each to map to 2-byte signed ints
      3) These are written in big-endian order to FITS

    This routine reverses this process to recover the original timestamp tuple by

      1) unpacking the timing bytes as big-endian 16-bit signed integers
      2) adding 32768 to each value
      3) packing these values as little-endian 16-bit unsigned integers
      4) unpacking the result as 32-bit, little-endian unsigned integers

    Parameters
    ----------

       tbytes: (bytes)
           a Python bytes object which contains the timestamp bytes as written in the
           FITS file

    Returns
    --------

      timestamp : (tuple)
          a tuple containing (frameCount, timeStampCount, years, day_of_year,
          hours, mins, seconds, nanoseconds, nsats, synced) values.

    [Straight from Stu Littlefair's original code.]

    """

    # format with 36 timing bytes (last two of which are ignored)
    buf = struct.pack('<HHHHHHHHHHHHHHHHH',
                      *(val + BZERO for val in struct.unpack('>hhhhhhhhhhhhhhhhh', tbytes[:-2])))
    return struct.unpack('<IIIIIIIIbb', buf)
