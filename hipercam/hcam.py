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

from hipercam import (CCD, Group, MCCD, Window, Winhead, HRAW, HipercamError)
from hipercam.utils import add_extension

__all__ = ['Rhead', 'Rdata', 'Rtime', 'HCM_NXTOT', 'HCM_NYTOT']

# Set the URL of the FileServer from the environment or default to localhost.
URL = 'ws://{:s}'.format(os.environ['HIPERCAM_DEFAULT_URL']) \
      if 'HIPERCAM_DEFAULT_URL' in os.environ else 'ws://localhost:8007/'

# FITS offset to convert from signed to unsigned 16-bit ints and vice versa
BZERO = (1 << 15)

# Maximum dimensions of HiPERCAM CCDs
HCM_NXTOT, HCM_NYTOT = 2148, 1040

# number of seconds in a day
DAYSEC = 86400.

class Rhead:
    """Reads an interprets header information from a HiPERCAM run (3D FITS file)
    and generates the window formats needed to interpret the data. The file is
    kept open while the Rhead exists. The file can be on the local disk or
    obtained via the fileserver.

    The following attributes are set::

       cheads    : list of astropy.io.fits.Header
           the headers for each CCD which contain any info needed per CCD, above
           all the 'skip' numbers, which are used when creating MCCD frames.

       clear, dummy, oscan, prescan  : bool
           flags to say whether or not the following were enabled:
               clears between exposures
               dummy output
               overscan (extra pixels in Y)
               prescan (extra pixels in X)

       fname     : string
           the file name

       header    : astropy.io.fits.Header
           the header of the FITS file.

       mode      : string
           the readout mode used, read from 'ESO DET READ CURNAME'

       nskips    : tuple
           the NSKIP parameters for each CCD (5 of them)

       ntbytes   : int
           number of timing bytes, deduced from the dimensions of the FITS cube
           (NAXIS1 contains all pixels and timing bytes) compared to the window
           dimensions.

       nwins     : tuple
           integer indices of the windows for any one quadrant. This is (0,) for
           one window modes, (0,1) for two window modes.

       tdead, tdelta, toff1, toff2, toff3, toff4 : float
           set of parameters used to calculate exact times. tdead is the
           dead time between exposures; tdelta is used to multiply NSKIP
           for each CCD; toff1 etc are offsets and lengths used to get the
           mid-exposure times and exposure durations.

       thead     : astropy.io.fits.Header
           the top-level header that will be used to create MCCDs.

       windows   : list of list of list of (Winhead, flip) tuples
           windows[nccd][nquad][nwin] returns the Winhead and set of axes for
           window nwin of quadrant nquad of CCD nccd. nwin indexes the windows
           of a given quadrant (at most there are 2, so nwin = 0 or 1). nquad
           indexes the quadrants in EFGH order. nccd index the 5 CCDs. 'flip'
           is a tuple of axes to be flipped using numpy.flip to account for
           where the CCD is read and reflections.

       wforms    : tuple of strings
           formats of each set of windows as strings designed to match the input
           expectd for hdriver.

       xbin      : int
           X-binning factor

       ybin      : int
           Y-binning factor

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
            status = json.loads(self._ws.recv())['status']
            if status == 'no such run':
                raise HipercamError('Run not found: {}'.format(fname))

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

        # set the mode, one or two windows per quadrant, drift etc.  This is
        # essential.
        self.mode = self.header['ESO DET READ CURNAME']
        self.thead['MODE'] = (self.mode, 'HiPERCAM readout mode')

        # check for overscan / prescan / clear
        self.clear = self.header['ESO DET CLRCCD']
        self.dummy = self.header['ESO DET DUMMY']
        self.oscan = self.header['ESO DET INCOVSCY']
        self.pscan = self.header['ESO DET INCPRSCX']

        # framesize parameters
        yframe = 520 if self.oscan else 512
        xframe = 1074 if self.pscan else 1024

        # store the NSKIP values for each CCD. These are needed for exact
        # timing for each CCD.
        self.nskips = (
            self.header['ESO DET NSKIPS1'],
            self.header['ESO DET NSKIPS2'],
            self.header['ESO DET NSKIPS3'],
            self.header['ESO DET NSKIPS4'],
            self.header['ESO DET NSKIPS5']
            )

        # set timing parameters toff1, toff2, toff3, toff4, tdelta, tdead
        # which are also needed for exact timing
        if self.clear:
            # Clear mode. NB the parameter 'TREAD' is actually a combination
            # of a frame-transfer and a read
            self.tdelta = 1e-3*(
                self.header['ESO DET TREAD'] + self.header['ESO DET TDELAY'] + \
                self.header['ESO DET TCLEAR']
            )

            self.toff1 = self.toff2 = 1.e-3*self.header['ESO DET TDELAY']/2.

            self.toff3 = self.toff4 = 1.e-3*self.header['ESO DET TDELAY']

            self.tdead = 1.e-3*(self.header['ESO DET TREAD']+self.header['ESO DET TCLEAR'])

        else:
            # No clear mode, there are zero sec dummy 'wipes' so no 'TCLEAR'
            # here
            self.tdelta = 1.e-3*(
                self.header['ESO DET TREAD'] + self.header['ESO DET TDELAY']
            )

            self.toff1 = 1.e-3*self.header['ESO DET TDELAY']/2.

            self.toff2 = 1.e-3*(
                self.header['ESO DET TDELAY'] - self.header['ESO DET TREAD'] +
                self.header['ESO DET TFT'] ) / 2.

            self.toff3 = 1.e-3*self.header['ESO DET TDELAY']

            self.toff4 = 1.e-3*(
                self.header['ESO DET TDELAY'] + self.header['ESO DET TREAD'] - \
                self.header['ESO DET TFT']
            )

            self.tdead = 1.e-3*self.header['ESO DET TFT']

        # binning factors
        self.xbin = self.header['ESO DET BINX1']
        self.ybin = self.header['ESO DET BINY1']

        # This mapping is needed for the g and z channels due to reflections.
        # They need to be flipped around x=0 and be rotated 180 degrees so that
        # whatever falls on E in the other CCDs falls on F
        QUAD_MAPPING = {'E': 'F', 'F': 'E', 'G': 'H', 'H': 'G'}
        QUAD = ('E', 'F', 'G', 'H')

        # bottom-left coordinate of quadrant
        LLX = {'E': 1, 'F': xframe+1, 'G': xframe+1, 'H': 1}
        LLY = {'E': 1, 'F': 1, 'G': yframe+1, 'H': yframe+1}

        # direction increasing X- or Y-start moves window in quad
        X_DIRN = {'E': 1, 'F': -1, 'G': -1, 'H': 1}
        Y_DIRN = {'E': 1, 'F': 1, 'G': -1, 'H': -1}

        # whether we need to add size of window to find bottom-left
        # coordinate of window
        ADD_YSIZES = {'E': 0, 'F': 0, 'G': 1, 'H': 1}
        ADD_XSIZES = {'E': 0, 'F': 1, 'G': 1, 'H': 0}

        # whether to offset llx for prescan pixels (basically
        # inverse of ADD_XSIZES)
        OFF_PRESCAN = {'E': 1, 'F': 0, 'G': 0, 'H': 1}

        # some axes need flipping, since the data is read out such that the
        # bottom-left is not always the first pixel which axes need flipping
        # to get data in correct orientation using numpy C-convention, 1=y,
        # 0=x. These are fed to numpy.flip later
        FLIP_AXES = {'E': (), 'F': (1,), 'G': (1, 0), 'H': (0,)}

        if self.mode.startswith('TwoWindow'):
            self.nwins = (0,1)
        else:
            self.nwins = (0,)

        # build windows for each window of each quadrant of each CCDs
        self.windows = []
        for nwin in self.nwins:
            self.windows.append([])

            # first set the windows of each CCD
            for nccd in range(5):
                self.windows[-1].append([])
                for quad in QUAD:
                    if nccd in (1, 4):
                        # reflections for g (1) and z (4)
                        quad = QUAD_MAPPING[quad]

                    if self.mode.startswith('FullFrame'):
                        nx = xframe // self.xbin
                        ny = yframe // self.ybin

                        if self.pscan:
                            # pre-scan adds 50 columns to every window
                            # unbinned with some slightly tricky details
                            # when binned
                            ipoff = int(np.ceil(50 / self.xbin))
                        else:
                            ipoff = 0

                        llx = LLX[quad] - OFF_PRESCAN[quad]*ipoff
                        lly = LLY[quad]

                    elif self.mode.startswith('OneWindow') or \
                         self.mode.startswith('TwoWindow'):
                        winID = 'ESO DET WIN{} '.format(nwin+1)
                        nx = self.header[winID + 'NX'] // self.xbin

                        if self.pscan:
                            # pre-scan adds 50 columns to every window
                            # unbinned with some slightly tricky details
                            # when binned
                            ipoff = int(np.ceil(50 / self.xbin))
                            nx += ipoff
                        else:
                            ipoff = 0

                        ny = self.header[winID + 'NY'] // self.ybin
                        win_xs = self.header[winID + 'XS{}'.format(quad)]
                        win_nx = self.header[winID + 'NX']
                        win_ys = self.header[winID + 'YS']
                        win_ny = self.header[winID + 'NY']
                        llx = (
                            LLX[quad] + X_DIRN[quad] * win_xs +
                            ADD_XSIZES[quad] * (xframe - win_nx) -
                            OFF_PRESCAN[quad]*ipoff
                        )
                        lly = (
                            LLY[quad] + Y_DIRN[quad] * win_ys +
                            ADD_YSIZES[quad] * (yframe - win_ny)
                        )

                    elif self.mode.startswith('Drift'):
                        nx = self.header['ESO DET DRWIN NX'] // self.xbin
                        ny = self.header['ESO DET DRWIN NY'] // self.ybin

                        if self.pscan:
                            # pre-scan adds 50 columns to every window
                            # unbinned with some slightly tricky details
                            # when binned
                            ipoff = int(np.ceil(50 / self.xbin))
                            nx += ipoff
                        else:
                            ipoff = 0

                        win_xs = self.header['ESO DET DRWIN XS{}'.format(quad)]
                        win_nx = self.header['ESO DET DRWIN NX']
                        win_ys = self.header['ESO DET DRWIN YS']
                        win_ny = self.header['ESO DET DRWIN NY']
                        llx = (
                            LLX[quad] + X_DIRN[quad] * win_xs +
                            ADD_XSIZES[quad] * (xframe - win_nx) -
                            OFF_PRESCAN[quad]*ipoff
                        )
                        lly = (
                            LLY[quad] + Y_DIRN[quad] * win_ys +
                            ADD_YSIZES[quad] * (yframe - win_ny)
                        )

                    else:
                        msg = 'mode {} not currently supported'.format(mode)
                        raise ValueError(msg)

                    # store the window and the axes to flip
                    self.windows[-1][-1].append(
                        (Winhead(llx, lly, nx, ny, self.xbin, self.ybin),
                         FLIP_AXES[quad])
                    )


        # set strings for logging purposes giving the settings used in hdriver
        self.wforms = []

        if self.mode.startswith('FullFrame'):
            self.wforms.append('Full&nbsp;Frame')

        elif self.mode.startswith('OneWindow') or \
             self.mode.startswith('TwoWindow'):
            winID = 'ESO DET WIN1 '
            ys = self.header[winID + 'YS'] + 1
            nx = self.header[winID + 'NX']
            ny = self.header[winID + 'NY']

            try:
                xsll = self.header[winID + 'XSLL']
                xslr = self.header[winID + 'XSLR']
                xsul = self.header[winID + 'XSUL']
                xsur = self.header[winID + 'XSUR']
            except:
                # Fallback
                xsll = self.header[winID + 'XSE'] + 1
                xsul = self.header[winID + 'XSH'] + 1
                xslr = 2049 - self.header[winID + 'XSF'] - nx
                xsur = 2049 - self.header[winID + 'XSG'] - nx

            self.wforms.append(
                '{:d}|{:d}|{:d}|{:d}|{:d}|{:d}|{:d}'.format(
                    xsll, xslr, xsul, xsur, ys, nx, ny)
                )

            if self.mode.startswith('TwoWindow'):
                winID = 'ESO DET WIN2 '
                ys = self.header[winID + 'YS'] + 1
                nx = self.header[winID + 'NX']
                ny = self.header[winID + 'NY']

                try:
                    xsll = self.header[winID + 'XSLL']
                    xslr = self.header[winID + 'XSLR']
                    xsul = self.header[winID + 'XSUL']
                    xsur = self.header[winID + 'XSUR']
                except:
                    # fallback option
                    xsll = self.header[winID + 'XSE'] + 1
                    xsul = self.header[winID + 'XSH'] + 1
                    xslr = 2049 - self.header[winID + 'XSF'] - nx
                    xsur = 2049 - self.header[winID + 'XSG'] - nx

                self.wforms.append(
                    '{:d}|{:d}|{:d}|{:d}|{:d}|{:d}|{:d}'.format(
                        xsll, xslr, xsul, xsur, ys, nx, ny)
                    )

        elif self.mode.startswith('Drift'):
            winID = 'ESO DET DRWIN '
            ys = self.header[winID + 'YS'] + 1
            nx = self.header[winID + 'NX']
            ny = self.header[winID + 'NY']

            try:
                xsl = self.header[winID + 'XSL']
                xsr = self.header[winID + 'XSR']
            except:
                # Fallback option
                xsl = self.header[winID + 'XSE'] + 1
                xsr = 2049 - self.header[winID + 'XSF'] - nx

            self.wforms.append(
                '{:d}|{:d}|{:d}|{:d}|{:d}'.format(xsl, xsr, ys, nx, ny)
                )

        # compute the total number of pixels (making use of symmetry between
        # quadrants and CCDs, hence factor 20)
        npixels = 0
        for nwin in self.nwins:
            win = self.windows[nwin][0][0][0]
            npixels += 20*win.nx*win.ny

        # any extra bytes above those needed for the data are timing
        self.ntbytes = self._framesize - (bitpix*npixels) // 8

        # Build (more) header info
        if 'DATE' in self.header:
            self.thead['DATE'] = (
                self.header['DATE'], self.header.comments['DATE']
            )

        if full and 'ESO DET GPS' in self.header:
            self.thead['GPS'] = (self.header['ESO DET GPS'],
                                 self.header.comments['ESO DET GPS'])

        if full:
            if 'EXPTIME' in self.header:
                self.thead['EXPTIME'] = (
                    self.header['EXPTIME'], self.header.comments['EXPTIME']
                )
            self.thead['XBIN'] = (self.xbin,
                                  self.header.comments['ESO DET BINX1'])
            self.thead['YBIN'] = (self.ybin,
                                  self.header.comments['ESO DET BINY1'])
            self.thead['SPEED'] = (self.header['ESO DET SPEED'],
                                   self.header.comments['ESO DET SPEED'])

        # Header per CCD. These are modified per CCD in Rdata
        self.cheads = []
        for n in range(5):
            chead = fits.Header()

            # Essential items
            hnam = 'ESO DET NSKIPS{:d}'.format(n+1)
            if hnam in self.header:
                # used to identify blank frames
                chead['NCYCLE'] = (
                    self.header[hnam]+1, 'readout cycle period (NSKIP+1)'
                )

            # Nice-if-you-can-get-them items
            if full:

                # whether this CCD has gone through a reflection
                hnam = 'ESO DET REFLECT{:d}'.format(n+1)
                if hnam in self.header:
                    chead['REFLECT'] = (self.header[hnam], 'is image reflected')

                # readout noise
                hnam = 'ESO DET CHIP{:d} RON'.format(n+1)
                if hnam in self.header:
                    chead['RONOISE'] = (self.header[hnam],
                                        self.header.comments[hnam])

                # gain
                hnam = 'ESO DET CHIP{:d} GAIN'.format(n+1)
                if hnam in self.header:
                    chead['GAIN'] = (self.header[hnam],
                                     self.header.comments[hnam])

            self.cheads.append(chead)

        # end of constructor / initialiser

    def ntotal(self):
        """
        Returns the total number of complete frames.
        """
        if self.server:
            data = json.dumps(dict(action='get_nframes'))
            self._ws.send(data)
            raw_bytes = self._ws.recv()
            if len(raw_bytes) == 0:
                raise hcam.HipercamError(
                    'failed to get total number of frames'
                    ' from server; 0 bytes returned'
                )
            d = eval(raw_bytes)
            ntot = d['nframes']
        else:
            """
            # rather than read NAXIS3, this is designed to allow for the
            # file being incremented. First record where we are so we can
            # return
            ctell = self._ffile.tell()

            # seek the end of file
            self._ffile.seek(0,2)

            # compute byte offset at this point and thus the total number of
            # complete frames (this is what we return)
            ntot = (self._ffile.tell() - self._hbytes) // self._framesize

            print( ctell, self._ffile.tell(), self._hbytes, self._framesize,
                   (self._ffile.tell() - self._hbytes) / self._framesize)

            # return pointer to starting position
            self._ffile.seek(ctell)

            """
            # the commented out section was an attempt to work out the number
            # of frames even when the file was growing. This fails owing the
            # 2880-byte FITS standard when one has small framesizes, thus I am
            # using NAXIS3 instead
            ntot = self.header['NAXIS3']

        return ntot

    def tinfo(self):
        """Returns with timing info for the run used to create the Rhead, namely the
        exposure times for each CCD appropriate for a "standard" exposure
        (i.e. not the first), the offsets to get to the mid-exposure from a
        timestamp, again for a "standard" exposure, the cycle length for each
        CCD (NSUB+1) and the dead time between frames. This is intended for use
        in generating logs.

        Returns (texps, toffs, ncycs, tdead) where texps, toffs and ncycs are
        5-element tuples representing the "standard" exposure times, the
        offsets from the timestamp, and the integer cycle lengths for each,
        and tdead is the deadtime. The times are all in seconds. For standard
        exposure in which a CCD is read out every time it can be, NCYC=1.

        """

        # we just pretend we are on a suitable frame for each CCD
        texps = []
        toffs = []
        ncycs = []
        for nccd, nskip in enumerate(self.nskips):
            toff, texp, flag = self.timing(0, 2*(nskip+1), nccd)
            toffs.append(DAYSEC*toff)
            texps.append(texp)
            ncycs.append(nskip+1)

        return (tuple(texps), tuple(toffs), tuple(ncycs), self.tdead)

    def timing(self, mjd, nframe, nccd):
        """Returns timing data for a particular CCD and frame.

        Arguments::

           mjd    : float
              the MJD of the timestamp

           nframe : int
              the frame number starting at 1.

           nccd : int
              the CCD number, starting at 0

        Returns: (mjd, exptime, flag) where 'mjd' is the MJD at mid-exposure for the CCD, 'exptime'
        is the exposure time in seconds, and 'flag' is a bool to indicate whether the time is thought
        good or not.
        """

        nskip = self.nskips[nccd]

        flag = nframe % (nskip+1) == 0
        if nframe == nskip + 1:
            tmid = mjd + (self.toff1 - self.tdelta*nskip/2)/DAYSEC
            texp = self.tdelta*nskip + self.toff3
        else:
            tmid = mjd + (self.toff2 - self.tdelta*nskip/2)/DAYSEC
            texp = self.tdelta*nskip + self.toff4
        return (tmid, texp, flag)

    def __del__(self):
        """Destructor closes the file or web socket"""
        if self.server:
            if hasattr(self, '_ws'): self._ws.close()
        else:
            if hasattr(self, '_ffile'): self._ffile.close()

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
    (=5 in all cases).

    Note: this casts all incoming data into 32-bit floats. This is perhaps
    inefficient memory-wise but it is safer than leaving a 2-bytes ints which
    can bite you.

    """

    def __init__(self, fname, nframe=1, server=False):
        """Connects to a raw HiPERCAM FITS file for reading. The file is kept
        open.  The Rdata object can then generate MCCD objects through being
        called as a function or iterator.

        Arguments::

           fname : (string)
              run name, e.g. 'run036'.

           nframe : (int)
              the frame number to read first [1 is the first]. This
              initialises an attribute of the same name that is used when
              reading frames sequentially. nframe=0 is an indication to set an
              attribute 'last' = True to indicate that it should always try to
              access the last frame.

           server : (bool)
              True/False for server vs local disk access. Server access goes
              through a websocket.  It uses a base URL taken from the
              environment variable "HIPERCAM_DEFAULT_URL", or, if that is not
              set, "ws://localhost:8007/".  The server here is Stu Littlefair's
              Python-based server that defaults to port 8007.

        """

        # read the header
        Rhead.__init__(self, fname, server)

        # flag to indicate should always try to get the last frame
        self.last = (nframe == 0)

        # set flag to indicate first time through
        self.first = True

        if server:
            self.nframe = nframe
        else:
            # local file: go to correct location
            if self.last:
                self.seek_last()
            else:
                self.seek_frame(nframe)

    # Rdata objects functions as iterators.
    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self.__call__()
        except (HendError):
            raise StopIteration

    def __call__(self, nframe=None):
        """Reads one exposure from the run the :class:`Rdata` is attached to. If
        `nframe` is None, then it will read the frame it is positioned at. If
        nframe is an integer > 0, it will try to read that particular frame;
        if nframe == 0, it reads the last complete frame. nframe == 1 gives
        the first frame. If all works, an MCCD object is returned. It raises
        an HendError if it reaches the end of a local disk file. If it fails
        to read from the server it returns None, but does not raise an
        Exception in order to allow continued attempts as reading from what
        might be a growing file.

        It maintains an attribute 'nframe' corresponding to the frame that the
        file pointer is positioned to read next (most relevant to reading of a
        local file). If it is set to read the last file and that file doesn't
        change it returns None as a sign of failure.

        Arguments::

           nframe : (int)
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted, unless
              self.nframe = 0 in which case it will try to get the most recent
              frame, whatever that is.

        Apart from reading the raw bytes, the main job of this routine is to
        divide up and re-package the bytes read into Windows suitable for
        constructing CCD objects.

        If access via a server is requested, it is assumed that the file being
        accessed could be being added to. In
        """

        if self.server:

            # access frames via the server

            # build the request
            if nframe == 0 or self.last:
                # just want the last complete frame. Note that further
                # down we check the frame counter against what we were
                # hoping for (at least 1 further on) and if it isn't
                # we actually return with None to indicate no progress.
                request = json.dumps(dict(action='get_last'))

            elif nframe is None:
                # in this case, on the first time through we need to
                # explicitly request the frame, otherwise we can just read
                # from where we are
                if self.first:
                    request = json.dumps(
                        dict(action='get_frame', frame_number=self.nframe))
                else:
                    request = json.dumps(dict(action='get_next'))

            else:
                # a particular frame number is being requested. Check whether
                # we have to request it explicitly or whether we can just get
                # the next one
                if self.nframe != nframe:
                    request = json.dumps(
                        dict(action='get_frame', frame_number=nframe))
                    self.nframe = nframe
                else:
                    request = json.dumps(dict(action='get_next'))

            # send the request
            self._ws.send(request)
            raw_bytes = self._ws.recv()

            if len(raw_bytes) == 0:
                # if we are trying to access a file that has not yet been
                # written, 0 bytes will be returned. In this case we return
                # with None. NB we do not update self.nframe in this case.
                return None

            # separate into the frame and timing data, correcting the frame
            # bytes for the FITS BZERO offset
            frame  = np.fromstring(raw_bytes[:-self.ntbytes], '>u2')
            frame += BZERO
            tbytes = raw_bytes[-self.ntbytes:]

        else:

            # access frames in local file

            if nframe == 0 or self.last:

                # get the number of frames
                nold = self.nframe
                ntot = self.seek_last()

                if ntot == 0:
                    # nothing ready yet
                    return None

                if self.first:
                    # first time through
                    self.nframe = ntot

                elif ntot < nold:
                    # indicates no new frame has been added since
                    # previous call after which the frame pointer
                    # was incremented to the value stored in nold
                    self.nframe += 1
                    return None

            elif nframe is not None:
                # move to correct place
                self.seek_frame(nframe)

            # read in frame and then the timing data, correcting the frame for
            # the standard FITS BZERO offset. At this stage we have the data
            # as unsigned 2-byte ints
            frame = np.fromfile(
                self._ffile, '>u2', (self._framesize - self.ntbytes) // 2)
            frame += BZERO
            tbytes = self._ffile.read(self.ntbytes)
            if len(tbytes) != self.ntbytes:
                raise HendError('failed to read frame from disk file')

        ##############################################################
        #
        # We now have all the pixels (in 'frame') and timing bytes (in
        # 'tbytes') of the frame of interest read into a 1D ndarray of
        # unsigned, 2-byte integers. Now interpret them.
        #
        ##############################################################

        # First the timing bytes. The frameCount starts from 0 so we
        # we add one to it
        frameCount, timeStampCount, years, day_of_year, hours, mins, \
            seconds, nanoseconds, nsats, synced = decode_timing_bytes(tbytes)
        frameCount += 1

        if self.server and (nframe == 0 or self.last) and \
           not self.first and self.nframe > frameCount:
            # server access tring to get the last complete frame. If the frame
            # just read in is the same as the one before (i.e. frameCount <
            # self.nframe), we return None to indicate that no progress is
            # taking place. The calling routine then needs to wait for a new
            # frame to come in. See rtplot for an example of this.
            return None

        # set the internal frame pointer to the frame just read
        self.nframe = frameCount

        if nsats == -1 and synced == -1:
            # invalid time; pretend we are in 2000-01-01 taking one frame per
            # second, just so we can get something.
            tstamp = Time(51544 + self.nframe/DAYSEC, format='mjd', precision=9)
        else:
            time_string = '{}:{}:{}:{}:{:.7f}'.format(
                years, day_of_year, hours, mins, seconds+nanoseconds/1e9
                )
            tstamp = Time(time_string, format='yday', precision=9)

        # copy over the top-level header to avoid it becoming a reference
        # common to all MCCDs produced by the routine
        thead = self.thead.copy()
        thead['TIMSTAMP'] = (tstamp.isot, 'Raw frame timestamp, UTC')
        thead['MJDUTC'] = (tstamp.mjd, 'MJD(UTC) equivalent')
        if (nsats == -1 and synced == -1) or synced == 0:
            thead['GOODTIME'] = (False, 'Is TIMSTAMP thought to be OK?')
        else:
            thead['GOODTIME'] = (True, 'Is TIMSTAMP thought to be OK?')
        thead['NFRAME'] = (frameCount, 'Frame number')

        # second, the data bytes

        # build Windows-->CCDs-->MCCD
        CNAMS = ('1', '2', '3', '4', '5')
        QNAMS = ('E', 'F', 'G', 'H')


        # update the headers, initialise the CCDs
        cheads = {}
        ccds = Group(CCD)
        for nccd, (cnam, chead) in enumerate(zip(CNAMS, self.cheads)):

            # explicitly copy each header to avoid propagation of references
            ch = chead.copy()

            # use ncycle to determine whether this frame is really data as
            # opposed to intermediate blank frame
            isdata = self.nframe % ch['NCYCLE'] == 0 \
                     if 'NCYCLE' in ch else True
            ch['DSTATUS'] = (isdata, 'Data status when NCYCLE > 1')

            # Get time at centre of exposure
            tmid, texp, flag = self.timing(tstamp.mjd, frameCount, nccd)
            ch['MJDUTC'] = (tmid, 'MJD(UTC) at centre of exposure')
            ch['GOODTIME'] = (flag and thead['GOODTIME'], 'Is MJDUTC reliable?')
            ch['EXPTIME'] = (texp, 'Exposure time (secs)')

            # store for later recovery when creating the Windows
            cheads[cnam] = ch

            # Create the CCDs
            ccds[cnam] = CCD(Group(Window), HCM_NXTOT, HCM_NYTOT)

        # npixel points to the start pixel of the set of windows under
        # consideration
        npixel = 0
        for nwin in self.nwins:

            # get an example window for this set
            win = self.windows[nwin][0][0][0]
            nchunk = 20*win.nx*win.ny

            # allwins contains data of all windows 1 or 2
            allwins = frame[npixel:npixel+nchunk]

            # re-format the data as a 4D array indexed by (ccd,window,y,x)

            # get number of samples per pixel, with a default of 4
            # pixel data order depends on number of samples, and
            # the data prior to implementing sampling is the same as
            # nsamps = 4
            nsamps = self.header.get('ESO DET NSAMP', 4)
            if nsamps == 4:
                data = as_strided(
                    allwins, strides=(8, 2, 40*win.nx, 40),
                    shape=(5, 4, win.ny, win.nx)
                )
            else:
                # with 1 sample per pixel the data cannot be simply
                # re-viewed but a copy must be made
                data = as_strided(
                    allwins, strides=(32, 2, 160, 8),
                    shape=(5, 4, len(frame)//80, 4)
                ).reshape(5, 4, 512, 1024)

            # now build the Windows
            for nccd, cnam in enumerate(CNAMS):
                for nquad, qnam in enumerate(QNAMS):
                    wnam = '{:s}{:d}'.format(qnam,nwin+1)

                    # recover the window and flip parameters
                    win, flip_axes = self.windows[nwin][nccd][nquad]

                    # get the window data and apply flips. using older
                    # flipud and fliplr rather than more modern flip
                    # to avoid a need to upgrade numpy as flip only
                    # came in version 1.12 and many linux distros
                    # are behind.
                    windata = data[nccd, nquad]
                    for ax in flip_axes:
                        # windata = np.flip(windata, ax)
                        if ax == 0:
                            windata = np.flipud(windata)
                        elif ax == 1:
                            windata = np.fliplr(windata)
                        else:
                            raise ValueError(
                                'can only flip on axis 0 or 1, but got = {:d}'.format(ax)
                            )

                    if nwin == 0 and nccd == 0 and nquad == 0:
                        # store CCD header in the first Winhead
                        win.update(cheads[cnam])

                    # finally store as a Window, converting to 32-bit
                    # floats to avoid problems down the line.
                    ccds[cnam][wnam] = Window(
                        win, windata.astype(np.float32)
                        )

            # move pointer on for next set of windows
            npixel += nchunk

        # create the MCCD
        mccd = MCCD(ccds, thead)

        # update the frame counter for the next call
        self.nframe += 1

        # show that we have read something
        self.first = False

        return mccd

    def seek_frame(self, n):
        """
        Moves pointer position to the start of frame n and
        updates the 'nframe' attribute appropriately.

        This is not implemented for the server.
        """
        if self.server:
            raise NotImplementedError('no seek_frame in the case of server access')
        else:
            # move pointer to the start of the last complete frame
            # and adjust nframe
            self._ffile.seek(self._hbytes+self._framesize*(n-1))
            self.nframe = n

    def seek_last(self):
        """
        In the case of a local file, this sets the pointer position
        to the start of the final complete, and updates the 'nframe'
        attribute appropriately. This is not implemented for the server.
        It returns the number of the frame it moves to
        """
        if self.server:
            raise NotImplementedError('no seek_last in the case of server access')
        else:
            # First seek the end of file
            self._ffile.seek(0,2)

            # compute byte offset at this point and thus the total number of
            # complete frames
            ntot = (self._ffile.tell()-self._hbytes) // self._framesize

            # move pointer
            self.seek_frame(ntot)

        return ntot

class Rtime (Rhead):
    """Callable, iterable object to generate timing data only from HiPERCAM raw data files.

    This is similar to Rdata except it only reads the timing data of each
    frame in the case of local disk files and it does no processing of raw
    data other than the timing in any circumstance.
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
        nframe == 1 gives the first frame.

        Arguments::

           nframe : int | none
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted.

        Returns:: (tstamp, tinfo) or None if you try to read a non-existent
        frame (useful if you are waiting from frames to be added). Here
        'tstamp' is an astropy.time.Time equivalent to the GPS timestamp
        associated with the frame, with no correction, while 'tinfo' is a
        5-elements tuple, one for each CCD listing for each one the MJD at
        mid-exposure, an exposure time in seconds and a flag to indicate
        whether the time is thought reliable, which is only the case if
        real data comes with the frame.
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
                # if we are trying to access a file that has not yet been
                # written, 0 bytes will be returned. In this case we return
                # with None. NB we do not update self.nframe in this case.
                return None

            # extract the timing data
            tbytes = raw_bytes[-self.ntbytes:]

        elif self.nframe > self.ntotal():
            # We have attempted to access a non-existent frame
            return None

        else:
            # find the start of the timing data if necessary
            if reset:
                self._ffile.seek(self._hbytes+self._framesize*self.nframe-self.ntbytes)

            # read the timing data, skip forward to next set
            tbytes = self._ffile.read(self.ntbytes)
            self._ffile.seek(self._framesize-self.ntbytes, 1)

        # Interpret the timing bytes
        frameCount, timeStampCount, years, day_of_year, \
            hours, mins, seconds, nanoseconds, nsats, synced = decode_timing_bytes(tbytes)

        if nsats == -1 and synced == -1:
            # invalid time; pretend we are on 2000-01-01 taking one frame per second.
            tstamp = Time(51544 + self.nframe/86400., format='MJD',precision=9)
        else:
            time_string = '{}:{}:{}:{}:{:.7f}'.format(
                years, day_of_year, hours, mins, seconds+nanoseconds/1e9
                )
            tstamp = Time(time_string, format='yday',precision=9)

        tinfo = []
        for nccd in range(5):
            tinfo.append(self.timing(tstamp.mjd,self.nframe,nccd))

        # update the internal frame counter
        self.nframe += 1

        # Return timing data
        return (tstamp, tuple(tinfo))

def decode_timing_bytes(tbytes):
    """Decode the timing bytes tacked onto the end of every HiPERCAM frame in the
    3D FITS file.

    The timing bytes are encoded as 8x4-byte little-endian unsigned integers
    followed by 2x2-byte little-endian . On saving to FITS they are mangled
    because:

      1) they are split into 2-byte uints
      2) 2**15 is subtracted from each to map to 2-byte signed ints
      3) These are written in big-endian order to FITS

    This routine reverses this process to recover the original timestamp tuple
    by

      1) unpacking the timing bytes as big-endian 16-bit signed integers
      2) adding 32768 to each value
      3) packing these values as little-endian 16-bit unsigned integers
      4) unpacking the result as 32-bit, little-endian unsigned integers

    Parameters
    ----------

       tbytes: (bytes)
           a Python bytes object which contains the timestamp bytes as written
           in the FITS file

    Returns
    --------

      timestamp : (tuple)
          a tuple containing (frameCount, timeStampCount, years, day_of_year,
          hours, mins, seconds, nanoseconds, nsats, synced) values.

    [Straight from Stu Littlefair's original code.]

    """

    # format with 36 timing bytes (last two of which are ignored).
    # The -36 is to try to cope with cases where more than 36 bytes come through.
    buf = struct.pack('<HHHHHHHHHHHHHHHHH',
                     *(val + BZERO for val in struct.unpack('>hhhhhhhhhhhhhhhhh',
                                                             tbytes[-36:-2])))
    return struct.unpack('<IIIIIIIIbb', buf)

class HendError(HipercamError):
    """
    Exception for the standard way to reach the end of a HiPERCAM raw data
    file (attempt to read out-of-range frame). This allows the iterator to die
    silently in this case while raising exceptions for less anticipated cases.
    """
    pass
