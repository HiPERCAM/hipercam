"""
Section for handling the raw HiPERCAM data files. These are FITS files
with a single HDU containing header information on the window format amonst
other things with the data in a "FITS cube" of unusual nature in that it has
dimensions of nframe x 1 x lots, where "lots" contains all pixels of all
windows of all CCDs plus timing bytes at the end. The pixels come in an odd
order (pixel 1 of each window and then the same for each CCD etc), and the
main work of the code is in disentangling them to produce :class:MCCD objects.

"""

import json
import os
import struct
import sys
import warnings

# temporary
import matplotlib.pyplot as plt
import numpy as np
import websocket
from astropy.io import fits
from astropy.time import Time, TimeDelta
from numpy.lib.stride_tricks import as_strided

from hipercam import (
    CCD,
    HRAW,
    MCCD,
    Group,
    Header,
    HipercamError,
    Window,
    Winhead,
    fday_to_hms,
    gregorian_to_mjd,
    mjd_to_gregorian,
)
from hipercam.utils import add_extension

__all__ = [
    "Rhead",
    "Rdata",
    "Rtime",
    "Rtbytes",
    "HCM_NXTOT",
    "HCM_NYTOT",
    "HCM_NPSCAN",
    "HCM_NOSCAN",
]

# Set the URL of the FileServer from the environment or default to localhost.
URL = "ws://" + os.environ.get("HIPERCAM_DEFAULT_URL", "localhost:8007/")

# FITS offset to convert from signed to unsigned 16-bit ints and vice versa
BZERO = 1 << 15

# number of seconds in a day
DAYSEC = 86400.0

# Maximum dimensions of HiPERCAM CCDs, imaging area
HCM_NXTOT, HCM_NYTOT = 2048, 1024

# Number of unbinned columns added to each window
# when a prescan is present
HCM_NPSCAN = 50

# Number of unbinned rows added to each window
# when an overscan is present
HCM_NOSCAN = 8

# bottom-left coordinate of quadrants
LLX = {"E": 1, "F": HCM_NXTOT // 2 + 1, "G": HCM_NXTOT // 2 + 1, "H": 1}
LLY = {"E": 1, "F": 1, "G": HCM_NYTOT // 2 + 1, "H": HCM_NYTOT // 2 + 1}

# the direction increasing X- or Y-start moves the window in a quadrant
X_DIRN = {"E": 1, "F": -1, "G": -1, "H": 1}
Y_DIRN = {"E": 1, "F": 1, "G": -1, "H": -1}

# whether we need to add size of window to find bottom-left coordinate of
# window
ADD_YSIZES = {"E": 0, "F": 0, "G": 1, "H": 1}
ADD_XSIZES = {"E": 0, "F": 1, "G": 1, "H": 0}

# some axes need flipping, since the data is read out such that the
# bottom-left is not always the first pixel. Which axes need flipping
# to get data in correct orientation uses numpy C-convention, 1=y,
# 0=x. These are fed to numpy.flip later
FLIP_AXES = {"E": (), "F": (1,), "G": (1, 0), "H": (0,)}

# Code to account for reflections of the g [CCD 1, zero offset] and z [CCD 4]
# CCDs, and potential 180 degree rotations (In Feb 2018 GTC run, the u-band
# CCD was 180 degrees rotated cf intended).
REFLECTED = (1, 4)
QNAMS = ("E", "F", "G", "H")
QNAMS_REFLECT = {"E": "F", "F": "E", "G": "H", "H": "G"}
QNAMS_ROTATED = {"E": "G", "F": "H", "G": "E", "H": "F"}

ESSENTIAL = "MODE"


class Rhead:
    """Reads an interprets header information from a HiPERCAM run (3D FITS file)
    and generates the window formats needed to interpret the data. The file is
    kept open while the Rhead exists. The file can be on the local disk or
    obtained via the fileserver.

    The following attributes are set::

       cheads : list of Header objects
           the headers for each CCD which contain any info needed per CCD, above
           all the 'skip' numbers, which are used when creating MCCD frames.

       clear, dummy, oscan, prescan  : bool
           flags to say whether or not the following were enabled:
               clears between exposures
               dummy output
               overscan (extra pixels in Y)
               prescan (extra pixels in X)

       fname : string
           the file name

       header : astropy.io.fits.Header
           the header of the FITS file.

       mode : string
           the readout mode used, read from 'ESO DET READ CURNAME'

       nskips : tuple
           the NSKIP parameters for each CCD (5 of them)

       ntbytes : int
           number of timing bytes, deduced from the dimensions of the FITS cube
           (NAXIS1 contains all pixels and timing bytes) compared to the window
           dimensions.

       nwins : tuple
           integer indices of the windows for any one quadrant. This is (0,) for
           one window modes, (0,1) for two window modes.

       tdead, tdelta, toff1, toff2, toff3, toff4 : float
           set of parameters used to calculate exact times. tdead is the
           dead time between exposures; tdelta is used to multiply NSKIP
           for each CCD; toff1 etc are offsets and lengths used to get the
           mid-exposure times and exposure durations.

       thead : Header
           the top-level header that will be used to create MCCDs.

       windows : list of list of list of (Winhead, flip) tuples
           windows[nccd][nquad][nwin] returns the Winhead and set of axes for
           window nwin of quadrant nquad of CCD nccd. nwin indexes the windows
           of a given quadrant (at most there are 2, so nwin = 0 or 1). nquad
           indexes the quadrants in EFGH order. nccd index the 5 CCDs. 'flip'
           is a tuple of axes to be flipped using numpy.flip to account for
           where the CCD is read and reflections. These are not in the end the
           final windows that will end in the MCCD because if there are pre-scan
           or over-scan present, we need to split them off making more windows.

       wforms : tuple of strings
           formats of each set of windows as strings designed to match the input
           expected for hdriver.

       xbin : int
           X-binning factor

       ybin : int
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
              environment variable "HIPERCAM_DEFAULT_URL" (pre-fixed by
              'ws://'), or, if that is not set, "ws://localhost:8007/".  The
              server here is Stu Littlefair's Python-based server that
              defaults to port 8007.

           server : (bool)
              True to access the data from a server. It uses a URL set in an
              environment variable "HIPERCAM_DEFAULT_URL" in this instance.

           full   : bool
              Flag controlling the amount of header detail to gather (as a
              time saver) True for detail, False for the bare minimum, the
              latter might be useful for large numbers of very small format
              images where you don't to waste effort. The full headers are
              extensive (~2300 lines, mostly voltages == 184kB) and it's hard
              to see why one would ever want them all.

        """

        # store the file name and whether server being used
        self.fname = fname
        self.server = server
        self.full = full

        if server:
            # open socket connection to server
            self._ws = websocket.create_connection(URL + fname)
            status = json.loads(self._ws.recv())["status"]
            if status == "no such run":
                raise HipercamError("Run not found: {}".format(fname))

            # read the header from the server
            data = json.dumps(dict(action="get_hdr"))
            self._ws.send(data)
            hd = self.header = Header(fits.Header.fromstring(self._ws.recv()))

            nsamps = self.header.get("ESO DET NSAMP", 1)
            self._framesize = (
                18
                + (
                    self.header["ESO DET ACQ1 WIN NX"]
                    * self.header["ESO DET ACQ1 WIN NY"]
                )
                // nsamps
            )
        else:
            # open the file
            self._ffile = open(add_extension(fname, HRAW), "rb")

            # read the header
            hd = self.header = Header(fits.Header.fromfile(self._ffile))

            # store number of bytes read from header. used to offset later
            # requests for data
            self._hbytes = self._ffile.tell()
            self._framesize = self.header["NAXIS1"]

        # calculate the framesize in bytes
        bitpix = abs(hd["BITPIX"])
        self._framesize *= bitpix // 8

        # store the scaling and offset
        self._bscale = hd["BSCALE"]
        self._bzero = hd["BZERO"]

        # create the top-level header
        self.thead = Header()

        # set the mode, one or two windows per quadrant, drift etc.  This is
        # essential.
        self.mode = hd["ESO DET READ CURNAME"]
        if self.full:
            self.thead["MODE"] = (self.mode, "HiPERCAM readout mode")
        self.drift = self.mode.startswith("Drift")
        if self.drift:
            self.ndwins = hd.get("ESO DET DRIFT NWINS", 0)
        else:
            self.ndwins = 0

        # check for overscan / prescan / clear
        self.clear = hd["ESO DET CLRCCD"]
        self.dummy = hd["ESO DET DUMMY"]
        self.oscan = hd["ESO DET INCOVSCY"]
        self.pscan = hd["ESO DET INCPRSCX"]
        self.thead["CLEAR"] = (self.clear, "Clear enabled")

        # store the NSKIP values for each CCD. These are needed for exact
        # timing for each CCD.
        self.nskips = (
            hd["ESO DET NSKIPS1"],
            hd["ESO DET NSKIPS2"],
            hd["ESO DET NSKIPS3"],
            hd["ESO DET NSKIPS4"],
            hd["ESO DET NSKIPS5"],
        )

        # set timing parameters toff1, toff2, toff3, toff4, tdelta, tdead
        # which are also needed for exact timing. Meanings:
        #
        # toff1  -- constant offset from timestamp for first frame
        # toff2  -- constant offset from timestamp for all other frames
        # toff3  -- constant part of exposure time for first frame
        # toff4  -- constant part of exposure time for all other frames
        # tdelta -- multiplier of nskip (same for all frames)
        # tdead  -- dead time between frames.

        E = hd["ESO DET TDELAY"]

        if self.clear:
            # Clear mode. NB the parameter 'TREAD' is actually a combination
            # of a frame-transfer and a read (F+R in the docs), hence it
            # is called FR here
            W = hd["ESO DET TCLEAR"]
            FR = hd["ESO DET TREAD"]

            self.tdelta = 1e-3 * (FR + W + E)
            self.toff1 = self.toff2 = 1.0e-3 * E / 2.0
            self.toff3 = self.toff4 = 1.0e-3 * E
            self.tdead = 1.0e-3 * (FR + W)

        elif self.drift:
            # Drift mode
            LD = hd["ESO DET DRIFT TLINEDUMP"]
            LS = hd["ESO DET DRIFT TLINESHIFT"]
            ND = self.ndwins
            R = hd["ESO DET TREAD"]

            # No NSKIP here, so tdelta irrelevant. in drift
            # mode, TREAD is the read time (no frame transfer)
            self.tdelta = 0.0
            self.toff1 = self.toff2 = 1.0e-3 * (
                E + LS + (E + LD + R) / 2.0 - ND * (E + R + LS + LD)
            )
            self.toff3 = self.toff4 = 1.0e-3 * (E + R + LD)
            self.tdead = 1.0e-3 * LS

        else:
            # No clear mode, there are zero sec dummy
            # 'wipes' so no 'W' appears here. TREAD
            # is back to being the sum of the read
            # and frame transfer times
            FR = hd["ESO DET TREAD"]
            F = hd["ESO DET TFT"]
            R = FR - F

            self.tdelta = 1.0e-3 * (FR + E)
            self.toff1 = 1.0e-3 * E / 2.0
            self.toff2 = 1.0e-3 * (E - R) / 2.0
            self.toff3 = 1.0e-3 * E
            self.toff4 = 1.0e-3 * (R + E)
            self.tdead = 1.0e-3 * F

        # binning factors
        self.xbin = hd["ESO DET BINX1"]
        self.ybin = hd["ESO DET BINY1"]

        if self.pscan:
            if HCM_NPSCAN % self.xbin == 0:
                # number of binned prescan pixels.
                self.npscan = HCM_NPSCAN // self.xbin

            else:
                raise HipercamError(
                    (
                        "Use of X-binning factors which are not divisors of"
                        " {:d} with overscan enabled is not supported"
                    ).format(HCM_NPSCAN)
                )

        if self.oscan:
            if HCM_NOSCAN % self.ybin == 0:
                # number of binned prescan pixels
                self.noscan = HCM_NOSCAN // self.ybin

            else:
                raise HipercamError(
                    (
                        "Use of Y-binning factors which are not divisors of"
                        " {:d} with overscan is not supported"
                    ).format(HCM_NOSCAN)
                )

        if self.mode.startswith("TwoWindow"):
            self.nwins = (0, 1)
        else:
            self.nwins = (0,)

        # Rotation: in the Feb 2018 GTC run, CCD 1 was physically
        # rotated 180 degrees. Stu accounted for this in the software.
        # This part is to account for it in the pipeline, defaulting
        # to no rotation by default.
        self.rotate = []
        for nccd in range(5):
            self.rotate.append(hd.get("ESO DET ROTATE{:d}".format(nccd + 1), False))

        # Build window parameters for each window of each quadrant of each
        # CCDs. At this stage, windows with pre- and/or over-scans are treated
        # as single blocks, data and scans. They are split up in Rdata.
        self.windows = []
        for nwin in self.nwins:
            self.windows.append([])

            # first set the windows of each CCD
            for nccd in range(5):
                self.windows[-1].append([])
                for qnam in QNAMS:
                    if nccd in REFLECTED:
                        # reflections for g (1) and z (4)
                        qnam = QNAMS_REFLECT[qnam]

                    if self.rotate[nccd]:
                        # if CCD rotated, swap E<-->G, F<-->H, etc
                        qnam = QNAMS_ROTATED[qnam]

                    if self.mode.startswith("FullFrame"):
                        # Full frame mode, account for pre- and over-scans
                        # in the full frame.
                        nx = (
                            HCM_NXTOT // 2 + HCM_NPSCAN
                            if self.pscan
                            else HCM_NXTOT // 2
                        ) // self.xbin

                        ny = (
                            HCM_NYTOT // 2 + HCM_NOSCAN
                            if self.oscan
                            else HCM_NYTOT // 2
                        ) // self.ybin

                        llx = LLX[qnam]
                        lly = LLY[qnam]

                    elif self.mode.startswith("OneWindow") or self.mode.startswith(
                        "TwoWindow"
                    ):
                        # windowed mode
                        winID = "ESO DET WIN{} ".format(nwin + 1)
                        nx = hd[winID + "NX"] // self.xbin
                        ny = hd[winID + "NY"] // self.ybin

                        if self.pscan:
                            # account for extra columns when a pre-scan is
                            # present
                            nx += self.npscan

                        if self.oscan:
                            # account for extra rows when an over-scan is
                            # present
                            ny += self.noscan

                        win_xs = hd[winID + "XS{}".format(qnam)]
                        win_nx = hd[winID + "NX"]
                        win_ys = hd[winID + "YS"]
                        win_ny = hd[winID + "NY"]
                        llx = (
                            LLX[qnam]
                            + X_DIRN[qnam] * win_xs
                            + ADD_XSIZES[qnam] * (HCM_NXTOT // 2 - win_nx)
                        )
                        lly = (
                            LLY[qnam]
                            + Y_DIRN[qnam] * win_ys
                            + ADD_YSIZES[qnam] * (HCM_NYTOT // 2 - win_ny)
                        )

                    elif self.mode.startswith("Drift"):
                        nx = hd["ESO DET DRWIN NX"] // self.xbin
                        ny = hd["ESO DET DRWIN NY"] // self.ybin

                        if self.pscan:
                            raise HipercamError("prescans with drift mode undefined")

                        win_xs = hd["ESO DET DRWIN XS{}".format(qnam)]
                        win_nx = hd["ESO DET DRWIN NX"]
                        win_ys = hd["ESO DET DRWIN YS"]
                        win_ny = hd["ESO DET DRWIN NY"]
                        llx = (
                            LLX[qnam]
                            + X_DIRN[qnam] * win_xs
                            + ADD_XSIZES[qnam] * (HCM_NXTOT // 2 - win_nx)
                        )
                        lly = (
                            LLY[qnam]
                            + Y_DIRN[qnam] * win_ys
                            + ADD_YSIZES[qnam] * (HCM_NYTOT // 2 - win_ny)
                        )

                    else:
                        msg = "mode {} not currently supported".format(mode)
                        raise ValueError(msg)

                    # store the window and the axes to flip
                    self.windows[-1][-1].append(
                        (
                            Winhead(llx, lly, nx, ny, self.xbin, self.ybin, ""),
                            FLIP_AXES[qnam],
                        )
                    )

        # set strings for logging purposes giving the settings used in hdriver
        self.wforms = []

        if self.mode.startswith("FullFrame"):
            self.wforms.append("Full&nbsp;Frame")

        elif self.mode.startswith("OneWindow") or self.mode.startswith("TwoWindow"):
            winID = "ESO DET WIN1 "
            ys = hd[winID + "YS"] + 1
            nx = hd[winID + "NX"]
            ny = hd[winID + "NY"]

            try:
                xsll = hd[winID + "XSLL"]
                xslr = hd[winID + "XSLR"]
                xsul = hd[winID + "XSUL"]
                xsur = hd[winID + "XSUR"]
            except:
                # Fallback
                xsll = hd[winID + "XSE"] + 1
                xsul = hd[winID + "XSH"] + 1
                xslr = HCM_NXTOT + 1 - hd[winID + "XSF"] - nx
                xsur = HCM_NXTOT + 1 - hd[winID + "XSG"] - nx

            self.wforms.append(
                "{:d},{:d},{:d},{:d},{:d},{:d},{:d}".format(
                    xsll, xslr, xsul, xsur, ys, nx, ny
                )
            )

            if self.mode.startswith("TwoWindow"):
                winID = "ESO DET WIN2 "
                ys = hd[winID + "YS"] + 1
                nx = hd[winID + "NX"]
                ny = hd[winID + "NY"]

                try:
                    xsll = hd[winID + "XSLL"]
                    xslr = hd[winID + "XSLR"]
                    xsul = hd[winID + "XSUL"]
                    xsur = hd[winID + "XSUR"]
                except:
                    # fallback option
                    xsll = hd[winID + "XSE"] + 1
                    xsul = hd[winID + "XSH"] + 1
                    xslr = HCM_NXTOT + 1 - hd[winID + "XSF"] - nx
                    xsur = HCM_NXTOT + 1 - hd[winID + "XSG"] - nx

                self.wforms.append(
                    "{:d},{:d},{:d},{:d},{:d},{:d},{:d}".format(
                        xsll, xslr, xsul, xsur, ys, nx, ny
                    )
                )

        elif self.mode.startswith("Drift"):
            winID = "ESO DET DRWIN "
            ys = hd[winID + "YS"] + 1
            nx = hd[winID + "NX"]
            ny = hd[winID + "NY"]

            try:
                xsl = hd[winID + "XSL"]
                xsr = hd[winID + "XSR"]
            except:
                # Fallback option
                xsl = hd[winID + "XSE"] + 1
                xsr = HCM_NXTOT + 1 - hd[winID + "XSF"] - nx

            self.wforms.append("{:d},{:d},{:d},{:d},{:d}".format(xsl, xsr, ys, nx, ny))

        # compute the total number of pixels (making use of symmetry between
        # quadrants and CCDs, hence factor 20)
        npixels = 0
        for nwin in self.nwins:
            win = self.windows[nwin][0][0][0]
            npixels += 20 * win.nx * win.ny

        # any extra bytes above those needed for the data are timing
        self.ntbytes = self._framesize - (bitpix * npixels) // 8

        # Build (more) header info
        if self.full:
            if "DATE" in hd:
                self.thead["DATE"] = hd.get_full("DATE")

            if "OBJECT" in hd:
                self.thead["OBJECT"] = hd.get_full("OBJECT")

            if "FILTERS" in hd:
                self.thead["FILTERS"] = hd.get_full("FILTERS")

            if "PI" in hd:
                self.thead["PI"] = hd.get_full("PI")

            if "OBSERVER" in hd:
                self.thead["OBSERVER"] = hd.get_full("OBSERVER")

            if "RA" in hd:
                self.thead["RA"] = hd.get_full("RA")

            if "DEC" in hd:
                self.thead["DEC"] = hd.get_full("DEC")

            if "RADEG" in hd:
                self.thead["RADEG"] = hd.get_full("RADEG")

            if "DECDEG" in hd:
                self.thead["DECDEG"] = hd.get_full("DECDEG")

            if "INSTRPA" in hd:
                self.thead["INSTRPA"] = hd.get_full("INSTRPA")

            if "IMAGETYP" in hd:
                self.thead["IMAGETYP"] = hd.get_full("IMAGETYP")

            if "ESO DET GPS" in hd:
                self.thead["GPS"] = hd.get_full("ESO DET GPS")

            if "INSTRUME" in hd:
                self.thead["INSTRUME"] = hd.get_full("INSTRUME")

            if "TELESCOP" in hd:
                self.thead["TELESCOP"] = hd.get_full("TELESCOP")

            self.thead["XBIN"] = hd.get_full("ESO DET BINX1")
            self.thead["YBIN"] = hd.get_full("ESO DET BINY1")
            self.thead["SPEED"] = hd.get_full("ESO DET SPEED")

        if "EXPTIME" in hd:
            self.thead["EXPTIME"] = hd.get_full("EXPTIME")

        # Header per CCD. These are modified per CCD in Rdata
        self.cheads = []
        for n in range(5):
            chead = Header()

            # Essential items
            hnam = "ESO DET NSKIPS{:d}".format(n + 1)
            chead["NCYCLE"] = (hd.get(hnam, 0) + 1, "readout cycle period (NSKIP+1)")

            # Nice-if-you-can-get-them items
            if self.full:
                # whether this CCD has gone through a reflection
                hnam = "ESO DET REFLECT{:d}".format(n + 1)
                if hnam in hd:
                    chead["REFLECT"] = (hd[hnam], "is image reflected")

                # readout noise
                hnam = "ESO DET CHIP{:d} RON".format(n + 1)
                if hnam in hd:
                    chead["RONOISE"] = self.header.get_full(hnam)

                # gain
                hnam = "ESO DET CHIP{:d} GAIN".format(n + 1)
                if hnam in self.header:
                    chead["GAIN"] = self.header.get_full(hnam)

            self.cheads.append(chead)

        # end of constructor / initialiser

    def ntotal(self):
        """
        Returns the total number of complete frames.
        """
        if self.server:
            data = json.dumps(dict(action="get_nframes"))
            self._ws.send(data)
            raw_bytes = self._ws.recv()
            if len(raw_bytes) == 0:
                raise HipercamError(
                    "failed to get total number of frames from server; 0 bytes returned"
                )
            d = eval(raw_bytes)
            ntot = d["nframes"]
        else:
            # the commented out section was an attempt to work out the number
            # of frames even when the file was growing. This fails owing the
            # 2880-byte FITS standard when one has small framesizes, thus I am
            # using NAXIS3 instead
            ntot = self.header["NAXIS3"]

        return ntot

    def tinfo(self):
        """Returns with timing info for the run used to create the Rhead, namely the
        exposure times for each CCD appropriate for a "standard" exposure
        (i.e. not the first), the offsets to get to the mid-exposure from a
        timestamp, again for a "standard" exposure, the cycle length for each
        CCD (NSUB+1) and the dead time between frames. This is intended for
        use in generating logs.

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
            toff, texp, flag = self.timing(2 * (nskip + 1), nccd)
            toffs.append(DAYSEC * toff)
            texps.append(texp)
            ncycs.append(nskip + 1)

        return (tuple(texps), tuple(toffs), tuple(ncycs), self.tdead)

    def timing(self, nframe, nccd):
        """Returns timing data for a particular CCD and frame.

        Arguments::

           nframe : int
              the frame number starting at 1.

           nccd : int
              the CCD number, starting at 0

        Returns: (tmid, exptime, flag) where 'tmid' is the offset to
        be applied to the time stamp (in days) to get to the
        mid-exposure pointof the CCD, 'exptime' is the exposure time
        in seconds, and 'flag' is a bool to indicate whether the time
        is thought good or not.

        """

        nskip = self.nskips[nccd]

        flag = nframe % (nskip + 1) == 0
        if nframe == nskip + 1:
            # special case for first frame of data
            tmid = (self.toff1 - self.tdelta * nskip / 2) / DAYSEC
            texp = self.tdelta * nskip + self.toff3
        else:
            # regular times
            tmid = (self.toff2 - self.tdelta * nskip / 2) / DAYSEC
            texp = self.tdelta * nskip + self.toff4

        return (tmid, texp, flag)

    def __del__(self):
        """Destructor closes the file or web socket"""
        if self.server:
            if hasattr(self, "_ws"):
                self._ws.close()
        else:
            if hasattr(self, "_ffile"):
                self._ffile.close()

    # Want to run this as a context manager
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.__del__()


class Rdata(Rhead):
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

    def __init__(self, fname, nframe=1, server=False, full=True):
        """Connects to a raw HiPERCAM FITS file for reading. The file is kept
        open.  The Rdata object can then generate MCCD objects through being
        called as a function or iterator.

        Arguments::

           fname : string
              run name, e.g. 'run036'.

           nframe : int
              the frame number to read first [1 is the first]. This
              initialises an attribute of the same name that is used when
              reading frames sequentially. nframe=0 is an indication to set an
              attribute 'last' = True to indicate that it should always try to
              access the last frame. nframe=-10 means try to get the
              10th-from-last frame, a cludge to get round an acquisition bug.

           server : bool
              True/False for server vs local disk access. Server access goes
              through a websocket.  It uses a base URL taken from the
              environment variable "HIPERCAM_DEFAULT_URL", or, if that is not
              set, "ws://localhost:8007/".  The server here is Stu Littlefair's
              Python-based server that defaults to port 8007.

           full : bool
              Controls the amount of header info returned with the data. True
              for as much detail as possible; False for minimal headers which
              require less resources downstream.

        """

        # read the header
        Rhead.__init__(self, fname, server, full)

        # flag to indicate should always try to get the last frame
        self.last = nframe == 0

        # store initial value of first
        self.ffirst = nframe

        # set flag to indicate first time through
        self.first = True

        if server:
            # get the specified frame, or - if last is True - the frame after the current last
            self.nframe = nframe if not self.last else self.ntotal() + 1
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
        except HendError:
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

        Parameters:

           nframe : int
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted, unless
              self.nframe = 0 in which case it will try to get the most recent
              frame, whatever that is.

        Apart from reading the raw bytes, the main job of this routine is to
        divide up and re-package the bytes read into Windows suitable for
        constructing CCD objects.

        If access via a server is requested, it is assumed that the file being
        accessed could be being added to.
        """

        if self.server:
            # access frames via the server

            # build the request
            if nframe == 0 or self.last:
                # just want the last complete frame.
                # If this is not at least one further on than the latest
                # frame we have, we return with None to indicate no progress.
                if self.ntotal() < self.nframe:
                    return None
                # Note that further down we check the frame counter against what we were
                # hoping for (at least 1 further on) to check for issues.
                request = json.dumps(dict(action="get_last"))

            elif self.ffirst < 0:
                # this case we are trying to get a frame -self.ffirst from
                # being the last frame.
                nget = self.ntotal() + self.ffirst
                if nget < 1:
                    return None

                if nget == self.nframe:
                    request = json.dumps(dict(action="get_next"))
                elif nget < self.nframe:
                    # we have already read frame nget, so return None for now
                    return None
                else:
                    request = json.dumps(dict(action="get_frame", frame_number=nget))
                    self.nframe = nget

            elif nframe is None:
                # in this case, on the first time through we need to
                # explicitly request the frame, otherwise we can just read
                # from where we are
                if self.first:
                    request = json.dumps(
                        dict(action="get_frame", frame_number=self.nframe)
                    )
                else:
                    request = json.dumps(dict(action="get_next"))

            else:
                # a particular frame number is being requested. Check whether
                # we have to request it explicitly or whether we can just get
                # the next one
                if self.nframe == nframe:
                    request = json.dumps(dict(action="get_next"))
                else:
                    request = json.dumps(dict(action="get_frame", frame_number=nframe))
                    self.nframe = nframe

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
            frame = np.fromstring(raw_bytes[: -self.ntbytes], ">u2")
            frame += BZERO
            tbytes = raw_bytes[-self.ntbytes :]

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
                ntot = self.seek_last()

                # move to correct place
                self.seek_frame(nframe)

            # read in frame and then the timing data, correcting the frame for
            # the standard FITS BZERO offset. At this stage we have the data
            # as unsigned 2-byte ints
            frame = np.fromfile(
                self._ffile, ">u2", (self._framesize - self.ntbytes) // 2
            )
            frame += BZERO
            tbytes = self._ffile.read(self.ntbytes)
            if len(tbytes) != self.ntbytes:
                raise HendError("failed to read frame from disk file")

        ##############################################################
        #
        # We now have all the pixels (in 'frame') and timing bytes (in
        # 'tbytes') of the frame of interest read into a 1D ndarray of
        # unsigned, 2-byte integers. Now interpret them.
        #
        ##############################################################

        # First the timing bytes. The frameCount starts from 0 so we
        # we add one to it
        (
            frameCount,
            timeStampCount,
            years,
            day_of_year,
            hours,
            mins,
            seconds,
            nanoseconds,
            nsats,
            synced,
        ) = htimer(tbytes)
        frameCount += 1

        if (
            self.server
            and (nframe == 0 or self.last or self.ffirst < 0)
            and self.nframe > frameCount
        ):
            # server access, trying to get the last complete frame or some
            # frames from the end. If the frame just read in is the same as
            # the one before (i.e. self.nframe > frameCount), we return None
            # to indicate that no progress is taking place. The calling
            # routine then needs to wait for a new frame to come in. See
            # rtplot for an example of this.

            # note that this should never be reached, since we should have returned
            # None above if the file was not ready yet.
            print(
                "  >>> WARNING: no progress in frame number, this should not happen, contact SL",
                file=sys.stderr,
            )
            return None

        elif not self.server and frameCount != self.nframe:
            if frameCount == self.nframe + 1:
                print(
                    "   >>> WARNING: frame count mis-match; a frame seems to have been dropped",
                    file=sys.stderr,
                )
            else:
                print(
                    "   >>> WARNING: frame count mis-match; {:d} frames seem to have been dropped".format(
                        frameCount - self.nframe
                    ),
                    file=sys.stderr,
                )

        # set the internal frame pointer to the frame just read
        self.nframe = frameCount

        time_string = "{}:{}:{}:{}:{:.9f}".format(
            years, day_of_year, hours, mins, seconds + nanoseconds / 1e9
        )
        try:
            imjd = gregorian_to_mjd(years, 1, 1) + day_of_year - 1
            fday = (hours + mins / 60 + (seconds + nanoseconds / 1e9) / 3600) / 24
        except ValueError:
            warnings.warn("Bad timestamp: " + time_string)
            imjd = 51544
            fday = self.nframe / DAYSEC

            # this will mark it as unreliable down below
            synced = 0

        # copy over the top-level header to avoid it becoming a reference
        # common to all MCCDs produced by the routine
        thead = self.thead.copy()

        year, month, day = mjd_to_gregorian(imjd)
        hour, minute, second = fday_to_hms(fday)
        thead["TIMSTAMP"] = (
            "{:4d}-{:02d}-{:02d}T{:02d}:{:02d}:{:012.9f}".format(
                year, month, day, hour, minute, second
            ),
            "Raw frame timestamp, UTC",
        )
        thead["MJDINT"] = (imjd, "Integer part of MJD(UTC), raw timestamp")
        thead["MJDFRAC"] = (fday, "Fractional part of MJD(UTC), raw timestamp")
        thead["MJDUTC"] = (imjd + fday, "MJD(UTC)")
        if (nsats == -1 and synced == -1) or synced == 0:
            thead["GOODTIME"] = (False, "TIMSTAMP OK?")
        else:
            thead["GOODTIME"] = (True, "TIMSTAMP OK?")
        thead["NFRAME"] = (frameCount, "Frame number")

        # second, the data bytes

        # build Windows-->CCDs-->MCCD
        CNAMS = ("1", "2", "3", "4", "5")

        # update the headers, initialise the CCDs
        cheads = {}
        ccds = Group(CCD)
        for nccd, (cnam, chead) in enumerate(zip(CNAMS, self.cheads)):
            # Explicitly copy each header to avoid propagation of references
            ch = chead.copy()

            # Determine whether the frame really contains data as
            # opposed to initial junk frames in the case of drift
            # mode or intermediate junk frames in the case of NSKIP>0
            ch["DSTATUS"] = (
                (self.nframe > self.ndwins) and (self.nframe % ch["NCYCLE"] == 0),
                "Valid data or not",
            )

            # Get time at centre of exposure. Some care here to store a
            # super-precise version of the time as a string along with the MJD
            # -- the latter suffers round-off at the 0.5 microsecond level in
            # double precision. The precise one allows better downstream
            # precision by offsetting ('toffset' in reduce).
            tmid, texp, flag = self.timing(frameCount, nccd)
            fdaymid = fday + tmid
            if fdaymid >= 1.0:
                fdaymid -= 1
                imjdmid = imjd + 1
            else:
                imjdmid = imjd

            ch["MJDINT"] = (imjdmid, "Integer part of MJD(UTC), mid-exposure")
            ch["MJDFRAC"] = (fdaymid, "Fractional part of MJD(UTC), mid-exposure")
            ch["MJDUTC"] = (imjdmid + fdaymid, "MJD(UTC) mid-exposure")
            ch["GOODTIME"] = (flag and thead["GOODTIME"], "MJDs OK?")
            ch["EXPTIME"] = (texp, "Exposure time (secs)")

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
            nchunk = 20 * win.nx * win.ny

            # allwins contains data of all windows 1 or 2
            allwins = frame[npixel : npixel + nchunk]

            # re-format the data as a 4D array indexed by (ccd,window,y,x)

            # Get number of samples per pixel, with a default of 4.  Pixel
            # data order depends on the number of samples, and the data prior
            # to implementing sampling (i.e. October 2017 WHT run) is the same
            # as nsamps = 4. Use 'as_strided' to perform the de-multiplexing.
            nsamps = self.header.get("ESO DET NSAMP", 4)
            if nsamps == 4:
                data = as_strided(
                    allwins,
                    strides=(8, 2, 40 * win.nx, 40),
                    shape=(5, 4, win.ny, win.nx),
                )
            else:
                # with 1 sample per pixel the data cannot be simply
                # re-viewed but a copy must be made.
                data = as_strided(
                    allwins, strides=(32, 2, 160, 8), shape=(5, 4, len(frame) // 80, 4)
                ).reshape(5, 4, win.ny, win.nx)

            # now build the Windows. This is where we chop off any pre- and
            # over-scan

            for nccd, cnam in enumerate(CNAMS):
                for nquad, qnam in enumerate(QNAMS):
                    wnam = "{:s}{:d}".format(qnam, nwin + 1)
                    if nccd in REFLECTED:
                        # reflections for g (1) and z (4)
                        qnam = QNAMS_REFLECT[qnam]

                    # recover the window and flip parameters
                    win, flip_axes = self.windows[nwin][nccd][nquad]

                    # get the window data and apply flips.
                    windata = data[nccd, nquad]
                    for ax in flip_axes:
                        windata = np.flip(windata, ax)

                    # at this point the data array has the appropriate
                    # orientation. Now finally store in a Window or
                    # two Windows, splitting off the pre-scan section
                    # in the latter case. We also convert to 32-bit
                    # floats to avoid problems down the line.

                    if self.pscan:
                        # pre-scan present. Reduce the size of the data window
                        # in preparation for removal of prescan data.
                        winh = Winhead(
                            win.llx,
                            win.lly,
                            win.nx - self.npscan,
                            win.ny,
                            win.xbin,
                            win.ybin,
                            win.outamp,
                            win,
                        )

                        if nwin == 0 and nquad == 0:
                            # store CCD header in the first Winhead
                            winh.update(cheads[cnam])

                        # Generate name for the prescan Window
                        wpnam = "{:s}P".format(wnam)

                        # work out lly-value of bottom of pre-scan windows
                        # which needs adjusting if the overscan is also on
                        if self.oscan and (qnam == "E" or qnam == "F"):
                            plly = 1 - HCM_NOSCAN
                        else:
                            plly = win.lly

                        if qnam == "E" or qnam == "H":
                            # Prescans are on the left of quadrants E and H
                            # (outputs are on the left). Create the Window
                            # with the image data, stripping off the prescan
                            ccds[cnam][wnam] = Window(
                                winh, windata[:, self.npscan :].astype(np.float32), True
                            )

                            # Store the prescan itself
                            winp = Winhead(
                                1 - HCM_NPSCAN,
                                plly,
                                self.npscan,
                                win.ny,
                                win.xbin,
                                win.ybin,
                                win.outamp,
                                win,
                            )

                            # Create the Window with the pre-scan
                            ccds[cnam][wpnam] = Window(
                                winp, windata[:, : self.npscan].astype(np.float32), True
                            )

                        else:
                            # Prescan on the right. Create the Window with the
                            # image data, stripping off the prescan
                            ccds[cnam][wnam] = Window(
                                winh,
                                windata[:, : -self.npscan].astype(np.float32),
                                True,
                            )

                            # Now save the prescan itself
                            winp = Winhead(
                                HCM_NXTOT + 1,
                                plly,
                                self.npscan,
                                win.ny,
                                win.xbin,
                                win.ybin,
                                win.outamp,
                                win,
                            )

                            # Create the Window with the pre-scan
                            ccds[cnam][wpnam] = Window(
                                winp,
                                windata[:, -self.npscan :].astype(np.float32),
                                True,
                            )

                    else:
                        # No prescan: just store the Window
                        if nwin == 0 and nquad == 0:
                            # store CCD header in the first Winhead
                            win.update(cheads[cnam])

                        ccds[cnam][wnam] = Window(win, windata.astype(np.float32), True)

                    if self.oscan:
                        # over-scan present. Leave any stripped pre-scans
                        # untouched, but remove rows off data windows.

                        win = ccds[cnam][wnam]

                        # Generate name for the overscan Window
                        wonam = "{:s}O".format(wnam)

                        if qnam == "E" or qnam == "F":
                            # Overscans are at the top of quadrants E and F
                            # (outputs are on the bottom).
                            wino = Winhead(
                                win.llx,
                                1 - HCM_NOSCAN,
                                win.nx,
                                self.noscan,
                                win.xbin,
                                win.ybin,
                                win.outamp,
                                win,
                            )

                            # Create the Window with the over-scan from top
                            # part of data frame
                            ccds[cnam][wonam] = Window(
                                wino, win.data[-self.noscan :, :], True
                            )

                            # Chop overscan off data array
                            win.data = win.data[: -self.noscan, :]

                        else:
                            # Overscans at the bottom of quadrants G and H
                            # (outputs are on the top).
                            wino = Winhead(
                                win.llx,
                                HCM_NYTOT + 1,
                                win.nx,
                                self.noscan,
                                win.xbin,
                                win.ybin,
                                win.outamp,
                                win,
                            )

                            # Create the Window with the over-scan
                            ccds[cnam][wonam] = Window(
                                wino, win.data[: self.noscan, :], True
                            )

                            # Chop overscan off data array
                            win.data = win.data[self.noscan :, :]

            # move pointer on for next set of windows
            npixel += nchunk

        # create the MCCD
        mccd = MCCD(ccds, thead)

        # update the frame counter for the next call
        self.nframe += 1

        # show that we have read something
        self.first = False

        # at last, return with the MCCD
        return mccd

    def seek_frame(self, n):
        """
        Moves pointer position to the start of frame n and
        updates the 'nframe' attribute appropriately.

        This is not implemented for the server.
        """
        if self.server:
            raise NotImplementedError("no seek_frame in the case of server access")
        else:
            # move pointer to the start of the last complete frame
            # and adjust nframe
            self._ffile.seek(self._hbytes + self._framesize * (n - 1))
            self.nframe = n

    def seek_last(self):
        """
        In the case of a local file, this sets the pointer position
        to the start of the final complete frame, and updates the 'nframe'
        attribute appropriately. This is not implemented for the server.
        It returns the number of the frame it moves to
        """
        if self.server:
            raise NotImplementedError("no seek_last in the case of server access")
        else:
            # compute byte offset at this point and thus the total number of
            # complete frames
            ntot = self.ntotal()

            # move pointer to this file
            self.seek_frame(ntot)

        return ntot


class Rtbytes(Rhead):
    """Callable, iterable object to returns timing bytes from HiPERCAM raw
    data files.

    This is similar to Rdata except it only reads the timing data
    and does no processing at all, just returning the bytes.

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

        # This sets whether on the next call we should regard the file
        # position as requiring a reset. Used by 'set' to flag up different
        # actions in __call__
        self.reset = False

        if not server:
            self._ffile.seek(
                self._hbytes + self._framesize * self.nframe - self.ntbytes
            )

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self.__call__()
        except HendError:
            raise StopIteration

    def __call__(self, nframe=None):
        """Reads the timing data of one frame from the run the :class:`Rtime`
        is attached to. If `nframe` is None, then it will read the frame it is
        positioned at. If nframe is an integer > 0, it will try to read that
        particular frame; if nframe == 0, it reads the last complete frame.
        nframe == 1 gives the first frame.

        Arguments::

           nframe : int | None
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted.

        Returns:: the timing bytes of the frame
        """

        self.set(nframe)

        if self.server:
            # define command to send to server
            if self.reset:
                # requires a seek as well as a read
                data = json.dumps(dict(action="get_frame", frame_number=self.nframe))
            else:
                # just read next set of bytes
                data = json.dumps(dict(action="get_next"))

            # get data
            self._ws.send(data)
            raw_bytes = self._ws.recv()

            if len(raw_bytes) == 0:
                # if we are trying to access a file that has not yet been
                # written, 0 bytes will be returned. In this case we return
                # with None. NB we do not update self.nframe in this case.
                return None

            # extract the timing data
            tbytes = raw_bytes[-self.ntbytes :]

        elif self.nframe > self.ntotal():
            # We have attempted to access a non-existent frame
            raise HendError("no more frames to access")

        else:
            # find the start of the timing data if necessary
            if self.reset:
                self._ffile.seek(
                    self._hbytes + self._framesize * self.nframe - self.ntbytes
                )

            # read the timing data, skip forward to next set
            tbytes = self._ffile.read(self.ntbytes)
            if len(tbytes) != self.ntbytes:
                raise HendError("failed to read timing bytes")
            self._ffile.seek(self._framesize - self.ntbytes, 1)

        # update the internal frame counter
        self.nframe += 1

        # set the reset status to save effort if next time
        # we simply want the next frame
        self.reset = False

        return tbytes

    def set(self, nframe=1):
        """Resets the position so that we are just about to read
        the timing data of nframe, but does not read them unlike
        __call__. Set the internal "reset" flag if needed.

        Arguments::

           nframe : int | None
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates just read whatever is the next frame.

        Returns:: the frame number of entry
        """

        old_frame = self.nframe

        if nframe is not None:
            if nframe == 0:
                # go for the last one
                nframe = self.ntotal()

            # record whether we need to be in reset mode
            # update frame counter
            self.reset = self.nframe != nframe
            self.nframe = nframe

        return old_frame


class Rtime(Rtbytes):
    """Callable, iterable object to generate timing data only from HiPERCAM raw data files.

    This is the same as Rtbytes except it interprets the bytes.
    """

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

        Returns:: (tstamp, tinfo, tflag) or None if you try to read a non-existent
        frame (useful if you are waiting from frames to be added). Here
        'tstamp' is an astropy.time.Time equivalent to the GPS timestamp
        associated with the frame, with no correction, while 'tinfo' is a
        5-elements tuple, one for each CCD listing for each one the MJD at
        mid-exposure, an exposure time in seconds and a flag to indicate
        whether the time is thought reliable, which is only the case if
        real data comes with the frame. The final tflag reflects whether
        there were any satellites and the GPS had been synced
        """

        tbytes = super().__call__(nframe)

        if tbytes is None:
            return None

        # Interpret the timing bytes
        (
            frameCount,
            timeStampCount,
            years,
            day_of_year,
            hours,
            mins,
            seconds,
            nanoseconds,
            nsats,
            synced,
        ) = htimer(tbytes)

        tflag = nsats != -1 or synced != -1
        try:
            time_string = "{}:{}:{}:{}:{:.7f}".format(
                years, day_of_year, hours, mins, seconds + nanoseconds / 1e9
            )
            tstamp = Time(time_string, format="yday", precision=7)
        except ValueError:
            tflag = False
            # invalid time; pretend we are on 2000-01-01 taking one frame per second.
            tstamp = Time(
                51544 + (self.nframe - 1) / 86400.0, format="mjd", precision=7
            )

        tinfo = []
        for nccd in range(5):
            tmid, texp, tiflag = self.timing(self.nframe - 1, nccd)
            tinfo.append((tstamp.mjd + tmid, texp, tiflag))

        # Return timing data
        return (tstamp, tuple(tinfo), tflag)


def htimer(tbytes):
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
    buf = struct.pack(
        "<HHHHHHHHHHHHHHHHH",
        *(val + BZERO for val in struct.unpack(">hhhhhhhhhhhhhhhhh", tbytes[-36:-2])),
    )
    return struct.unpack("<IIIIIIIIbb", buf)


class HendError(HipercamError):
    """
    Exception for the standard way to reach the end of a HiPERCAM raw data
    file (attempt to read out-of-range frame). This allows the iterator to die
    silently in this case while raising exceptions for less anticipated cases.
    """

    pass
