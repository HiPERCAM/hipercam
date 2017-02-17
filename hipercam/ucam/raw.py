"""
Section for handling the raw ULTRACAM data and xml files
"""

import struct
import warnings
import xml.dom.minidom
import numpy as np

from .constants import *
from .server import get_nframe_from_server, URL
from .time import Time
from .uhead import Uhead

class Rwin:
    """Trivial container class for basic window info to help readability of
    code. Sets attributes llx, lly, nx, ly representing the lower-left pixel
    of a window and its binned dimensions.  Used in the 'wins' attribute of
    Rhead

    """
    def __init__(self, llx, lly, nx, ny):
        self.llx = llx
        self.lly = lly
        self.nx  = nx
        self.ny  = ny

    def __str__(self):
        return str(self.llx) + ',' + str(self.lly) + ',' + str(self.nx) + ',' + str(self.ny)

class Rhead:
    """Represents essential header info of Ultracam/Ultraspec data read from a
    run###.xml file.

    Attributes set::

       run : (string)
          The ULTRACAM run e.g. 'run005'

       server : (bool)
          True/False to indicate server access

       application : (string)
          Data acquisition application template name.

       framesize : (int)
          Total number of bytes per frame.

       headerwords : (int)
          Number of words (2-bytes/word) in timing info at start of a frame.

       instrument : (string)
          Instrument name.

       nccd : (int)
          Number of CCDs

       mode : (string)
          A standardised summary of the readout mode derived from the application name.

       user : (dict)
          Dictionary of user information. Set to None if there was none found.

       win : (list)
          A list of Rwin objects, one per window. ULTRACAM data is multi-CCD
          but the windows of each CCD are identical so the information is only stored
          once for all CCDs.

       xbin : (int)
          binning factor in X direction

       ybin : (int)
          binning factor in Y direction

       nxmax : (int)
          maximum X dimension, unbinned pixels

       nymax : (int)
          maximum Y dimension, unbinned pixels.

       exposeTime : (float)
          exposure delay, secs

       numexp : (int)
          integer exposure number; can be -1 meaning unlimited

       gainSpeed :
          readout speed (ULTRACAM)

       v_ft_clk :
          vertical frame transfer parameter

       nblue : (int)
          how often to read the blue CCD (ULTRACAM)

       speed :
          readout speed (ULTRASPEC)

       en_clr :
          enable clear flag (ULTRASPEC)

       hv_gain :
          HV gain (ULTRASPEC)

       output :
          which output (ULTRASPEC)

       version :
          version number

       target :
          target name (None if not found)

       filters :
          filter names (None if not found)

    """

    def __init__(self, run, server=False):
        """Reads a run###.xml file. UltracamErrors are thrown if some items are not
        found.  In some case it will carry on and corresponding attributes are
        returned as None.

        Arguments::

           run : (string)
              run name e.g. 'run002'. Can include path to disk file.

           server : (bool)
              True to attempt to access the ATC FileServer. It uses a URL set in an
              environment variable "ULTRACAM_DEFAULT_URL" in this instance.

        """

        self.run = run
        self.server = server
        if server:
            # get from server
            full_url = URL + run + '?action=get_xml'
            sxml = urllib.request.urlopen(full_url).read()
            udom = xml.dom.minidom.parseString(sxml)
        else:
            # local disk file
            udom = xml.dom.minidom.parse(run + '.xml')

        # Find framesize and headerwords.
        node = udom.getElementsByTagName('data_status')[0]
        self.framesize = int(node.getAttribute('framesize'))
        self.headerwords = int(node.getElementsByTagName('header_status')[0].getAttribute(
            'headerwords'))

        # Frame format and other detail.
        node = udom.getElementsByTagName('instrument_status')[0]
        self.instrument = node.getElementsByTagName('name')[0].childNodes[0].data
        if self.instrument == 'Ultracam':
            self.instrument = 'ULTRACAM'
            self.nccd = 3
            self.nxmax, self.nymax = 1080, 1032
        elif self.instrument == 'Ultraspec':
            self.instrument = 'ULTRASPEC'
            self.nccd = 1
            self.nxmax, self.nymax = 1056, 1072
        else:
            raise ValueError(
                'Rhead.__init__: run = {:s}, failed to identify instrument.'.format(self.run)
            )

        self.application = [nd for nd in node.getElementsByTagName('application_status') \
                            if nd.getAttribute('id') == 'SDSU Exec'][0].getAttribute('name')

        # gather together majority of values
        param = {}
        for nd in node.getElementsByTagName('parameter_status'):
            param[nd.getAttribute('name')] = nd.getAttribute('value')

        # get user info, if present
        try:
            nlist = udom.getElementsByTagName('user')
            if len(nlist):
                user = {}
                node = nlist[0]
                for nd in node.childNodes:
                    if nd.nodeType == xml.dom.Node.ELEMENT_NODE and nd.hasChildNodes():
                        user[nd.tagName] = nd.childNodes[0].data
            else:
                user = None
        except Exception as err:
            user = None

        # Translate applications into meaningful mode names
        app = self.application
        if app == 'ap8_250_driftscan' or app == 'ap8_driftscan' or app == 'ap_drift_bin2' or \
           app == 'appl8_driftscan_cfg':
            self.mode    = 'DRIFT'
        elif app == 'ap5_250_window1pair' or app == 'ap5_window1pair' or app == 'ap_win2_bin8' or \
             app == 'ap_win2_bin2' or app == 'appl5_window1pair_cfg':
            self.mode    = '1-PAIR'
        elif app == 'ap5b_250_window1pair' or app == 'appl5b_window1pair_cfg':
            self.mode    = '1-PCLR'
        elif app == 'ap6_250_window2pair' or app == 'ap6_window2pair' or \
             app == 'ap_win4_bin1' or app == 'ap_win4_bin8' or app == 'appl6_window2pair_cfg':
            self.mode    = '2-PAIR'
        elif app == 'ap7_250_window3pair' or app == 'ap7_window3pair' or \
             app == 'appl7_window3pair_cfg':
            self.mode    = '3-PAIR'
        elif app == 'ap3_250_fullframe' or app == 'ap3_fullframe' or \
             app == 'appl3_fullframe_cfg':
            self.mode    = 'FFCLR'
        elif app == 'appl4_frameover_cfg' or app == 'ap4_frameover':
            self.mode    = 'FFOVER'
        elif app == 'appl10_frameover_mindead_cfg':
            self.mode    = 'FFOVNC'
        elif app == 'ap9_250_fullframe_mindead' or app == 'ap9_fullframe_mindead' or \
             app == 'appl9_fullframe_mindead_cfg':
            self.mode    = 'FFNCLR'
        elif app == 'ccd201_winbin_con' or app == 'ccd201_winbin_cfg':
            if int(param['X2_SIZE']) == 0:
                self.mode    = 'USPEC-1'
            elif int(param['X3_SIZE']) == 0:
                self.mode    = 'USPEC-2'
            elif int(param['X4_SIZE']) == 0:
                self.mode    = 'USPEC-3'
            else:
                self.mode    = 'USPEC-4'
        elif app == 'ccd201_driftscan_cfg':
            self.mode    = 'UDRIFT'
        elif app == 'ap1_poweron' or app == 'ap1_250_poweron' or \
             app == 'ap2_250_poweroff' or app == 'appl1_pon_cfg' or \
             app == 'appl2_pof_cfg' or app == 'ccd201_pon_cfg':
            self.mode = 'PONOFF'
            return
        else:
            raise ValueError(
                'Rhead.__init__: file = {:s} failed to identify application = {:s}'.format(
                    self.run,app)
            )

        # binning factors
        self.xbin = int(param['X_BIN_FAC']) if 'X_BIN_FAC' in param \
                    else int(param['X_BIN'])
        self.ybin = int(param['Y_BIN_FAC']) if 'Y_BIN_FAC' in param \
                    else int(param['Y_BIN'])

        # Windows. For each one store: x & y coords of lower-left pixel,
        # binned dimensions
        self.win = []
        fsize = 2*self.headerwords
        if self.instrument == 'ULTRACAM':
            try:
                self.exposeTime = float(param['EXPOSE_TIME'])
            except ValueError as err:
                raise ValueError(
                    'Rhead.__init__: file = {:s} failed to interpret EXPOSE_TIME'.format(self.run)
                )

            self.numexp = int(param['NO_EXPOSURES'])
            self.gainSpeed = hex(int(param['GAIN_SPEED']))[2:] \
                              if 'GAIN_SPEED' in param else None

            if 'V_FT_CLK' in param:
                self.v_ft_clk = six.indexbytes(struct.pack('I',int(param['V_FT_CLK'])),2)
            elif app == 'appl7_window3pair_cfg':
                self.v_ft_clk = 140;
            else:
                self.v_ft_clk = 0

            self.nblue = int(param['NBLUE']) if 'NBLUE' in param else 1

            if self.mode == 'FFCLR' or self.mode == 'FFNCLR':
                self.win.append(Rwin(  1, 1, 512//self.xbin, 1024//self.ybin))
                self.win.append(Rwin(513, 1, 512//self.xbin, 1024//self.ybin))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

            elif self.mode == 'FFOVER' or self.mode == 'FFOVNC':
                # In overscan mode the extra pixels are clocked out after the
                # data has been read which effectively means they should
                # appear in the centre of the chip.  However, this would ruin
                # the location of the physical pixels relative to all other
                # formats (unless they always included a centre gap). Thus the
                # extra sections are placed off to the right-hand and top
                # sides where they do not affect the pixel registration. This
                # code requires some corresponding jiggery-pokery in Rdata
                # because the actual data comes in in just two windows. The 6
                # windows come in 3 pairs of equal sizes, hence the single
                # fsize increment line per pair.

                # first set the physical data windows
                self.win.append(Rwin(  1, 1, 512//self.xbin, 1024//self.ybin))
                self.win.append(Rwin(513, 1, 512//self.xbin, 1024//self.ybin))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

                # left overscan (place on right)
                self.win.append(Rwin(1025, 1, 28//self.xbin, 1032//self.ybin))

                # right overscan
                self.win.append(Rwin(1053, 1, 28//self.xbin, 1032//self.ybin))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

                # top left overscan
                self.win.append(Rwin(1, 1025, 512//self.xbin, 8//self.ybin))

                # top right overscan
                self.win.append(Rwin(513, 1025, 512//self.xbin, 8//self.ybin))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

            else:

                ystart = int(param['Y1_START'])
                xleft = int(param['X1L_START'])
                xright = int(param['X1R_START'])
                nx = int(param['X1_SIZE']) // self.xbin
                ny = int(param['Y1_SIZE']) // self.ybin
                self.win.append(Rwin(xleft, ystart, nx, ny))
                self.win.append(Rwin(xright, ystart, nx, ny))

                fsize += 12*self.win[-1].nx*self.win[-1].ny

            if self.mode == '2-PAIR' or self.mode == '3-PAIR':
                ystart = int(param['Y2_START'])
                xleft = int(param['X2L_START'])
                xright = int(param['X2R_START'])
                nx = int(param['X2_SIZE']) // self.xbin
                ny = int(param['Y2_SIZE']) // self.ybin
                self.win.append(Rwin(xleft, ystart, nx, ny))
                self.win.append(Rwin(xright, ystart, nx, ny))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

            if self.mode == '3-PAIR':
                ystart = int(param['Y3_START'])
                xleft = int(param['X3L_START'])
                xright = int(param['X3R_START'])
                nx = int(param['X3_SIZE']) // self.xbin
                ny = int(param['Y3_SIZE']) // self.ybin
                self.win.append(Rwin(xleft,ystart,nx,ny))
                self.win.append(Rwin(xright,ystart,nx,ny))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

        elif self.instrument == 'ULTRASPEC':

            self.exposeTime = float(param['DWELL'])
            self.numexp = int(param['NUM_EXPS'])
            self.speed = ('S' if param['SPEED'] == '0' else \
                          ('M' if param['SPEED'] == '1' else 'F')) \
                if 'SPEED' in param else None
            self.en_clr = (True if param['EN_CLR'] == '1' else False) \
                          if 'EN_CLR' in param else None
            self.hv_gain = int(param['HV_GAIN']) if 'HV_GAIN' \
                           in param else None
            self.output = ('N' if param['OUTPUT'] == '0' else 'A') \
                          if 'OUTPUT' in param else None

            xstart = int(param['X1_START'])
            ystart = int(param['Y1_START'])
            nx = int(param['X1_SIZE'])
            ny = int(param['Y1_SIZE'])
            self.win.append(Rwin(xstart,ystart,nx,ny))
            fsize += 2*self.win[-1].nx*self.win[-1].ny

            if self.mode == 'USPEC-2' or self.mode == 'USPEC-3' or \
               self.mode == 'USPEC-4' or self.mode == 'UDRIFT':
                xstart = int(param['X2_START'])
                ystart = ystart if self.mode == 'UDRIFT' else \
                         int(param['Y2_START'])
                nx = int(param['X2_SIZE'])
                ny = ny if self.mode == 'UDRIFT' else int(param['Y2_SIZE'])
                self.win.append(Rwin(xstart,ystart,nx,ny))
                fsize += 2*self.win[-1].nx*self.win[-1].ny

            if self.mode == 'USPEC-3' or self.mode == 'USPEC-4':
                xstart = int(param['X3_START'])
                ystart = int(param['Y3_START'])
                nx = int(param['X3_SIZE'])
                ny = int(param['Y3_SIZE'])
                self.win.append(Rwin(xstart,ystart,nx,ny))
                fsize += 2*self.win[-1].nx*self.win[-1].ny

            if self.mode == 'USPEC-4':
                xstart = int(param['X4_START'])
                ystart = int(param['Y4_START'])
                nx = int(param['X4_SIZE'])
                ny = int(param['Y4_SIZE'])
                self.win.append(Rwin(xstart,ystart,nx,ny))
                fsize += 2*self.win[-1].nx*self.win[-1].ny

        if fsize != self.framesize:
            raise ValueError(
                'Rhead.__init__: file = {:s}. Framesize = {:s} clashes with calculated value = {:s}'.format(self.run,self.framesize,fsize)
            )

        # nasty stuff coming up ...
        self.version   = int(user['revision']) if user is not None and \
                         'revision' in user else \
                         (int(param['REVISION']) if 'REVISION' in param \
                          else int(param['VERSION']) if 'VERSION' in param \
                          else -1)

        if 'REVISION' in param or 'VERSION' in param:
            vcheck = int(param['REVISION']) if 'REVISION' in param else \
                     int(param['VERSION'])
            if vcheck != self.version:
                raise ValueError(
                    'Rhead.__init__: clashing version numbers: {:s} vs {:s}'.format(
                        self.version,vcheck)
                )

        if self.headerwords == 16:
            VERSIONS = [100222, 111205, 120716, 120813, 130307, 130317, 140331]
            if self.version not in VERSIONS:
                raise ValueError(
                    'Rhead.__init__: could not recognise version = {:s}'.format(self.version)
                )

        self.whichRun = ''
        if self.instrument == 'ULTRACAM':
            if user is None:
                self.timeUnits = 0.001
            else:
                self.timeUnits = 0.0001
            if 'REVISION' not in param and 'VERSION' not in param and 'V_FT_CLK' not in param:
                self.whichRun = 'MAY2002'
        else:
            if user is not None and self.headerwords == 16 and  self.version >= 120813:
                self.timeUnits = 0.0001
            else:
                self.timeUnits = 0.001

        # convert to seconds
        self.exposeTime *= self.timeUnits

        # conditional loading of headers
        if user and 'target' in user: self.target = user['target']
        if user and 'filters' in user: self.filters = user['filters']
        if user and 'PI' in user: self.pi = user['PI']
        if user and 'ID' in user: self.id = user['ID']
        if user and 'Observers' in user: self.observers = user['Observers']
        if user and 'flags' in user: self.dtype = user['flags']
        if user and 'ccd_temp' in user: self.ccdtemp = user['ccd_temp']
        if user and 'SlidePos' in user: self.slidepos = user['SlidePos'].split()[0]
        if user and 'RA' in user: self.RA  = user['RA']
        if user and 'Dec' in user: self.Dec = user['Dec']
        if user and 'Tracking' in user: self.track = user['Tracking']
        if user and 'TTflag' in user: self.ttflag = user['TTflag']
        if user and 'Focus' in user: self.focus = user['Focus']
        if user and 'PA' in user: self.PA = user['PA']
        if user and 'Eng_PA' in user: self.engpa = user['Eng_PA']
        if user and 'ccd_temp' in user: self.ccdtemp = user['ccd_temp']
        if user and 'finger_temp' in user: self.fingertemp = user['finger_temp']
        if user and 'finger_pcent' in user: self.fingerpcent = user['finger_pcent']

    def npix(self):
        """
        Returns number of (binned) pixels per CCD
        """
        np = 0
        for win in self.win:
            np += win.nx*win.ny
        return np

    def isPonoff(self):
        """
        Is the run a power on / off? (no data)
        """
        return self.mode == 'PONOFF'

class Ahead(Uhead):
    """
    Sub-class of Uhead to allow checked addition by attribute
    from an Rhead into a Uhead.
    """
    def __init__(self, rhead):
        Uhead.__init__(self)
        self.rhead = rhead

    def add_attr(self, key, attr, itype, comment):
        """
        attr only added in if it exists in self.rhead
        """
        if hasattr(self.rhead, attr):
            value = getattr(self.rhead, attr)
            Uhead.add_entry(self, key, value, itype, comment)

class Rdata (Rhead):
    """
    Callable, iterable object to represent ULTRACAM/SPEC raw data files.

    The idea is to open the file, and then the object generated can be
    used to deliver frames by specifying a frame number. Frames can be
    read individually e.g.::

      rdat = Rdata('run045')
      fr10 = rdat(10)
      fr11 = rdat()

    reads frame numbers 10 and then 11 in from 'run045', or sequentially::

      for frm in Rdata('run045'):
         print 'nccd = ',frm.nccd()

    Rdata maintains an internal file object that is always at the start of a
    frame. This enables sequential reads to be swift. If an attempt is made
    to access a frame that does not exist, it defaults to the start of the
    file.

    The above code returns :class:`trm.ultracam.MCCD` objects for ULTRACAM data. By default
    the frames are always read as :class:`trm.ultracam.MCCD` objects, however there is
    always a flag

    :class:`trm.ultracam.Rdata` can talk to the ATC FileServer if it is running, e.g.::

      rdat = Rdata('run045',server=True)

    """

    def __init__(self, run, nframe=1, flt=True, server=False, ccd=False):
        """
        Connects to a raw data file for reading. The file is kept open.
        The file pointer is set to the start of frame nframe. The Rdata
        object can then generate MCCD or CCD objects through being called
        as a function or iterator.

        Args:
          run (string) : run name, e.g. 'run036'.

          nframe (int) : frame number for first read, starting at 1 as the first.

          flt (bool ) : True for reading data in as floats. This is the default for
               safety, however the data are stored on disk as unsigned 2-byte
               ints. If you are not doing much to the data, and wish to keep them
               in this form for speed and efficiency, then set flt=False.  This
               parameter is used when iterating through an Rdata. The __call__
               method can override it.

          server (bool) : True/False for server vs local disk access

          ccd (bool) : flag to read data as a :class:`trm.ultracam.CCD` rather
                       than a :class:`trm.ultracam.MCCD` object if only one CCD
                       per frame. Default is always to read as an MCCD.
        """

        Rhead.__init__(self, run, server)
        if self.isPonoff():
            raise PowerOnOffError('Rdata.__init__: attempted to read a power on/off')

        # Attributes set are:
        #
        # _fobj   -- file object opened on data file (set to None if using server)
        # _nf     -- next frame to be read
        # _run    -- name of run
        # _flt    -- whether to read as float (else uint16)
        # _tstamp -- list of immediately preceding times
        if server:
            self._fobj   = None
        else:
            if six.PY3:
                """
                This exists because in Python 3, `open()` returns an
                `io.BufferedReader` by default.  This is bad, because
                `io.BufferedReader` doesn't support random access, which we may need in
                some cases.  In the Python 3 case (implemented in the py3compat module)
                we must call open with buffering=0 to get a raw random-access file
                """
                self._fobj   = open(run + '.dat', 'rb', buffering=0)
            else:
                self._fobj   = open(run + '.dat', 'rb')

        self._nf     = nframe
        self._run    = run
        self._flt    = flt
        self._tstamp = []
        self._ccd    = ccd
        if not server and nframe != 1:
            self._fobj.seek(self.framesize*(nframe-1))

    def __iter__(self):
        """
        Generator to allow Rdata to function as an iterator.
        This produces the same type of object as __call__ does.
        """
        try:
            while 1:
                yield self.__call__(flt=self._flt)
        except UendError:
            pass
        except urllib.error.HTTPError:
            pass

    def set(self, nframe=1):
        """
        Sets the internal file pointer to point at frame nframe.

        Args:
          nframe (int) : frame number to get, starting at 1. 0 for the
                         last (complete) frame. A value of 'None' will be
                         ignored. A value < 0 will cause an exception. A
                         value greater than the number of frames in the
                         file will work, but will cause an exception to be
                         raised on the next attempted read.
        """

        # position read pointer
        if nframe is not None:
            if nframe < 0:
                raise UltracamError('Rdata.set: nframe < 0')
            elif nframe == 0:
                if self.server:
                    self._nf = get_nframe_from_server(self.run)
                else:
                    # wind to end of file, work out where we are
                    # to deduce number of frames
                    self._fobj.seek(0,2)
                    fp = self._fobj.tell()
                    nf = fp // self.framesize
                    self._fobj.seek(self.framesize*(nf-1)-fp,2)
                    self._nf = nf

            elif self._nf != nframe:
                if not self.server:
                    self._fobj.seek(self.framesize*(nframe-1))
                self._nf = nframe

    def __call__(self, nframe=None, flt=None):
        """
        Reads the data of frame nframe (starts from 1) and returns a
        corresponding CCD or UCAM object, depending upon the type of data. If
        nframe is None, just reads whatever frame we are on. Raises an
        exception if it fails to read data.  Resets to start of the file in
        this case. The data are stored internally as either 4-byte floats or
        2-byte unsigned ints.

        nframe -- frame number to get, starting at 1. 0 for the last
                  (complete) frame. None just returns the next frame.

        flt    -- Set True to reading data in as floats. The data are stored on
                  disk as unsigned 2-byte ints. If you are not doing much to
                  the data, and wish to keep them in this form for speed and
                  efficiency, then set flt=False. If None then the value used when
                  constructing the MCCD will be used.

        Returns a UCAM object for ULTRACAM, CCD for ULTRASPEC.
        """

        if flt is None: flt = self._flt

        # position read pointer
        self.set(nframe)

        if self.server:
            # read timing and data in one go from the server
            full_url = URL + self.run + '?action=get_frame&frame=' + str(self._nf-1)
            buff     = urllib.request.urlopen(full_url).read()
            if len(buff) != self.framesize:
                self._nf = 1
                raise UltracamError('Rdata.__call__: failed to read frame ' + str(self._nf) +
                                    ' from FileServer. Buffer length vs expected = '
                                    + str(len(buff)) + ' vs ' + str(self.framesize) + ' bytes.')

            # have data. Re-format into the timing bytes and unsigned 2 byte
            # int data buffer
            tbytes = buff[:2*self.headerwords]
            buff   = np.fromstring(buff[2*self.headerwords:],dtype='uint16')
        else:
            # read timing bytes
            tbytes = self._fobj.read(2*self.headerwords)
            if len(tbytes) != 2*self.headerwords:
                self._fobj.seek(0)
                self._nf = 1
                raise UendError('Rdata.__call__: failed to read timing bytes')

            # read data
            buff = np.fromfile(self._fobj,'<u2',int(self.framesize/2-self.headerwords))
            if len(buff) != self.framesize/2-self.headerwords:
                self._fobj.seek(0)
                self._nf = 1
                raise UltracamError('Rdata.__call__: failed to read frame ' + str(self._nf) +
                                    '. Buffer length vs attempted = '
                                    + str(len(buff)) + ' vs ' + str(self.framesize/2-self.headerwords))

        # OK from this point, both server and local disk methods are the same
        if self.instrument == 'ULTRACAM':
            time,info,blueTime,badBlue = utimer(tbytes, self, self._nf)
        elif self.instrument == 'ULTRASPEC':
            time,info = utimer(tbytes, self, self._nf)

        # move frame counter on by one
        self._nf += 1

        # build header
        head = Ahead(self)

        head.add_entry('User','Data entered by user at telescope')
        head.add_attr('User.target', 'target', ITYPE_STRING, 'Object name')
        head.add_attr('User.pi','pi',ITYPE_STRING,'Principal investigator')
        head.add_attr('User.id','id',ITYPE_STRING,'Programme ID')
        head.add_attr('User.observers','observers',ITYPE_STRING,'Observers')
        head.add_attr('User.dtype','dtype',ITYPE_STRING,'Data type')

        head.add_entry('Instrument','Instrument setup information')
        head.add_attr('Instrument.instrument','instrument',ITYPE_STRING,
                      'Instrument identifier')
        head.add_attr('Instrument.headerwords','headerwords',ITYPE_INT,
                      'Number of 2-byte words in timing')
        head.add_attr('Instrument.framesize','framesize',ITYPE_INT,
                       'Total number of bytes per frame')

        head.add_entry('Run', 'Run specific information')
        head.add_attr('Run.run','_run',ITYPE_STRING,'run the frame came from')
        head.add_attr('Run.mode','mode',ITYPE_STRING,'readout mode used')
        head.add_entry('Run.ntmin',info['ntmin'],ITYPE_INT,
                       'number of sequential timestamps needed')
        head.add_attr('Run.filters','filters',ITYPE_STRING,
                      'Filter name or names')
        head.add_attr('Run.expose','exposeTime',ITYPE_FLOAT,'exposure time')
        if self.instrument == 'ULTRASPEC':
            head.add_attr('Run.output','output',ITYPE_STRING,'CCD output used')
            head.add_attr('Run.speed','speed',ITYPE_STRING,'Readout speed')
        elif self.instrument == 'ULTRACAM':
            head.add_attr('Run.speed','gainSpeed',ITYPE_STRING,'Readout speed')

        head.add_attr('Run.focus','focus',ITYPE_FLOAT,
                      'Telescope focus')
        head.add_attr('Run.ccdtemp','ccdtemp',ITYPE_FLOAT,
                      'CCD temperature (K)')
        head.add_attr('Run.fingtemp','fingertemp',ITYPE_FLOAT,
                      'Cold finger temperature (K)')
        head.add_attr('Run.fingpcen','fingerpcent',ITYPE_FLOAT,
                      'Cold finger percentage')
        head.add_attr('Run.slidepos','slidepos',ITYPE_STRING,
                      'Slide position (pixels)')
        head.add_attr('Run.RA','RA',ITYPE_STRING,'Right Ascension (J2000)')
        head.add_attr('Run.Dec','Dec',ITYPE_STRING,'Declination (J2000)')
        head.add_attr('Run.PA','PA',ITYPE_FLOAT,'Position angle (degrees)')
        head.add_attr('Run.EngPA','engpa',ITYPE_FLOAT,
                      'Engineering position angle (degrees)')
        head.add_attr('Run.track','track',ITYPE_STRING,
                      'Telescope judged to be tracking by usdriver')
        head.add_attr('Run.ttflag','ttflag',ITYPE_STRING,
                      'Telescope judged to be tracking by TCS')

        head.add_entry('Frame', 'Frame specific information')
        head.add_attr('Frame.frame','_nf',ITYPE_INT,'frame number within run')
        head.add_entry('Frame.midnight',info['midnightCorr'],ITYPE_BOOL,
                       'midnight bug correction applied')
        head.add_entry('Frame.ferror',info['frameError'],ITYPE_BOOL,
                       'problem with frame numbers found')

        # interpret data
        xbin, ybin = self.xbin, self.ybin
        if self.instrument == 'ULTRACAM':
            # 3 CCDs. Windows come in pairs. Data from equivalent windows come out
            # on a pitch of 6. Some further jiggery-pokery is involved to get the
            # orientation of the frames correct.
            wins1, wins2, wins3 = [],[],[]

            if self.mode != 'FFOVER' and self.mode != 'FFOVNC':
                # Non-overscan modes:
                # flag indicating that outer pixels will be removed. This is because of a readout bug
                # that affected all data taken prior to the VLT run of May 2007 spotted via the
                # lack of a version number in the xml file
                strip_outer = self.version == -1
                noff = 0
                for wl, wr in zip(self.win[::2],self.win[1::2]):
                    npix = 6*wl.nx*wl.ny
                    if flt:
                        if strip_outer:
                            wins1.append(
                                Window(np.reshape(buff[noff:noff+npix:6].astype(np.float32),
                                                  (wl.ny,wl.nx))[:,1:],wl.llx,wl.lly,xbin,ybin))
                            wins1.append(
                                Window(np.reshape(buff[noff+1:noff+npix:6].astype(np.float32),
                                                  (wr.ny,wr.nx))[:,-2::-1],wr.llx+xbin,wr.lly,xbin,ybin))
                            wins2.append(
                                Window(np.reshape(buff[noff+2:noff+npix:6].astype(np.float32),
                                                  (wl.ny,wl.nx))[:,1:],wl.llx,wl.lly,xbin,ybin))
                            wins2.append(
                                Window(np.reshape(buff[noff+3:noff+npix:6].astype(np.float32),
                                                  (wr.ny,wr.nx))[:,-2::-1],wr.llx+xbin,wr.lly,xbin,ybin))
                            wins3.append(
                                Window(np.reshape(buff[noff+4:noff+npix:6].astype(np.float32),
                                                  (wl.ny,wl.nx))[:,1:],wl.llx,wl.lly,xbin,ybin))
                            wins3.append(
                                Window(np.reshape(buff[noff+5:noff+npix:6].astype(np.float32),
                                                  (wr.ny,wr.nx))[:,-2::-1],wr.llx+xbin,wr.lly,xbin,ybin))
                        else:
                            wins1.append(
                                Window(np.reshape(buff[noff:noff+npix:6].astype(np.float32),
                                                  (wl.ny,wl.nx)),wl.llx,wl.lly,xbin,ybin))
                            wins1.append(
                                Window(np.reshape(buff[noff+1:noff+npix:6].astype(np.float32),
                                                  (wr.ny,wr.nx))[:,::-1],wr.llx,wr.lly,xbin,ybin))
                            wins2.append(
                                Window(np.reshape(buff[noff+2:noff+npix:6].astype(np.float32),
                                                  (wl.ny,wl.nx)),wl.llx,wl.lly,xbin,ybin))
                            wins2.append(
                                Window(np.reshape(buff[noff+3:noff+npix:6].astype(np.float32),
                                                  (wr.ny,wr.nx))[:,::-1],wr.llx,wr.lly,xbin,ybin))
                            wins3.append(
                                Window(np.reshape(buff[noff+4:noff+npix:6].astype(np.float32),
                                                  (wl.ny,wl.nx)),wl.llx,wl.lly,xbin,ybin))
                            wins3.append(
                                Window(np.reshape(buff[noff+5:noff+npix:6].astype(np.float32),
                                                  (wr.ny,wr.nx))[:,::-1],wr.llx,wr.lly,xbin,ybin))
                    else:
                        if strip_outer:
                            wins1.append(Window(np.reshape(buff[noff:noff+npix:6],(wl.ny,wl.nx))[:,1:],
                                                wl.llx,wl.lly,xbin,ybin))
                            wins1.append(Window(np.reshape(buff[noff+1:noff+npix:6],(wr.ny,wr.nx))[:,-2::-1],
                                                wr.llx+xbin,wr.lly,xbin,ybin))
                            wins2.append(Window(np.reshape(buff[noff+2:noff+npix:6],(wl.ny,wl.nx))[:,1:],
                                                wl.llx,wl.lly,xbin,ybin))
                            wins2.append(Window(np.reshape(buff[noff+3:noff+npix:6],(wr.ny,wr.nx))[:,-2::-1],
                                                wr.llx+xbin,wr.lly,xbin,ybin))
                            wins3.append(Window(np.reshape(buff[noff+4:noff+npix:6],(wl.ny,wl.nx))[:,1:],
                                                wl.llx,wl.lly,xbin,ybin))
                            wins3.append(Window(np.reshape(buff[noff+5:noff+npix:6],(wr.ny,wr.nx))[:,-2::-1],
                                                wr.llx+xbin,wr.lly,xbin,ybin))
                        else:
                            wins1.append(Window(np.reshape(buff[noff:noff+npix:6],(wl.ny,wl.nx)),
                                                wl.llx,wl.lly,xbin,ybin))
                            wins1.append(Window(np.reshape(buff[noff+1:noff+npix:6],(wr.ny,wr.nx))[:,::-1],
                                                wr.llx,wr.lly,xbin,ybin))
                            wins2.append(Window(np.reshape(buff[noff+2:noff+npix:6],(wl.ny,wl.nx)),
                                                wl.llx,wl.lly,xbin,ybin))
                            wins2.append(Window(np.reshape(buff[noff+3:noff+npix:6],(wr.ny,wr.nx))[:,::-1],
                                                wr.llx,wr.lly,xbin,ybin))
                            wins3.append(Window(np.reshape(buff[noff+4:noff+npix:6],(wl.ny,wl.nx)),
                                                wl.llx,wl.lly,xbin,ybin))
                            wins3.append(Window(np.reshape(buff[noff+5:noff+npix:6],(wr.ny,wr.nx))[:,::-1],
                                                wr.llx,wr.lly,xbin,ybin))
                    noff += npix
            else:
                # Overscan modes need special re-formatting. See the description under Rhead
                # for more on this. The data come in the form of two windows 540 by 1032
                # (divided by binning factors). The first thing we do is read these windows
                # into 6 numpy arrays.
                nxb  = 540 // xbin
                nyb  = 1032 // ybin
                npix = 6*nxb*nyb
                if flt:
                    winl1 = np.reshape(buff[:npix:6].astype(np.float32),(nyb,nxb))
                    winr1 = np.reshape(buff[1:npix:6].astype(np.float32),(nyb,nxb))[:,::-1]
                    winl2 = np.reshape(buff[2:npix:6].astype(np.float32),(nyb,nxb))
                    winr2 = np.reshape(buff[3:npix:6].astype(np.float32),(nyb,nxb))[:,::-1]
                    winl3 = np.reshape(buff[4:npix:6].astype(np.float32),(nyb,nxb))
                    winr3 = np.reshape(buff[5:npix:6].astype(np.float32),(nyb,nxb))[:,::-1]
                else:
                    winl1 = np.reshape(buff[:npix:6],(nyb,nxb))
                    winr1 = np.reshape(buff[1:npix:6],(nyb,nxb))[:,::-1]
                    winl2 = np.reshape(buff[2:npix:6],(nyb,nxb))
                    winr2 = np.reshape(buff[3:npix:6],(nyb,nxb))[:,::-1]
                    winl3 = np.reshape(buff[4:npix:6],(nyb,nxb))
                    winr3 = np.reshape(buff[5:npix:6],(nyb,nxb))[:,::-1]

                # For the reasons outlined in Rhead, we actually want to chop up
                # these 2 "data windows" into 6 per CCD. This is what we do next:

                # overscan is arranged as
                # 24 columns on LH of LH window
                #  4 columns on RH of LH window
                #  4 columns on LH of RH window
                # 24 columns on RH of RH window
                #  8 rows along top of LH and RH windows

                # Window 1 of the six comes from lower-left of left-hand data window
                w = self.win[0]
                xoff = 24 // xbin
                wins1.append(Window(winl1[:w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins2.append(Window(winl2[:w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins3.append(Window(winl3[:w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))

                # Window 2 comes from lower-right of right-hand data window
                w = self.win[1]
                xoff = 4 // xbin
                wins1.append(Window(winr1[:w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins2.append(Window(winr2[:w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins3.append(Window(winr3[:w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))

                # Window 3 is bias associated with left-hand data window (leftmost 24 and rightmost 4)
                w = self.win[2]
                lh = 24 // xbin
                rh = 4 // xbin
                wins1.append(Window(np.concatenate( (winl1[:w.ny,:lh],
                                                     winl1[:w.ny,-rh:]),axis=1),w.llx,w.lly,xbin,ybin))
                wins2.append(Window(np.concatenate( (winl2[:w.ny,:lh],
                                                     winl2[:w.ny,-rh:]),axis=1),w.llx,w.lly,xbin,ybin))
                wins3.append(Window(np.concatenate( (winl3[:w.ny,:lh],
                                                     winl3[:w.ny,-rh:]),axis=1),w.llx,w.lly,xbin,ybin))

                # Window 4 is bias associated with right-hand data window (leftmost 4 and rightmost 24)
                w = self.win[3]
                lh = 4 // xbin
                rh = 24 // xbin
                wins1.append(Window(np.concatenate( (winr1[:w.ny,:lh],
                                                     winr1[:w.ny,-rh:]),axis=1),w.llx,w.lly,xbin,ybin))
                wins2.append(Window(np.concatenate( (winr2[:w.ny,:lh],
                                                     winr2[:w.ny,-rh:]),axis=1),w.llx,w.lly,xbin,ybin))
                wins3.append(Window(np.concatenate( (winr3[:w.ny,:lh],
                                                     winr3[:w.ny,-rh:]),axis=1),w.llx,w.lly,xbin,ybin))

                # Window 5 comes from top strip of left-hand data window
                w = self.win[4]
                xoff = 24 // xbin
                yoff = 1024 // ybin
                wins1.append(Window(winl1[yoff:yoff+w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins2.append(Window(winl2[yoff:yoff+w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins3.append(Window(winl3[yoff:yoff+w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))

                # Window 6 comes from top of right-hand data window
                w = self.win[5]
                xoff = 4 // xbin
                yoff = 1024 // ybin
                wins1.append(Window(winr1[yoff:yoff+w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins2.append(Window(winr2[yoff:yoff+w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))
                wins3.append(Window(winr3[yoff:yoff+w.ny,xoff:xoff+w.nx],w.llx,w.lly,xbin,ybin))

            # Build the CCDs
            ccd1 = CCD(wins1, time, self.nxmax, self.nymax, True, None)
            ccd2 = CCD(wins2, time, self.nxmax, self.nymax, True, None)
            ccd3 = CCD(wins3, blueTime, self.nxmax, self.nymax, not badBlue, None)

            # Return a UCAM object
            return UCAM([ccd1,ccd2,ccd3], head)

        elif self.instrument == 'ULTRASPEC':

            wins = []
            noff = 0
            if self.mode.startswith('USPEC'):
                for w in self.win:

                    npix = w.nx*w.ny

                    # chop at left edge
                    nchop = max(0,17-w.llx)
                    nchop = nchop // xbin if nchop % xbin == 0 else nchop // xbin + 1

                    llx = max(1, w.llx + nchop*xbin - 16) if self.output == 'N' else \
                        max(1, 1074 - w.llx - w.nx*xbin)

                    if self.output == 'N':
                        # normal output, multi windows.
                        if flt:
                            wins.append(
                                Window(np.reshape(buff[noff:noff+npix].astype(np.float32),
                                                  (w.ny,w.nx))[:,nchop:],llx,w.lly,xbin,ybin))
                        else:
                            wins.append(
                                Window(np.reshape(buff[noff:noff+npix],(w.ny,w.nx))[:,nchop:],
                                       llx,w.lly,xbin,ybin))

                    elif self.output == 'A':
                        # avalanche output, multi windows.
                        if flt:
                            wins.append(
                                Window(np.reshape(buff[noff:noff+npix].astype(np.float32),
                                                  (w.ny,w.nx))[:,nchop::-1],llx,w.lly,xbin,ybin))
                        else:
                            wins.append(
                                Window(np.reshape(buff[noff:noff+npix],
                                                  (w.ny,w.nx))[:,nchop::-1],llx,w.lly,xbin,ybin))

                    noff += npix

            elif self.mode == 'UDRIFT':

                # drift mode Need to compute for left and right windows
                wl, wr = self.win
                npix = wl.nx*wl.ny + wr.nx*wr.ny

                # chop at left edge
                nchopl = max(0,17-wl.llx)
                nchopl = nchopl // xbin if nchopl % xbin == 0 else nchopl // xbin + 1

                nchopr = max(0,17-wr.llx)
                nchopr = nchopr // xbin if nchopr % xbin == 0 else nchopr // xbin + 1

                llxl = max(1, wl.llx + nchopl*xbin - 16) if self.output == 'N' else \
                    max(1, 1074 - wl.llx - wl.nx*xbin)
                llxr = max(1, wr.llx + nchopr*xbin - 16) if self.output == 'N' else \
                    max(1, 1074 - wr.llx - wr.nx*xbin)

                if self.output == 'N':
                    # normal output, drift
                    if flt:
                        comb = np.reshape(buff[:npix].astype(np.float32),(wl.ny,wl.nx+wr.nx))
                    else:
                        comb = np.reshape(buff[:npix],(wl.ny,wl.nx+wr.nx))

                elif self.output == 'A':
                    # avalanche output, drift
                    if flt:
                        comb = np.reshape(buff[:npix].astype(np.float32)[:,::-1],(wl.ny,wl.nx+wr.nx))
                    else:
                        comb = np.reshape(buff[:npix],(wl.ny,wl.nx+wr.nx)[:,::-1])

                wins.append(Window(comb[:,nchopl:wl.nx],llxl,wl.lly,xbin,ybin))
                wins.append(Window(comb[:,wl.nx+nchopr:],llxr,wl.lly,xbin,ybin))

            if self._ccd:
                return CCD(wins, time, self.nxmax, self.nymax, True, head)
            else:
                return MCCD([CCD(wins, time, self.nxmax, self.nymax, True, head),], head)

        else:
            raise UltracamError('Rdata.__init__: have not implemented anything for ' + self.instrument)

    def ntotal(self):
        """
        Returns total number of frames in data file
        """
        if self.server:
            ntot = get_nframe_from_server(self.run)
        else:
            self._fobj.seek(0,2)
            ntot = self._fobj.tell() // self.framesize
            self._fobj.seek(self.framesize*(self._nf-1))

        return ntot

    def nframe(self):
        """
        Returns next frame number to be read if reading
        sequentially (starts at 1)
        """
        return self._nf

    def time(self, nframe=None):
        """
        Returns timing information of frame nframe (starts from 1). This saves
        effort reading the data in some cases. See for example the ustats script.
        Once done, it moves to the next frame.

        nframe -- frame number to get, starting at 1. 0 for the last
                  (complete) frame.

        See utimer for what gets returned by this. See also Rtime for a class
        dedicated to reading times only.
        """

        # position read pointer
        self.set(nframe)

        if self.server:
            # somewhat inefficiently in the server case, we have to read the whole
            # frame because there are no options for timing data alone, although
            # at least no data re-formatting is required.
            full_url = URL + self.run + '?action=get_frame&frame=' + str(self._nf-1)
            buff     = urllib.request.urlopen(full_url).read()
            if len(buff) != self.framesize:
                self._nf = 1
                raise UltracamError('Rdata.time: failed to read frame ' + str(self._nf) +
                                    ' from FileServer. Buffer length vs expected = '
                                    + str(len(buff)) + ' vs ' + str(self.framesize) + ' bytes.')
            tbytes = buff[:2*self.headerwords]
        else:
            # read timing bytes alone
            tbytes = self._fobj.read(2*self.headerwords)
            if len(tbytes) != 2*self.headerwords:
                self._fobj.seek(0)
                self._nf = 1
                raise UendError('Rdata.time: failed to read timing bytes')

        tinfo = utimer(tbytes, self, self._nf)

        # step to start of next frame
        if not self.server:
            self._fobj.seek(self.framesize-2*self.headerwords,1)

        # move frame counter on by one
        self._nf += 1

        return tinfo

class Rtime (Rhead):
    """
    Iterator class to enable swift reading of Ultracam times.
    """
    def __init__(self, run, nframe=1, server=False):
        """
        Connects to a raw data file for reading. The file is kept open.
        The file pointer is set to the start of frame nframe.

        Arguments:

        run     -- as in 'run036'. Will try to access equivalent .xml and .dat
                   files

        nframe  -- frame to position for next read, starting at 1 as the first.

        server  -- True/False for server vs local disk access
        """
        Rhead.__init__(self, run, server)
        if self.isPonoff():
            raise PowerOnOffError('Rtime.__init__: attempted to read a power on/off')

        # Attributes set are:
        #
        # _fobj   -- file object opened on data file
        # _nf     -- next frame to be read
        # _run    -- name of run
        if server:
            self._fobj   = None
        else:
            self._fobj   = open(run + '.dat', 'rb')
        self._nf     = nframe
        self._run    = run
        if not server and nframe != 1:
            self._fobj.seek(self.framesize*(nframe-1))

    def __iter__(self):
        """
        Generator to allow Rtime to function as an iterator
        """
        try:
            while 1:
                yield self.__call__()
        except UendError:
            pass
        except urllib.error.HTTPError:
            pass

    def set(self, nframe=1):
        """
        Sets the internal file pointer to point at frame nframe.

        nframe  -- frame number to get, starting at 1. 0 for the last (complete) frame. 'None'
                   will be ignored. A value < 0 will cause an exception. A value greater than
                   the number of frames in the file will work, but will cause an exception to
                   be raised on the next attempted read.
        """

        # position read pointer
        if nframe is not None:
            if nframe < 0:
                raise UltracamError('Rtime.set: nframe < 0')
            elif nframe == 0:
                if self.server:
                    self._nf = get_nframe_from_server(self.run)
                else:
                    self._fobj.seek(0,2)
                    fp = self._fobj.tell()
                    nf = fp // self.framesize
                    self._fobj.seek(self.framesize*(nf-1)-fp,2)
                    self._nf = nf
            elif self._nf != nframe:
                if not self.server:
                    self._fobj.seek(self.framesize*(nframe-1))
                self._nf = nframe

    def __call__(self, nframe=None):
        """
        Returns timing information of frame nframe (starts from 1).

        nframe -- frame number to get, starting at 1. 0 for the last
                  (complete) frame.

        See utimer for gets returned by this.
        """

        # position read pointer
        self.set(nframe)

        if self.server:
            # have to read both timing and data in one go from the server
            # and just ignore the data
            full_url = URL + self.run + '?action=get_frame&frame=' + str(self._nf-1)
            buff     = urllib.request.urlopen(full_url).read()
            if len(buff) != self.framesize:
                self._nf = 1
                raise UltracamError('Rtime.__call__: failed to read frame ' + str(self._nf) +
                                    ' from FileServer. Buffer length vs expected = '
                                    + str(len(buff)) + ' vs ' + str(self.framesize) + ' bytes.')

            # have data. Re-format into the timing bytes and unsigned 2 byte int data buffer
            tbytes = buff[:2*self.headerwords]
        else:
            # read timing bytes
            tbytes = self._fobj.read(2*self.headerwords)
            if len(tbytes) != 2*self.headerwords:
                self._fobj.seek(0)
                self._nf = 1
                raise UendError('Data.get: failed to read timing bytes')

            # step to start of next frame
            self._fobj.seek(self.framesize-2*self.headerwords,1)

        tinfo = utimer(tbytes, self, self._nf)

        # move frame counter on by one
        self._nf += 1

        return tinfo

    def close_file(self):
        """
        Closes the file connected to the Rtime object to allow
        other reads to go ahead
        """
        self._fobj.close()

def utimer(tbytes, rhead, fnum):
    """
    Computes the Time corresponding of the most recently read frame,
    None if no frame has been read. For the Time to be reliable
    (Time.good == True), several frames might have needed to
    be read and their times passed through tstamp.

     tbytes  -- string of timing bytes

     rhead   -- Rhead header of data file

     fnum    -- frame number we think we are on.

    Returns (time,info,blueTime,badBlue) for ULTRACAM or (time,info) for ULTRASPEC

     time         : the Time as best as can be determined

     info         : a dictionary of extras to allow for possible future
                    updates without breaking code. Currently returns values for the
                    following keys:

         nsat         -- number of satellites (if set)
         format       -- the format integer used to translate the timing bytes
         vclock_frame -- vertical clocking time for a whole frame [ULTRACAM only]
         whichRun     -- special case run identifier. MAY2002 or nothing
         defTstamp    -- whether the "default" time stamping cycle was thought to apply
         gps          -- the raw GPS time associated with the frame, no corrections applied
         frameError   -- was there a frame numbering clash
         midnightCorr -- was the midnight bug correction applied

     ULTRACAM only:

     blueTime     : different time for the blue frames of ULTRACAM for nblue > 1

     badBlue      : blue frame is bad for nblue > 1 (nblue-1 out nblue are bad)

    """

    # This is an involved routine owing to the various hardware
    # changes and bugs that have cropped up over the years. The
    # pipeline equivalent routine is read_header.cc

    # return immediately if no bytes have been read.
    if tbytes is None: return None

    # start by assuming the time will be good, but many things
    # can wreck this. Once False, it can never return to True
    goodTime = True
    reason   = ''

    if rhead.instrument == 'ULTRASPEC' and rhead.version == -1:
        format = 2
    elif rhead.version == -1 or rhead.version == 70514 or \
            rhead.version == 80127:
        format = 1
    elif rhead.version == 100222 or rhead.version == 110921 or \
            rhead.version == 111205 or rhead.version == 120716 or \
            rhead.version == 120813 or rhead.version == 130303 or \
            rhead.version == 130317 or rhead.version == 140331:
        format = 2
    else:
        raise UltracamError('Rdata._timing: version = ' + str(rhead.version) + ' unrecognised.')

    frameNumber = struct.unpack('<I', tbytes[4:8])[0]
    if frameNumber != fnum:
        warnings.warn('ultracam.utimer: run ' + rhead._run +
                      ' expected frame number = ' +
                      str(fnum) + ' but found ' + str(frameNumber))
        frameError = True
    else:
        frameError = False

    # initialize some attributes of utimer which are used as the equivalent
    # of static variables in C++. They are:
    #
    #  previousFrameNumber  -- frame number of immediately preceding frame read, 0 if none
    #  tstamp               -- list of raw GPS times of preceding frames, [0] most recent
    #  blueTimes            -- list of modified times of preceding frames, [0] most recent

    if hasattr(utimer,'previousFrameNumber'):
        if frameNumber != utimer.previousFrameNumber + 1 or utimer.run != rhead._run:
            utimer.tstamp    = []
            utimer.blueTimes = []
    else:
        utimer.tstamp    = []
        utimer.blueTimes = []
    utimer.run = rhead._run

    utimer.previousFrameNumber = frameNumber

    if format == 1:
        nsec, nnsec = struct.unpack('<II', tbytes[9:17])
        nsat = struct.unpack('<h', tbytes[21:23])[0]
        if nsat <= 2:
            goodTime = False
            reason   = 'too few satellites (' + str(nsat) + ')'
        IMAX = struct.unpack('<I', '\xff\xff\xff\xff')[0]
        if nsec  == IMAX: nsec = 0
        if nnsec == IMAX: nnsec = 0

    elif format == 2:
        nexp  = struct.unpack('<I', tbytes[8:12])[0]
        if nexp*rhead.timeUnits != rhead.exposeTime:
            goodTime = False
            reason = 'XML expose time does not match time in timing bytes.'
        nsec, nnsec = struct.unpack('<II', tbytes[12:20])
        nnsec *= 100
        nsat   = None
        tstamp = struct.unpack('<H', tbytes[24:26])[0]

        if goodTime and (tstamp & PCPS_ANT_FAIL):
            goodTime = False
            reason   = 'GPS antenna failed'

        if goodTime and (tstamp & PCPS_INVT):
            goodTime = False
            reason   = 'GPS battery disconnected'

        if goodTime and not (tstamp & PCPS_SYNCD):
            goodTime = False
            reason = 'GPS clock not yet synced since power up'

        if goodTime and (tstamp & PCPS_FREER):
            goodTime = False
            reason = 'GPS receiver has not verified its position'

    else:
        raise UltracamError('Rdata.time: format = ' + str(format) + ' not recognised.')

    def tcon1(offset, nsec, nnsec):
        """
        Convert to MJD
        """
        return offset + float(nsec+nnsec/1.e9)/DSEC

    def tcon2(year, month, day, nsec, nnsec):
        """
        Convert to MJD
        """
        return (datetime.date(year,month,day).toordinal() - MJD0) + float(nsec+nnsec/1.e9)/DSEC

    if rhead.instrument == 'ULTRACAM':
        # One of the bits in the first byte is set if the blue frame is junk.
        # Unfortunately which bit was set changed hence the check of the format
        fbyte = six.indexbytes(tbytes,0)
        badBlue = rhead.nblue > 1 and \
            ((format == 1 and bool(fbyte & 1<<3)) or (format == 2 and bool(fbyte & 1<<4)))

    if format == 1 and nsat == -1:
        goodTime = False
        reason = 'no satellites.'
        mjd = tcon1(DEFDAT, nsec, nnsec)
        if rhead.instrument == 'ULTRACAM':
            if rhead.v_ft_clk > 127:
                vclock_frame = 6.e-9*(40+320*(rhead.v_ft_clk - 128))
            else:
                vclock_frame = 6.e-9*(40+40*rhead.v_ft_clk)
        day, month, year = struct.unpack('<BBH',tbytes[17:21])

    else:

        if rhead.whichRun == 'MAY2002' and format == 1:
            # Had no date info in the timestamps of this run
            # nsec was offset from 12 May 2002
            mjd = tcon1(MAY2002, nsec, nnsec)

            # Some of the nights ran over the end of the week. Spot
            # from dates prior to the start of the run on 16 May.
            if mjd < MAY2002+4: mjd += 7

            # Correct 10 second error that affected the May 2002 run.
            # Only possible if we have read the previous frames, with
            # time stored in an attribute called tstamp of this function
            if len(utimer.tstamp) and mjd < utimer.tstamp[0]: mjd += 10./DSEC

            # Fix problem with very first night
            if mjd < MAY2002+5.5:
                vclock_frame = 10.0e-6
            else:
                vclock_frame = 24.46e-6

        else:

            # OK now have date info, although we
            # need a stack of special case fixes
            if format == 1:
                day, month, year = struct.unpack('<BBH',tbytes[17:21])

                if month == 9 and year == 263: year = 2002
                if year < 2002:
                    mjd = tcon1(SEP2002, nsec, nnsec)

                elif month == 9 and year == 2002:
                    tdiff = datetime.date(year,month,day).toordinal()-MJD0-SEP2002
                    nweek = tdiff // 7
                    days  = tdiff - 7*nweek
                    if days > 3 and nsec < 2*DSEC:
                        nweek += 1
                    elif days <= 3 and nsec > 5*DSEC:
                        nweek -= 1
                    mjd = tcon1(SEP2002+7*nweek,nsec,nnsec)
                else:
                    mjd = tcon2(year, month, day, nsec % DSEC, nnsec)

            elif format == 2:
                mjd = tcon1(UNIX, nsec, nnsec)

            if rhead.instrument == 'ULTRACAM':
                if mjd > TSTAMP_CHANGE1:
                    if rhead.v_ft_clk > 127:
                        vclock_frame = 6.e-9*(40+320*(rhead.v_ft_clk - 128))
                    else:
                        vclock_frame = 6.e-9*(40+40*rhead.v_ft_clk)

                else:
                    if rhead.v_ft_clk > 127:
                        vclock_frame = 6.e-9*(80+160*(rhead.v_ft_clk - 128))
                    else:
                        vclock_frame = 6.e-9*(80+20*rhead.v_ft_clk)

    # 'midnight bug' correction
    if (int(mjd-3) % 7) == ((nsec // DSEC) % 7):
        warnings.warn('ultracam.utimer: run ' + rhead._run + ' midnight bug detected and corrected')
        mjd += 1
        midnightCorr = True
    else:
        midnightCorr = False

    # save this as the raw GPS time.
    gps = mjd

    # next variable determines when the timestamp is assumed to be taken within the
    # read cycle which has changed a few times owing to various small mishaps.
    defTstamp = mjd < TSTAMP_CHANGE1 or (mjd > TSTAMP_CHANGE2 and mjd < TSTAMP_CHANGE3)

    # Push time to front of tstamp list
    utimer.tstamp.insert(0,mjd)

    if rhead.instrument == 'ULTRACAM':
        VCLOCK_STORAGE = vclock_frame
        if rhead.gainSpeed == 'cdd':
            cds_time = CDS_TIME_CDD
        elif rhead.gainSpeed == 'fbb':
            cds_time = CDS_TIME_FBB
        elif rhead.gainSpeed == 'fdd':
            cds_time = CDS_TIME_FDD
        else:
            raise UltracamError('ultracam.utimer: did not recognize gain speed setting = ' + rhead.gainSpeed)
        VIDEO = SWITCH_TIME + cds_time

    elif rhead.instrument == 'ULTRASPEC':
        cdsTime = 0.
        # one extra parameter in addition to those from Vik's
        USPEC_FT_TIME  = 0.0067196 if mjd < USPEC_CHANGE else 0.0149818


    if rhead.instrument == 'ULTRACAM' and \
            (rhead.mode == 'FFCLR' or rhead.mode == 'FFOVER' or rhead.mode == '1-PCLR'):

        # never need more than two times
        if len(utimer.tstamp) > 2: utimer.tstamp.pop()
        ntmin = 2

        if defTstamp:
            mjdCentre  = utimer.tstamp[0]
            mjdCentre += rhead.exposeTime/DSEC/2.
            exposure   = rhead.exposeTime

        else:

            # Time taken to clear CCD
            clearTime   = (1033. + 1027)*vclock_frame

            # Time taken to read CCD (assuming cdd mode)
            if rhead.mode == 'FFCLR':
                readoutTime = (1024/rhead.ybin)*(VCLOCK_STORAGE*rhead.ybin +
                                                  536*HCLOCK + (512/rhead.xbin+2)*VIDEO)/1.e6
            elif rhead.mode == 'FFOVER':
                readoutTime = (1032/rhead.ybin)*(VCLOCK_STORAGE*rhead.ybin +
                                                  540*HCLOCK + (540/rhead.xbin+2)*VIDEO)/1.e6
            else:
                nxb          = rhead.win[1].nx
                nxu          = rhead.xbin*nxb
                nyb          = rhead.win[1].ny
                xleft        = rhead.win[0].llx
                xright       = rhead.win[1].llx + nxu - 1
                diff_shift   = abs(xleft - 1 - (1024 - xright) )
                num_hclocks  =  nxu + diff_shift + (1024 - xright) + 8 \
                    if (xleft - 1 > 1024 - xright) else nxu + diff_shift + (xleft - 1) + 8
                readoutTime = nyb*(VCLOCK_STORAGE*rhead.ybin +
                                    num_hclocks*HCLOCK + (nxb+2)*VIDEO)/1.e6

            # Frame transfer time
            frameTransfer = 1033.*vclock_frame

            if len(utimer.tstamp) == 1:

                # Case where we have not got a previous timestamp. Hop back over the
                # readout and frame transfer and half the exposure delay
                mjdCentre  = utimer.tstamp[0]
                mjdCentre -= (frameTransfer+readoutTime+rhead.exposeTime/2.)/DSEC
                if goodTime:
                    goodTime = False
                    reason = 'no previous GPS time found in non-default mode'

            else:

                # Case where we have got previous timestamp is somewhat easier and perhaps
                # more reliable since we merely need to step forward over the clear time and
                # half the exposure time.
                mjdCentre  = utimer.tstamp[1]
                mjdCentre += (clearTime + rhead.exposeTime/2.)/DSEC

            exposure = rhead.exposeTime

    elif rhead.instrument == 'ULTRACAM' and \
            (rhead.mode == 'FFNCLR' or rhead.mode == '1-PAIR' or rhead.mode == 'FFOVNC'
             or rhead.mode == '2-PAIR' or rhead.mode == '3-PAIR'):

        # never need more than three times
        if len(utimer.tstamp) > 3: utimer.tstamp.pop()
        ntmin = 3

        # Time taken to move 1033 rows.
        frameTransfer = 1033.*vclock_frame

        if rhead.mode == 'FFNCLR':
            readoutTime = (1024/rhead.ybin)*(VCLOCK_STORAGE*rhead.ybin +
                                              536*HCLOCK + (512/rhead.xbin+2)*VIDEO)/1.e6
        elif rhead.mode == 'FFOVNC':
                readoutTime = (1032/rhead.ybin)*(VCLOCK_STORAGE*rhead.ybin +
                                                  540*HCLOCK + (540/rhead.xbin+2)*VIDEO)/1.e6
        else:

            readoutTime = 0.
            xbin = rhead.xbin
            ybin = rhead.ybin
            ystart_old = -1
            for wl, wr in zip(rhead.win[::2],rhead.win[1::2]):

                nxu  = xbin*wl.nx
                nyu  = ybin*wl.ny

                ystart = wl.lly
                xleft  = wl.llx
                xright = wr.llx + nxu - 1

                if ystart_old > -1:
                    ystart_m = ystart_old
                    nyu_m    = nyu_old
                    y_shift  = (ystart-ystart_m-nyu_m)*VCLOCK_STORAGE
                else:
                    ystart_m = 1
                    nyu_m    = 0
                    y_shift  = (ystart-1)*VCLOCK_STORAGE

                # store for next time
                ystart_old = ystart
                nyu_old    = nyu

                # Number of columns to shift whichever window is further from
                # the edge of the readout to get ready for simultaneous
                # readout.
                diff_shift = abs(xleft - 1 - (1024 - xright) )

                # Time taken to dump any pixels in a row that come after the
                # ones we want.  The '8' is the number of HCLOCKs needed to
                # open the serial register dump gates If the left window is
                # further from the left edge than the right window is from the
                # right edge, then the diffshift will move it to be the same
                # as the right window, and so we use the right window
                # parameters to determine the number of hclocks needed, and
                # vice versa.
                num_hclocks = nxu + diff_shift + (1024 - xright) + 8 \
                    if (xleft - 1 > 1024 - xright) else nxu + diff_shift + (xleft - 1) + 8

                # Time taken to read one line. The extra 2 is required to fill the video pipeline buffer
                line_read = VCLOCK_STORAGE*ybin + num_hclocks*HCLOCK + (nxu/xbin+2)*VIDEO

                readoutTime += y_shift + (nyu/ybin)*line_read

            readoutTime /= 1.e6

        if defTstamp:
            if frameNumber == 1:
                mjdCentre  = utimer.tstamp[0]
                exposure   = rhead.exposeTime
                mjdCentre -= (frameTransfer+exposure/2.)/DSEC

            else:
                if len(utimer.tstamp) > 1:
                    texp = DSEC*(utimer.tstamp[0] - utimer.tstamp[1]) - frameTransfer
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += texp/2./DSEC
                    exposure   = texp

                else:
                    texp       = readoutTime + rhead.exposeTime
                    mjdCentre  = utimer.tstamp[0]
                    mjdCentre -= (frameTransfer+texp/2.)/DSEC
                    exposure   = texp

                    if goodTime:
                        goodTime = False
                        reason = 'could not establish an accurate time without previous GPS timestamp'

        else:
            if frameNumber == 1:
                mjdCentre  = utimer.tstamp[0]
                exposure   = rhead.exposeTime
                mjdCentre -= (frameTransfer+readoutTime+exposure/2.)/DSEC

                if goodTime:
                    goodTime = False
                    reason = 'cannot establish an accurate time for first frame in this mode'

            else:

                if len(utimer.tstamp) > 2:
                    texp       = DSEC*(utimer.tstamp[1] - utimer.tstamp[2]) - frameTransfer
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += (rhead.exposeTime - texp/2.)/DSEC
                    exposure   = texp

                elif len(utimer.tstamp) == 2:
                    texp = DSEC*(utimer.tstamp[0] - utimer.tstamp[1]) - frameTransfer
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += (rhead.exposeTime - texp/2.)/DSEC
                    exposure   = texp

                    if goodTime:
                        goodTime = False
                        reason = 'cannot establish an accurate time with only two prior timestamps'

                else:
                    texp       = readoutTime + rhead.exposeTime
                    mjdCentre  = utimer.tstamp[0]
                    mjdCentre += (rhead.exposeTime-texp-frameTransfer-texp/2.)/DSEC
                    exposure   = texp

                    if goodTime:
                        goodTime = False
                        reason = 'cannot establish an accurate time with only one prior timestamp'

    elif rhead.instrument == 'ULTRACAM' and rhead.mode == 'DRIFT':

        wl     = rhead.win[0]
        xbin   = rhead.xbin
        ybin   = rhead.ybin
        nxu    = xbin*wl.nx
        nyu    = ybin*wl.ny
        ystart = wl.lly
        xleft  = wl.llx
        wr     = rhead.win[1]
        xright = wr.llx + nxu -1

        # Maximum number of windows in pipeline
        nwins = int((1033./nyu+1.)/2.)
        pipe_shift = int(1033.-(((2.*nwins)-1.)*nyu))

        # Time taken for (reduced) frame transfer, the main advantage of drift mode
        frameTransfer = (nyu + ystart - 1)*vclock_frame

        # Number of columns to shift whichever window is further from the edge of the readout
        # to get ready for simultaneous readout.
        diff_shift = abs(xleft - 1 - (1024 - xright) )

        # Time taken to dump any pixels in a row that come after the ones we want.
        # The '8' is the number of HCLOCKs needed to open the serial register dump gates
        # If the left window is further from the left edge than the right window is from the
        # right edge, then the diffshift will move it to be the same as the right window, and
        # so we use the right window parameters to determine the number of hclocks needed, and
        # vice versa.
        num_hclocks  = nxu + diff_shift + (1024 - xright) + 8 \
            if (xleft - 1 > 1024 - xright) else nxu + diff_shift + (xleft - 1) + 8

        # Time taken to read one line. The extra 2 is required to fill the video pipeline buffer
        line_read = VCLOCK_STORAGE*ybin + num_hclocks*HCLOCK + (nxu/xbin+2)*VIDEO

        readoutTime = ((nyu/ybin)*line_read + pipe_shift*VCLOCK_STORAGE)/1.e6

        # Never need more than nwins+2 times
        if len(utimer.tstamp) > nwins+2: utimer.tstamp.pop()
        ntmin = nwins+2

        if defTstamp:

            # Pre board change or post-bug fix
            if len(utimer.tstamp) > nwins:
                texp = DSEC*(utimer.tstamp[nwins-1] - utimer.tstamp[nwins]) - frameTransfer
                mjdCentre  = utimer.tstamp[nwins]
                mjdCentre += texp/2./DSEC
                exposure   = texp

            else:

                # Set to silly value for easy checking
                mjdCentre = DEFDAT
                exposure  = rhead.exposeTime
                if goodTime:
                    goodTime = False
                    reason = 'too few stored timestamps for drift mode'

        else:

            if len(utimer.tstamp) > nwins+1:

                texp = DSEC*(utimer.tstamp[nwins] - utimer.tstamp[nwins+1]) - frameTransfer
                mjdCentre  = utimer.tstamp[nwins]
                mjdCentre += (rhead.exposeTime-texp/2.)/DSEC
                exposure   = texp

            elif len(utimer.tstamp) == nwins+1:

                texp       = DSEC*(utimer.tstamp[nwins-1] - utimer.tstamp[nwins]) - frameTransfer
                mjdCentre  = utimer.tstamp[nwins]
                mjdCentre += (rhead.exposeTime-texp/2.)/DSEC
                exposure   = texp

                if goodTime:
                    goodTime = False
                    reason = 'too few stored timestamps for drift mode'

            else:

                # Set to silly value for easy checking
                mjdCentre = DEFDAT
                exposure  = rhead.exposeTime

                if goodTime:
                    goodTime = False
                    reason = 'too few stored timestamps for drift mode'

    elif rhead.instrument == 'ULTRASPEC' and rhead.mode.startswith('USPEC'):

        # Avoid excessive accumulation of timestamps.
        if len(utimer.tstamp) > 3: utimer.tstamp.pop()
        ntmin = 3

        if utimer.tstamp[0] < USPEC_CHANGE:
            goodTime = False
            reason = 'timestamp too early'
            readoutTime = 0.
            texp = readoutTime + rhead.exposeTime
            mjdCentre = utimer.tstamp[0]
            if rhead.en_clr or frameNumber == 1:

                mjdCentre -= rhead.exposeTime/2./DSEC
                exposure   = rhead.exposeTime

            elif len(utimer.tstamp) > 1:

                texp = DSEC*(utimer.tstamp[0] - utimer.tstamp[1]) - USPEC_FT_TIME
                mjdCentre -= texp/2./DSEC
                exposure   = texp

            else:

                # Could be improved with an estimate of the read time
                mjdCentre -= rhead.exposeTime/2./DSEC
                exposure   = rhead.exposeTime

                if goodTime:
                    reason   = 'too few stored timestamps'
                    goodTime = False

        else:

            if rhead.en_clr or frameNumber == 1:
                # Special case for the first frame or if clears are enabled.
                exposure = rhead.exposeTime
                if len(utimer.tstamp) == 1:
                    mjdCentre = utimer.tstamp[0]
                    mjdCentre -= (-USPEC_FT_TIME-rhead.exposeTime/2.)/DSEC
                    if goodTime:
                        reason = 'cannot establish an accurate time without at least 1 prior timestamp'
                        goodTime = False
                else:
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += (USPEC_CLR_TIME+rhead.exposeTime/2.)/DSEC

            elif len(utimer.tstamp) > 2:

                # Can backtrack two frames to get a good exposure time.
                texp = DSEC*(utimer.tstamp[1] - utimer.tstamp[2]) - USPEC_FT_TIME
                mjdCentre  = utimer.tstamp[1]
                mjdCentre += (rhead.exposeTime-texp/2.)/DSEC
                exposure   = texp

            elif len(utimer.tstamp) == 2:

                # Can only back up one, so estimate of exposure time is
                # actually based on the exposure following the one of
                # interest. Probably not too bad, but technically unreliable
                # as a time.
                texp = DSEC*(utimer.tstamp[0] - utimer.tstamp[1]) - USPEC_FT_TIME
                mjdCentre  = utimer.tstamp[1]
                mjdCentre += (rhead.exposeTime-texp/2.)/DSEC
                exposure   = texp

                if goodTime:
                    reason = 'cannot establish an accurate time without at least 2 prior timestamps'
                    goodTime = False

            else:

                # Only one time
                mjdCentre  = utimer.tstamp[0]
                mjdCentre -= (rhead.exposeTime/2.+rhead.exposeTime)/DSEC
                exposure   = rhead.exposeTime

                if goodTime:
                    reason = 'too few stored timestamps'
                    goodTime = False

    elif rhead.instrument == 'ULTRASPEC' and rhead.mode == 'UDRIFT':

        ybin   = rhead.ybin
        nyu    = ybin*rhead.win[0].ny
        ystart = rhead.win[0].lly
        nwins  = int(((1037. / nyu) + 1.)/2.)
        frameTransfer = USPEC_FT_ROW*(ystart+nyu-1.)+USPEC_FT_OFF

        # Never need more than nwins+2 times
        if len(utimer.tstamp) > nwins+2: utimer.tstamp.pop()
        ntmin = nwins+2

        # Non-standard mode

        if len(utimer.tstamp) > nwins+1:

            texp       = DSEC*(utimer.tstamp[nwins] - utimer.tstamp[nwins+1]) - frameTransfer
            mjdCentre  = utimer.tstamp[nwins]
            mjdCentre += (rhead.exposeTime-texp/2.)/DSEC
            exposure   = texp

        elif len(utimer.tstamp) == nwins+1:

            texp          = DSEC*(utimer.tstamp[nwins-1] - utimer.tstamp[nwins]) - frameTransfer
            mjdCentre     = utimer.tstamp[nwins]
            mjdCentre     = (rhead.exposeTime-texp/2.)/DSEC
            exposure      = texp
            if goodTime:
                reason = 'too few stored timestamps'
                goodTime = False

        else:

            # Set to silly value for easy checking
            mjdCentre  = DEFDAT
            exposure   = rhead.exposeTime
            if goodTime:
                reason = 'too few stored timestamps'
                goodTime = False

    # Return some data
    time = Time(mjdCentre, exposure, goodTime, reason)

    if rhead.instrument == 'ULTRACAM':

        if rhead.nblue > 1:

            # The mid-exposure time for the OK blue frames in this case is computed by averaging the
            # mid-exposure times of all the contributing frames, if they are available.
            utimer.blueTimes.insert(0,time)

            if badBlue:

                # just pass through the standard time for the junk frames
                blueTime = time

            else:

                # if any of the contributing times is flagged as unreliable, then
                # so is the final time. This is also unreliable if any
                # contributing frame times are missing. Time is calculated as
                # half-way point between start of first and end of last
                # contributing exposure.  Corrections are made if there are too
                # few contributing exposures (even though the final value will
                # still be flagged as unreliable
                ncont  = min(rhead.nblue, len(utimer.blueTimes))
                start  = utimer.blueTimes[ncont-1].mjd - utimer.blueTimes[ncont-1].expose/2./DSEC
                end    = utimer.blueTimes[0].mjd       + utimer.blueTimes[0].expose/2./DSEC
                expose = DSEC*(end - start)

                # correct the times
                ok = ncont == rhead.nblue
                reason = ''
                if not ok:
                    expose *= rhead.nblue/float(ncont)
                    start   = end - expose/DSEC
                    reason  = 'not all contributing frames found'
                else:
                    ok = utimer.blueTimes[0].good and utimer.blueTimes[ncont-1].good
                    if not ok: reason  = 'time of start or end frame was unreliable'

                blueTime = Time((start+end)/2., expose, ok, reason)

            # Avoid wasting memory storing past times
            if len(utimer.blueTimes) > rhead.nblue: utimer.blueTimes.pop()

        else:
            blueTime = time

    # return lots of potentially useful extras in a dictionary
    info = {'nsat' : nsat, 'format' : format, 'whichRun' : rhead.whichRun,
            'defTstamp' : defTstamp, 'gps' : gps, 'frameError' : frameError,
            'midnightCorr' : midnightCorr, 'ntmin' : ntmin}

    if rhead.instrument == 'ULTRACAM':
        info['vclock_frame'] = vclock_frame
        return (time,info,blueTime,badBlue)

    elif rhead.instrument == 'ULTRASPEC':
        return (time,info)

    else:
        raise UltracamError('ultracam.utimer: did not recognize instrument = ' + rhead.instrument)
