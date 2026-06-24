.. Developer guide for third-party instrument plugins

.. include:: globals.rst

Writing an instrument plugin
****************************

The |hiper| pipeline supports third-party instruments through a plugin system
based on Python entry-points.  A plugin is a Python module that exposes a
small set of classes and a helper function, which the pipeline uses to read raw
data and present it as a stream of :class:`~hipercam.MCCD` objects.

This page explains how to write such a plugin.  As a concrete example
throughout, we assume an instrument called ``myinstrument`` that has a single
1024x1024 CCD and stores its data in a multi-extension FITS file: one primary
HDU (header only, carrying run-level metadata) followed by one image extension
per acquired frame.

Overview
========

The pipeline discovers plugins through the ``hipercam.instruments``
entry-point group.  When a data source string of the form
``'<instrument>:local'`` is passed to
:func:`~hipercam.spooler.data_source`, the pipeline calls
:func:`~hipercam.instruments.load_instrument` with the instrument name, which
looks up and imports the registered module.  The module must expose the
following names at module level:

==========================  ===========================================================
Name                        Purpose
==========================  ===========================================================
``Rdata``                   Data-reading class (iterable context manager)
``supports_server_access``  Boolean indicating if server access is supported
``get_ccd_pars``            Function returning CCD geometry as an ``OrderedDict``
==========================  ===========================================================

For a local-only plugin (i.e. one that does not support reading from a
network server) set supports_server_access to False and the ``server=True`` 
code paths should raise :exc:`NotImplementedError`.

The convenience base class :class:`~hipercam.instruments.RdataBase` provides
boilerplate ``__iter__``, ``__next__``, ``__enter__`` and ``__exit__``
implementations which allow subclasses to be used as iterators and context managers, 
and is the recommended base class for your ``Rdata``.

Package setup
=============

Create a Python package with the following minimal layout::

    myinst_plugin/
        src/
            myinst_plugin/
                __init__.py
                instrument.py
        pyproject.toml

In ``pyproject.toml`` declare the entry-point that registers your instrument
with the pipeline::

    [build-system]
    requires = ["setuptools>=61.0"]
    build-backend = "setuptools.build_meta"

    [project]
    name = "myinst_plugin"
    version = "0.1.0"
    dependencies = ["hipercam", "astropy", "numpy"]

    [project.entry-points."hipercam.instruments"]
    myinstrument = "myinst_plugin"

After installing your package the pipeline will automatically find the plugin by name.

.. note::
   If you are using ``pip install -e .`` to install your package in 
   editable (development) mode, you must re-run ``pip install -e .`` 
   each time you change the entry-point declaration, even in editable 
   (development) mode, so that the entry-point metadata is regenerated.

In ``__init__.py`` you can put any module-level code you like, but you must 
at least import the module level objects ``Rdata``, ``get_ccd_pars`` and 
``supports_server_access`` which we will define in the ``example.py`` file::

    # src/myinst_plugin/__init__.py
    from .instrument import Rdata, get_ccd_pars, supports_server_access

Module contract
===============

Your module must expose the names listed in the overview table.  Below is the
minimum skeleton::

    # src/myinst_plugin/instrument.py

    from hipercam.instruments import IendError, RdataBase

    supports_server_access = False

    class Rdata(RdataBase):  # data reader — the main class
        ...

    def get_ccd_pars(resource, server=False):
        ...


Implementing Rdata
==================

``Rdata`` is the central class.  The pipeline expect an initialiser with 
the signature below::

    Rdata(fname, nframe, server, **kwargs)

where ``nframe`` (1-based) is the first frame to read and ``server`` indicates
if the data should be read from a remote server or a local file.  The ``server``
parameter will always be ``False`` when accessed via ``'<instrument>:local'``.

Inheriting from :class:`~hipercam.instruments.RdataBase` provides
``__iter__``, ``__next__``, ``__enter__`` and ``__exit__`` for free.  You only
need to implement ``__call__``  and ``__del__``.

Constructor
-----------
In the snippet below we show how to write the constructor for a local-only plugin
that reads from a multi-extension FITS file.  The constructor opens the file and
stores the file handle and total number of frames as instance attributes.  The
``__del__`` method ensures that the file is closed when the object is destroyed.

::

    import numpy as np
    from astropy.io import fits as fits
    from hipercam.instruments import IendError, RdataBase

    class Rdata(RdataBase):
        """Reads frames from a myinstrument raw data file."""

        def __init__(self, fname, nframe=1, server=False, **kwargs):
            if server:
                raise NotImplementedError(
                    "myinstrument plugin does not support server access"
                )
            if not fname.endswith(".fits"):
                fname += ".fits"
            self._fname = fname
            self._hdul  = fits.open(self._fname, memmap=True)
            # HDU 0 is header-only; frames start at HDU index 1
            self._ntotal = len(self._hdul) - 1
            # set the nframe attribute to the first frame to read (1-based)
            self.nframe  = nframe

        def __del__(self):
            try:
                self._hdul.close()
            except Exception:
                pass


Reading a frame
---------------

``__call__`` reads the frame at ``self.nframe`` (or at a specific frame if
``nframe`` is supplied), builds an :class:`~hipercam.MCCD`, and advances the
internal counter. ``__call__`` should return a :class:`~hipercam.MCCD` if 
possible and raise :class:`~hipercam.instruments.IendError` when all frames 
have been read. For data which is still being written to disc, or for server
access, ``__call__`` should return ``None`` to indicate that we are still
waiting for the data::

        def __call__(self, nframe=None):
            """Read one frame and return it as an MCCD.

            Raises IendError when all frames have been read.
            """
            # seek the requested frame if nframe is supplied; 
            # otherwise read the next frame in sequence
            if nframe is not None:
                self.nframe = nframe

            # for data that is already written, indicate we've reached 
            # the end of the file with an exception
            if self.nframe > self._ntotal:
                raise IendError("end of file reached")

            # get the header and data for the current frame.
            hdu   = self._hdul[self.nframe]  # 1-based: HDU 0 is the primary
            data  = hdu.data.astype(np.float32)
            fhead = hdu.header               # raw FITS header for this extension

            # delegate construction of the MCCD object to a helper method.
            mccd = self._build_mccd(data, fhead, self.nframe)

            # increment the frame counter so a subsequent call
            # (with nframe=None) will read the next frame in sequence.
            self.nframe += 1

            # return the MCCD object to the caller
            return mccd


Constructing MCCD objects
=========================

The :class:`~hipercam.MCCD` is the standard data container used throughout the
pipeline.  A :class:`~hipercam.MCCD` consists of a collection of one or more
:class:`~hipercam.CCD` objects, each of which contains one or more 
:class:`~hipercam.Window` objects. 

A :class:`~hipercam.CCD` can contain multiple windows to deal with instruments
that have more than one output amplifier, each reading out a different region 
of the CCD.

It is assembled from the inside out::

    `~hipercam.header.Header`  →  Winhead  →  Window  →  CCD  →  MCCD

Top-level header
----------------

The top-level header of the :class:`~hipercam.MCCD` is an
``hipercam.header.Header`` or ``~astropy.io.fits.Header`` object.  
The pipeline requires the following keywords to be present:

=============  ======  =======================================================
Keyword        Type    Description
=============  ======  =======================================================
``INSTRUME``   str     Instrument name
``TIMSTAMP``   str     UTC timestamp of the mid-exposure in ISO 8601 format,
                       e.g. ``'2026-04-16T21:30:00.000000000'``
``MJDINT``     int     Integer part of the UTC mid-exposure MJD
``MJDFRAC``    float   Fractional part of the UTC mid-exposure MJD
``MJDUTC``     float   Full MJD(UTC), equal to ``MJDINT + MJDFRAC``
``GOODTIME``   bool    ``True`` if the timestamp is reliable
``EXPTIME``    float   Exposure duration in seconds
``NFRAME``     int     Frame number within the run (1-based)
=============  ======  =======================================================

In addition to this, any other FITS keywords may be included as desired. 

Example::

    from astropy.io import fits
    thead = fits.Header()
    thead['INSTRUME'] = ('myinstrument', 'Instrument name')
    thead['TIMSTAMP'] = (fhead['DATE-OBS'], 'UTC timestamp')
    mjd     = fhead['MJD-OBS']
    mjdint  = int(mjd)
    mjdfrac = mjd - mjdint
    thead['MJDINT']   = (mjdint,  'Integer part of MJD(UTC)')
    thead['MJDFRAC']  = (mjdfrac, 'Fractional part of MJD(UTC)')
    thead['MJDUTC']   = (mjd,     'MJD(UTC)')
    thead['GOODTIME'] = (True,    'Timestamp reliable')
    thead['EXPTIME']  = (float(fhead['EXPTIME']), 'Exposure time (s)')
    thead['NFRAME']   = (nframe,  'Frame number')

Per-CCD header
--------------

There is no separate header object for a :class:`~hipercam.CCD`.  Instead,
CCD-level information — in particular per-CCD timing data — is stored in the
header of the **first** :class:`~hipercam.Window` in that CCD.  This header is
a :class:`hipercam.Header` object (a lightweight FITS-compatible header class)
and must carry the same timing keywords as the top-level header:

=============  ===============================================================
Keyword        Description
=============  ===============================================================
``TIMSTAMP``   UTC timestamp (ISO 8601)
``MJDINT``     Integer part of MJD(UTC)
``MJDFRAC``    Fractional part of MJD(UTC)
``MJDUTC``     Full MJD(UTC)
``GOODTIME``   Timestamp reliability flag
``EXPTIME``    Exposure time (s)
``NFRAME``     Frame number
=============  ===============================================================

It can also contain any other CCD-level keywords as desired (e.g filter names, 
CCD temperature, etc).   

::

    from hipercam.header import Header

    ccd_head = Header()
    ccd_head['TIMSTAMP'] = (fhead['DATE-OBS'], 'UTC timestamp')
    ccd_head['MJDINT']   = (mjdint,  'Integer part of MJD(UTC)')
    ccd_head['MJDFRAC']  = (mjdfrac, 'Fractional part of MJD(UTC)')
    ccd_head['MJDUTC']   = (mjd,     'MJD(UTC)')
    ccd_head['GOODTIME'] = (True,    'Timestamp reliable')
    ccd_head['EXPTIME']  = (float(fhead['EXPTIME']), 'Exposure time (s)')
    ccd_head['NFRAME']   = (nframe,  'Frame number')

Building Winhead and Window
---------------------------

A :class:`~hipercam.Winhead` describes the position and binning of a
rectangular readout region of the CCD. Windows may only take up a subset
of the full CCD area, and may be binned. For a full-frame 1024x1024 readout
with 1x1 binning::

    from hipercam.window import Winhead, Window

    # Winhead(llx, lly, nx, ny, xbin, ybin, outamp, head=ccd_head)
    #
    #   llx, lly  — lower-left corner in **unbinned pixels**, 1-based
    #   nx, ny    — window dimensions in **binned** pixels
    #   xbin,ybin — binning factors
    #   outamp    — output amplifier location; 'UL', 'UR', 'LL' etc. '' means unknown
    #   head      — per-CCD header created above, attached to
    #               the first (and here only) window
    #
    winhd = Winhead(1, 1, 1024, 1024, 1, 1, '', head=ccd_head)
    win   = Window(winhd, data=data)

.. note::

   If a CCD has more than one readout window, you need to create
   a :class:`~hipercam.Window` object for each one. Only the **first** 
   window in the :class:`~hipercam.CCD` carries the per-CCD timing header; 
   all subsequent windows should omit the `head` parameter.

Building CCD and MCCD
---------------------

A :class:`~hipercam.CCD` is assembled from one or more :class:`~hipercam.Window` 
objects. The :class:`~hipercam.Group` object is a simple wrapper around
an ``OrderedDict`` that we can use to create a group of windows::

    from hipercam import Group, CCD

    # CCD with total size 1024x1024, no pre/over-scan
    ccd_obj = CCD(Group(Window), nxtot=1024, nytot=1024, nxpad=0, nypad=0)
    # add the single window to the CCD with label '1'
    ccd['1'] = win

Finally we can create a group of one or more CCDs and use it to construct the
:class:`~hipercam.MCCD` by providing the group and the top-level header::

    ccds = Group(CCD)
    # add the single CCD to the group with label '1'
    ccds['1'] = ccd_obj
    mccd = MCCD(ccds, head=thead)
    return mccd

Complete ``_build_mccd`` helper
--------------------------------

Putting all of the above together into a single private method::

    from astropy.io import fits as fits
    from hipercam.header import Header
    from hipercam.group import Group
    from hipercam.window import Winhead, Window
    from hipercam.ccd import CCD, MCCD

    def _build_mccd(self, data, fhead, nframe):
        """Build an MCCD from a 2-D pixel array and a raw FITS header."""

        # --- timing values extracted from the raw per-frame header ---
        mjd     = fhead['MJD-OBS']
        mjdint  = int(mjd)
        mjdfrac = mjd - mjdint
        tstamp  = fhead['DATE-OBS']
        exptime = float(fhead['EXPTIME'])

        # --- top-level MCCD header (astropy fits.Header) ---
        thead = fits.Header()
        thead['INSTRUME'] = ('myinstrument', 'Instrument name')
        thead['TIMSTAMP'] = (tstamp,   'UTC timestamp')
        thead['MJDINT']   = (mjdint,   'Integer part of MJD(UTC)')
        thead['MJDFRAC']  = (mjdfrac,  'Fractional part of MJD(UTC)')
        thead['MJDUTC']   = (mjd,      'MJD(UTC)')
        thead['GOODTIME'] = (True,     'Timestamp reliable')
        thead['EXPTIME']  = (exptime,  'Exposure time (s)')
        thead['NFRAME']   = (nframe,   'Frame number')

        # --- per-CCD header, stored in the first window (hipercam.Header) ---
        ccd_head = Header()
        ccd_head['TIMSTAMP'] = (tstamp,  'UTC timestamp')
        ccd_head['MJDINT']   = (mjdint,  'Integer part of MJD(UTC)')
        ccd_head['MJDFRAC']  = (mjdfrac, 'Fractional part of MJD(UTC)')
        ccd_head['MJDUTC']   = (mjd,     'MJD(UTC)')
        ccd_head['GOODTIME'] = (True,    'Timestamp reliable')
        ccd_head['EXPTIME']  = (exptime, 'Exposure time (s)')
        ccd_head['NFRAME']   = (nframe,  'Frame number')

        # --- window, CCD, MCCD ---
        winhd   = Winhead(1, 1, 1024, 1024, 1, 1, '', head=ccd_head)
        win     = Window(winhd, data=data)

        # CCD with total size 1024x1024, no pre/over-scan
        ccd_obj = CCD(Group(Window), nxtot=1024, nytot=1024, nxpad=0, nypad=0)
        # add the single window to the CCD with label '1'
        ccd['1'] = win

        ccds    = Group(CCD)
        ccds['1'] = ccd_obj
        return MCCD(ccds, head=thead)

Implementing get_ccd_pars
=========================

:func:`get_ccd_pars` is called by the pipeline to determine CCD geometry
before any frames are read (for example, to set up display windows).  It must
return an :class:`~collections.OrderedDict` mapping each CCD label to a
4-tuple ``(nxmax, nymax, nxpad, nypad)``:

==========  =====================================================================
Value       Meaning
==========  =====================================================================
``nxmax``   Full unbinned X dimension of the CCD imaging area
``nymax``   Full unbinned Y dimension of the CCD imaging area
``nxpad``   Extra X pixels to allow for pre/over-scan regions (0 if none)
``nypad``   Extra Y pixels to allow for pre/over-scan regions (0 if none)
==========  =====================================================================

For ``myinstrument`` with its single 1024x1024 CCD and no scan regions::

    from collections import OrderedDict

    def get_ccd_pars(resource, server=False):
        """Return CCD geometry for myinstrument.

        Arguments::

           resource : str
              Run name.  Not used; myinstrument has fixed geometry.

           server : bool
              Unused; present for API compatibility.

        """
        return OrderedDict([('1', (1024, 1024, 0, 0))])

Using the plugin
================

Once the package is installed the plugin is available from any of the command-line
scripts. For API use, the plugin is accessed with the
``'<instrument>:local'`` source string::

    from hipercam.spooler import data_source, get_ccd_pars

    # Iterate through all frames in a file
    with data_source('myinstrument:local', 'run001', first=1) as spool:
        for mccd in spool:
            # mccd is a hipercam.MCCD object
            ccd1 = mccd['1']   # the single CCD, labelled '1'
            ...

    # Query CCD geometry
    pars = get_ccd_pars('myinstrument:local', 'run001')
    # returns OrderedDict([('1', (1024, 1024, 0, 0))])

Complete example
================

Below is a complete ``instrument.py`` for a single-CCD instrument storing
full-frame 1024×1024 data in multi-extension FITS format.

.. code-block:: python

    """
    myinst_plugin.instrument
    ~~~~~~~~~~~~~~~~~~

    HiPERCAM instrument plugin for myinstrument.

    File format assumed::

        HDU 0  -- primary, header-only (run-level metadata)
        HDU 1  -- first frame  (image data + per-frame FITS header)
        HDU 2  -- second frame
        ...

    Access via hipercam::

        from hipercam.spooler import data_source
        with data_source('myinstrument:local', 'run001') as spool:
            for mccd in spool:
                ...

    """

    from collections import OrderedDict

    import numpy as np
    from astropy.io import fits as fits

    from hipercam.header import Header
    from hipercam.group import Group
    from hipercam.window import Winhead, Window
    from hipercam.ccd import CCD, MCCD
    from hipercam.instruments import IendError, RdataBase

    # tell the hipercam pipeline that this instrument does not support server access
    supports_server_access = False

    class Rdata(RdataBase):
        """Reads frames from a myinstrument data file.

        Usage::

            with Rdata('run001', nframe=1) as rdat:
                for mccd in rdat:
                    ...

        """

        def __init__(self, fname, nframe=1, server=False, **kwargs):
            if server:
                raise NotImplementedError(
                    "myinstrument plugin does not support server access"
                )
            self._fname  = fname + ".fits"
            self._hdul   = fits.open(self._fname, memmap=True)
            self._ntotal = len(self._hdul) - 1
            self.nframe  = nframe

        def __call__(self, nframe=None):
            """Read one frame and return it as an MCCD.

            Arguments::

               nframe : int or None
                  Frame number to read (1-based).  If None, reads the next
                  frame in sequence.

            Raises :exc:`~hipercam.instruments.IendError` when the end of
            the file is reached.
            """
            # seek the requested frame if nframe is supplied; 
            # otherwise read the next frame in sequence
            if nframe is not None:
                self.nframe = nframe

            # for data that is already written, indicate we've reached 
            # the end of the file with an exception
            if self.nframe > self._ntotal:
                raise IendError("end of file reached")

            # get the header and data for the current frame.
            hdu   = self._hdul[self.nframe]  # 1-based: HDU 0 is the primary
            data  = hdu.data.astype(np.float32)
            fhead = hdu.header               # raw FITS header for this extension

            # delegate construction of the MCCD object to a helper method.
            mccd = self._build_mccd(data, fhead, self.nframe)

            # increment the frame counter so a subsequent call
            # (with nframe=None) will read the next frame in sequence.
            self.nframe += 1

            # return the MCCD object to the caller
            return mccd

        def __del__(self):
            try:
                self._hdul.close()
            except Exception:
                pass

        # -------------------------------------------------------------------
        # Internal helper
        # -------------------------------------------------------------------

        def _build_mccd(self, data, fhead, nframe):
        """Build an MCCD from a 2-D pixel array and a raw FITS header."""

            # --- timing values extracted from the raw per-frame header ---
            mjd     = fhead['MJD-OBS']
            mjdint  = int(mjd)
            mjdfrac = mjd - mjdint
            tstamp  = fhead['DATE-OBS']
            exptime = float(fhead['EXPTIME'])

            # --- top-level MCCD header (astropy fits.Header) ---
            thead = fits.Header()
            thead['INSTRUME'] = ('myinstrument', 'Instrument name')
            thead['TIMSTAMP'] = (tstamp,   'UTC timestamp')
            thead['MJDINT']   = (mjdint,   'Integer part of MJD(UTC)')
            thead['MJDFRAC']  = (mjdfrac,  'Fractional part of MJD(UTC)')
            thead['MJDUTC']   = (mjd,      'MJD(UTC)')
            thead['GOODTIME'] = (True,     'Timestamp reliable')
            thead['EXPTIME']  = (exptime,  'Exposure time (s)')
            thead['NFRAME']   = (nframe,   'Frame number')

            # --- per-CCD header, stored in the first window (hipercam.Header) ---
            ccd_head = Header()
            ccd_head['TIMSTAMP'] = (tstamp,  'UTC timestamp')
            ccd_head['MJDINT']   = (mjdint,  'Integer part of MJD(UTC)')
            ccd_head['MJDFRAC']  = (mjdfrac, 'Fractional part of MJD(UTC)')
            ccd_head['MJDUTC']   = (mjd,     'MJD(UTC)')
            ccd_head['GOODTIME'] = (True,    'Timestamp reliable')
            ccd_head['EXPTIME']  = (exptime, 'Exposure time (s)')
            ccd_head['NFRAME']   = (nframe,  'Frame number')

            # --- Create single Window object (only one readout amplifier) ---
            # Winhead describes window position and binning.
            # For a full-frame 1024x1024 readout with 1x1 binning:
            #   llx=1, lly=1      lower-left unbinned pixel (1-based)
            #   nx=1024, ny=1024  window size in binned pixels
            #   xbin=1, ybin=1    binning factors
            #   outamp=''         output amplifier location unknown
            winhd   = Winhead(1, 1, 1024, 1024, 1, 1, '', head=ccd_head)
            win     = Window(winhd, data=data)

            # CCD with total size 1024x1024, no pre/over-scan
            ccd_obj = CCD(Group(Window), nxtot=1024, nytot=1024, nxpad=0, nypad=0)
            # add the single window to the CCD with label '1'
            ccd['1'] = win

            ccds    = Group(CCD)
            ccds['1'] = ccd_obj
            return MCCD(ccds, head=thead)

    def get_ccd_pars(resource, server=False):
        """Return CCD geometry for myinstrument.

        Arguments::

           resource : str
              File name.  Not used; myinstrument has fixed geometry.

           server : bool
              Unused; present for API compatibility.

        Returns an :class:`~collections.OrderedDict` mapping CCD label
        ``'1'`` to ``(nxmax, nymax, nxpad, nypad)`` = ``(1024, 1024, 0, 0)``.
        """
        return OrderedDict([('1', (1024, 1024, 0, 0))])
