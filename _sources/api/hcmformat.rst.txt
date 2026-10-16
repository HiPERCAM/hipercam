.. include:: ../globals.rst

|hiper|'s hcm file format
=========================

The |hiper| pipeline stores individual multi-CCD exposures in files with
extension '.hcm'. These are in fact FITS files that can be examined with
e.g. 'fv' or listed with 'fitshead', but the extension '.hcm' is used to
distinguish them clearly from the raw '.fits' files and to reduce the chances
of overwriting the latter. In order for these to be useful through the rest of
the software, some assumptions about their internal structure are necessary,
and they need to contain location information to show where sub-windows are
placed within the CCDs.  The purpose of this document is to specify this in as
much detail as possible.

The basic structure of an '.hcm' file is as follows::

  [General Header] [CCD 1] [CCD 2] [CCD 3] [CCD 4] [CCD 5]

'General Header' here is the first HDU [HDU0] which contains no data but has
general header information pertaining to the entire multi-CCD frame. Each CCD
in turn has the structure::

   Window 1, Window 2, Window 3, Window 4, ...

Each window has some header and data. The first has some extra headers general
to the CCD as a whole. Each window comes with a label (header item 'WINDOW').
For |hiper| data, these have values like 'E1', 'G2', 'E', 'F', 'G' and 'H'
denoting the quadrant of the CCD that the Window comes from. For a standard
windowed mode, with 1 window per quadrant, there will be 4 HDUs per CCD, so a
total of 5x4+1 = 25 HDUs. The Windows for a given CCD may come in any order,
but they should come as a contiguous block.

HDU structure and keywords
--------------------------

I now expand on this, detailing the header keywords expected. I concentrate on
those essential for the operation of pipeline commands, |reduce| above all.

1) HDU0. An hcm file starts with an HDU with no data but some header.

   Header items generated from HiPERACAM raw files:

      NUMCCD   : int
         specifies the number of CCDs.

      DATE     : string
         date the run file was written if from |hiper|

      GPS      : bool
         GPS status

      TIMSTAMP : string
         time stamp associated with the individual frame

      NFRAME   : int
         frame number

   Of these, 'TIMSTAMP' is essential.

2) Then each CCD comes as a block of HDUs, one per Window. The header in
   the first HDU *must* contain::

      NXTOT    : int
         maximum X-dimension, unbinned pixels to give a reference
         size for plots of this CCD.

      NYTOT    : int
         maximum Y-dimension, unbinned pixels.

      CCD      : string
         the label of the CCD

      MJDUTC   : float
         the MJD at mid-exposure of the CCD. This might be the same as
         TIMSTAMP of the general header, but it might not.

      GOODTIME : bool
         flag to indicate whether the time is thought to be reliable. Defaults
         to True if not present.

   It can also contain a keyword 'DSTATUS' giving a boolean that indicates
   whether the frame contains data or not (True means yes). This is needed
   because ULTRACAM and |hiper| have options to suspend readout ('NBLUE' and
   'NSKIP') which lead to blank dummy frames in between the real data
   frames. If 'DSTATUS' is not found, the CCD will be assumed to be OK.

   All the HDUs in a CCD *must* contain::

      CCD    : string
         string label matching the label of the first HDU of the CCD.
         This is used to identify which HDUs belong to a given CCD.

      LLX    : int
         X-coordinate of lower-left-most unbinned pixel that
         contributes to the window, starting with 1 on the far
         left.

      LLY    : int
         Y-coordinate of lower-left-most unbinned pixel that
         contributes to the window, starting with 1 on the bottom

      XBIN   : int
         binning factor in X

      YBIN   : int
         binning factor in Y

    Finally they *should* also contain::

      WINDOW : string
         a string label for the window. If it doesn't exist
         an attempt will be made to generate labels by counting.

Foreign data
------------

The format of '.hcm' files needs to be created if you want to apply the
pipeline to data from other instruments. Since such data is typically
FITS-format as are '.hcm' files, this might simply be a matter of adding new
keywords and re-arranging HDUs. The script for carrying this out is
|fits2hcm|. If you have some data, look at that script and submit new code
through github.

Example headers
---------------

Here is the result of running 'fitsheader' on a '.hcm' file which
may make the above easier to understand. Bear in mind there are many
extras associated with making it easier to display hcm files with 'ds9'.

.. literalinclude:: ../../data/example_headers

