HiPERCAM's hcm file format
==========================

Individual multi-CCD exposures are stored in files with extension
'.hcm'. These are FITS files that can be examined with e.g. 'fv'.  In order
for these to be useful through the rest of the software, some assumptions
about their internal structure are necessary, and they need to contain
location information to show where sub-windows are placed within the CCDs.
The purpose of this document is to specify this in as much detail as possible.

The basic structure is as follows::

  [General Header] [CCD 1] [CCD 2] [CCD 3] [CCD 4] [CCD 5]

'General Header' here is the first HDU [HDU0], and contains general
information pertaining to the entire multi-CCD frame. Each CCD
in turn has the structure::

   Header, Window E1, Window F1, Window H1, Window G1

in 5 HDUs, the first being header specific to the CCD, the rest being data
with one HDU per window. The labels (1 to 5, E1 to G1) are specific to
HiPERCAM in one window mode, and could differ in general. In this example
there would be a total of 1+(5x5) = 26 HDUs, from HDU0 to HDU25.

HDU structure and keywords
--------------------------

I now expand on this, detailing the header keywords expected.

1) HDU0. An hcm file starts with an HDU with no data but some header.

   Header items generated from HiPERACAM raw files::

      NUMCCD   : (int)
         specifies the number of CCDs.

      DATE     : (string)
         date the run file was written if from HiPERCAM

      GPS      : (bool)
         GPS status

      TIMSTAMP : (string)
         time stamp associated with the individual frame

      NFRAME   : (int)
         frame number

   Of these, 'TIMSTAMP' is essential.

2) Then each CCD has a structure of an HDU with header but no data, followed
   by an HDU per window. In this case with 4 windows, we would expect 5 HDUs
   per CCD therefore.

   The header in the first HDU *must* contain::

      NXTOT    : (int)
         maximum X-dimension, unbinned pixels to give a reference
         size for plots of this CCD.

      NYTOT    : (int)
         maximum Y-dimension, unbinned pixels. 

      CCD      : (string)
         the label of the CCD

      MJDUTC   : (float)
         the MJD at mid-exposure of the CCD. This might be the same as
         TIMSTAMP of the general header, but it might not.

      GOODTIME : (bool)
         flag to indicate whether the time is thought to be reliable

   It can also contain a keyword 'DSTATUS' giving a boolean that indicates
   whether the frame contains data or not (True means yes). This is needed
   because ULTRACAM and HiPERCAM have options to suspend readout ('NBLUE',
   'NSKIP') which lead to blank dummy frames in between the real data
   frames. If 'DSTATUS' is not found, the CCD will be assumed to be OK.

   All HDUs after the first in a CCD *must* contain::

      CCD    : (string)
         string label matching the label of the first HDU of the CCD.
         This is used to identify which HDUs belong to a given CCD.

      LLX    : (int)
         X-coordinate of lower-left-most unbinned pixel that
         contributes to the window, starting with 1 on the far
         left.

      LLY    : (int)
         Y-coordinate of lower-left-most unbinned pixel that
         contributes to the window, starting with 1 on the bottom

      XBIN   : (int)
         binning factor in X

      YBIN   : (int)
         binning factor in Y

    Finally they *should* also contain::

      WINDOW : (string)
         a string label for the window. If it doesn't exist
         an attempt will be made to generate labels by counting.

Example headers
---------------

Here is the result of running 'fitsheader' on a '.hcm' file which
may make the above easier to understand:

.. literalinclude:: ../data/example_headers

