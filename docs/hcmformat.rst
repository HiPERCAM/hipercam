HiPERCAM's hcm file format
==========================

Individual multi-CCD exposures are stored in files with extension
'.hcm'. These are FITS files that can be examined with e.g. 'fv'.  In order
for these to be useful through the rest of the software, some assumptions
about their internal structure are necessary. The purpose of this document is
to specify this in as much detail as possible.

HDU structure
-------------

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

   None of these are essential.

2) Then each CCD has a structure of an HDU with header but no data, followed
   by an HDU per window. In this case with 4 windows, we would expect 5 HDUs
   per CCD therefore.

   The first header-only HDU *must* contain::

      NXTOT  : (int)
         maximum X-dimension, unbinned pixels.

      NYTOT  : (int)
         maximum Y-dimension, unbinned pixels. 

      CCD    : (string)
         the label of the CCD

      MJDUTC : (float)
         the MJD at mid-exposure of the CCD

      GOODTIM: (bool)
         flag to indicate whether the time is thought to be reliable

   It can also contain a keyword 'DSTATUS' giving a boolean that indicates
   whether the frame contains data or not (Trur means yes). This is needed
   because ULTRACAM and HiPERCAM have options to suspend readout ('NBLUE')
   which lead to blank dummy frames in between the real data frames. If
   'DSTATUS' is not found, then it will assumed to be OK.

   All HDUs after the first belonging to a given CCD *must* also contain
   the 'CCD' keyword with a label that matches that of the first header.
   This is used to identify which HDUs belong to the CCD.

   All HDUs after the first in a CCD *must* contain::

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
         an attempt will be made to generate one just by counting.
