.. Links to usefule files 21/05/2018

.. include:: globals.rst

.. highlightlang:: rest

Useful files
************

This page documents some files that may come in useful while observing or
during reduction.

CCD defects
===========

The |hiper| CCDs have a number of defects, dust and the like, and, as of
21 May 2018, a hair in the lower-left quadrant of CCD 4. When acquiring
targets these may not be visible so you are **strongly advised** to plot a file
of defects. You can also create your own using the pipeline command
|setdefect| and you are encouraged to do so if you see a defect not already
marked in whatever file you use.

Defects are classified as "moderate" (plotted in yellow) or "severe"
(red). What one counts as one or the other is a matter of taste. If you are
too picky, you could end up with so many defects that youare afraid to
observe.

I have decided to define defect file according to the following prescription:
``moderate``: 10 to 25% deviation from the norm; ``severe``: greater than 25%
from the norm. I mark each pixel in these categories, so a really bad feature
might appear in multiple pixels (there is a very nasty feature in the
lower-right of the lower-left quadrant of CCD 2 for example). I plan to update
the file with time, and retain old ones for monitoring purposes.  I
identified the bad pixels using a flat field divided by a smoothed version of
itself (smoothed with a 20 pixel FWHM 2D gaussian filter). The file used for
this is also supplied for reference.

When observing, you should move the target and essential comparison stars away
from any such defects, particularly those that appear in red.

Files:

  1. 2018-05-21 (:download:`hipercam defect file
     <hipercam_defects_2018_05_21.dft>`)

  2. 2018-06-02 (:download:`ultracam defect file
     <ultracam_defects_2018_06_02.dft>`)

.. Warning::
   For ULTRACAM I need also to add hot pixel defects which need to be regarded
   differently to flat field and poor charge transfer defects as a poor hot
   pixel could be very bad for a faint target but not matter much for a bright
   one.

