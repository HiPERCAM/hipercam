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
targets these may not be visible so you are strongly advised to plot a file
of defects. You can also create your own using the pipeline command
|setdefect| and you are encouraged to do so if you see a defect not already
marked in whatever file you use.

Defects are classified as "moderate" (plotted in yellow) or "severe"
(red). What one counts as one or the other is a matter of taste. If you are
too picky, you could end up with so many defects that youare afraid to
observe.

I have defined an example defect file (see link below) according to the
following prescription: ``moderate``: 10 to 25% deviation from the norm;
``severe``: greater than 25% from the norm. I marked each pixel in these
categories, so a really bad feature might appear in multiple pixels. I plan
to update this file with time, and retain old ones for monitoring purposes.
I identified the bad pixels using a flat field divided by a smoothed version
of itself (smoothed with a 20 pixel FWHM 2D gaussian filter). The file used
for this is also supplied for reference.

Defect files:

  1. Date 2018-05-21. :download:`Defect file <hipercam_defects_2018_05_21.dft>`,
     :download:`Filtered flat field <filtered_flat_2018_05_21.hcm>`.

