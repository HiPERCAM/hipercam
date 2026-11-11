.. documentation of HiPERCAM's photometry

.. include:: globals.rst

Stacking data with shiftadd
***************************

Very often one wants to combine a number of frames together to improve the
signal-to-noise of faint sources. Whilst |averun| can be used to do this, 
it simply combines the frames without any attempt to align them. For very 
deep images, usually one wants to align the frames first, and this is what 
|shiftadd| does. This page discusses the usage of |shiftadd| and provides
some tips on how to get the best results. 

.. Note::
   TLDR: For those in a hurry, stacking frames with |shiftadd| requires the following steps:

      1. Create a low signal-to-noise stack using |averun|.
      2. Create an aperture file using |setaper| and define reference apertures.
      3. Set parameters for aperture positioning and image calibration using |genred|.
      4. Run |shiftadd| to align and stack the frames.

.. contents:: Contents
   :local:

.. _aligning_images:

Aligning images
===============

The first step in stacking images is to find the offsets between the frames. 
Rather than re-write the machinery for doing this from scratch, |shiftadd|
uses the same code used for finding the positions of stars during the photometry.
Therefore, the first step is to proceed as if you are performing a  
:doc:`data reduction <reduction>`. 

In particular, you will want to create an :ref:`aperture file <aperture_files>`
using |setaper| and define some **reference apertures** using the `r` keyboard 
shortcut. These positions of these reference apertures in each frame will be 
used to determine the offsets between the frames. You need to place at least one
reference aperture for every CCD that you want to align. |shiftadd| does not
rotate frames (only applying linear shifts), so the primary advantage of using
multiple apertures is to improve the accuracy of the offsets. 

.. _reduction_file_for_parameter_setting:
Reduction file for parameter setting
====================================

There are many parameters that control how apertures are :ref:`positioned <aperture_positioning>`.
Because these are also common with the parameters used in photometry, we create a 
:ref:`reduction file <reduce_files>` to hold these parameters. This is done using the |genred|
command.

.. _calibrating_images:

Calibrating images
==================

Before stacking images, we usually want to apply the usual calibration steps to 
remove instrumental effects. This includes bias subtraction, flat fielding and fringe removal. 
|genred| will prompt you for the relevant calibration files to use in these steps. 


.. _stacking_images:

Stacking images
===============

The actual stacking of images is done using `reproject <https://reproject.readthedocs.io/en/stable/>`_
as an engine. The documentation for the reproject package is the best place to look for the details
of the various algorithms used to resample the data. Here we only outline a few issues which are
worth considering for different use cases.

To limit memory use for large runs, |shiftadd| stores reprojected windows in temporary
disk-backed arrays and combines them in spatial tiles. Processing each CCD requires
temporary space for one full-frame, double-precision array per input window in the
accepted frames. Temporary files are removed on normal completion, errors, or a
keyboard interrupt; an abrupt process termination may leave them behind. Set
the ``TMPDIR`` environment variable before running |shiftadd| to choose a
filesystem with enough free space. For small datasets, the hidden option
``memmap=False`` keeps reprojected windows in RAM instead, avoiding disk-backed
stack storage and blocked reprojection. This can be faster, but requires RAM
for one full-frame, double-precision array per input window on a CCD, in
addition to combination and reprojection working memory.
As with other hidden |shiftadd| parameters, the chosen value is remembered for
later runs; use ``memmap=True`` to switch back to disk-backed storage.

Exact pixel shifts
------------------
For very faint sources, you may want to avoid any resampling of the data, which can wash out the signal. 
In this case, you can use set the `reprmethod=interp` parameter with `reprorder=0` to apply integer pixel 
shifts only.  For almost all other cases you want to use `reprmethod=adaptive` or `reprmethod=exact` 
as these methods guarantee conservation of flux and can provide better resampling for brighter sources.

.. _instrument_specific_issues:
Instrument specific issues
==========================
There are some peculiarities of |shiftadd| that should be noted when used on |ultra| or |hiper| data.
These instruments have multiple CCDs, and each CCD has two or more output amplifiers with different
gains, bias levels and readout noise. The output amplifiers naturally divide the CCD into multiple
regions (2 in the case of ULTRACAM and 4 in the case of |hiper|). In addition, the observer can define 
*windows* which fill some, or all, of these regions.

The |shiftadd| command reprojects each region/window separately and masks data outside
its shifted footprint. However, the reprojected windows are combined on a common
full-frame grid. If shifted footprints of neighbouring windows from different frames
overlap, their pixels can be combined at the overlap, including at a window edge.
