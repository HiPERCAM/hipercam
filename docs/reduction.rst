.. Reduction guide created on 25/10/2017

.. include:: globals.rst

Reducing |hiper| data
*********************

This guide assumes that you have got started with the |hiper| pipeline (see
:doc:`telescope` for a quick start to observing), but now would like a more
detailed guide to reducing your data. It covers the following steps:

  1. Creation of bias frames,
  2. Creation of flat field frames,
  3. Setting up aperture files,
  4. Defining the reduction file,
  5. Plotting the results.

Bias frames
===========

Bias frames can be quickly taken with |hiper|. All dome lights should be off,
the focal plane slide should be in to block light, and ideally the telescope
mirrors closed. Bias frames should be taken in clear mode with the shortest
possible exposures to minimize the time spent accumulating photons.

We standardly take 50 or 100 bias exposures. These can be combined by
averaging pixel-by-pixel with rejection of outliers to remove cosmic rays, or
by taking the pixel-by-pixel median. These operations can be carried out with
|combine| once individual exposures have been extracted with |grab|. It is
always advisable to inspect the frames visually with |rtplot| to check for
problems, e.g. with the readout, or light getting onto the images. So long as
several tens of bias exposures are available, then clipped mean combination of
frames is preferable to median combination because it leads to a lower level
statistical noise. Median combination works better for small number of frames
where the clipped mean can be less effective at removing outliers. Usually one
should combine bias frames with offsets to correct for any drift in the mean
level which could otherwise affect the action of |combine|.

The two operations of |grab| followed by |combine|, along with clean-up of the
temporary files can be carried out with the single command |makebias|. This
also saves the frames to a temporary location to avoid polluting the working
diectory with lots of files. Thus assuming all frames in bias run
:file:`run0002.fits` are OK, the following command will make the combined
bias frame::

  makebias run0002 1 0 3.0 yes run0002

rejecting pixels deviating by more that 3.0 sigma from the mean.  You might
also want to create a more memorable link to the output hcm file, depending
upon the particular type of bias::

  ln -s run0002.hcm bias-ff-slow.hcm

for example, for a full frame bias in slow readout mode.

.. Note::
   One should take bias calibration images in all the output formats used
   to guard against small changes that can occur as one changes output
   window formats. We have not yet established the significance of this
   for |hiper|.

.. Warning::
   Do not take bias frames too soon after (within less than 20 minutes)
   powering on the CCDs to avoid higher than normal dark current.

Flat fields
===========

Flat fields with |hiper| are taken at twilight in areas of sky with as few
stars as possible. One should offset the telescope between images to allow
stars to be removed. |hiper| takes images at high-speed so one can normally
aquire more than 100 frames in a single run, but the different CCDs will have
different count levels on any one frame, and will come out of saturation at
different times. The count levels will also be falling or rising according to
whether the flats were taken at evening or morning twilight.

The task of making the flat fields is thus to combine a series of frames with
differing count levels while removing features that vary between images. In
order to do this, one must normalise the images by their mean levels, but
weight them appropriately in the final combination to avoid giving too much
weight to under-exposed images. This is tedious by hand, and therefore the
command |makeflat| was written to carry out all the necessary tasks.

As with the biases, it is strongly recommended that you inspect the frames to
be combined using |rtplot| to avoid including any disastrous
ones. Saturated frames are spotted using user-defined mean levels at which to
reject frames. The documentation of |makeflat| has details of how it
works, and you are referred to this mor more information.

.. Warning::
   It is highly advisable to compute multiple versions of the flat field
   using different values of the parameter ``ngroup`` which can have a
   significant effect on the removal of stars from the final frame. See
   |makeflat|.

Aperture files
==============

The pipeline photometry is straightforward aperture photometry. Many of the
details can be defined when setting the apertures using |setaper|. A key
decision to be made at this stage is whether you think your target will
remain detectable on each frame throughout the run. Detectable means that
it's position can be measured and thus the relevant aperture re-positioned.
If not, then |setaper| gives you the option of ``linking`` any target to
another, with the idea that a brighter target can define the position shifts
which are applied to the fainter target.


Reduction files
===============

Once you have defined the apertures for all CCDs that you can, you can create
a reduction file using |genred|. This reads which CCDs and apertures have been
set and tries to create a file that at will work with |reduce|, even if it may
well not be ideal. Once you have this file, then you should expect to go
through a sequence of running |reduce|, adjusting the reduce file, re-running
|reduce| until you have a useful plot going.

Plotting results
================

TBD |plog|
