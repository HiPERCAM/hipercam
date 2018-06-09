.. Reduction guide created on 25/10/2017

.. include:: globals.rst

.. |fig-1| replace:: :numref:`fig-1`

Reducing |hiper| data
*********************

This guide assumes that you have got started with the |hiper| pipeline (see
:doc:`telescope` for a quick start to observing), but now would like a more
detailed guide to reducing your data.  It covers the following steps:

.. contents::
   :local:

.. Note::

  I am aware that these pages need more illustrations. This is work in
  progress.


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

Once you have a bias frame, then it can be used by editing in its name in the
calibration section of the reduce file.

.. Note::
   One should take bias calibration images in all the output formats used
   to guard against small changes that can occur as one changes output
   window formats. We have not yet established the significance of this
   for |hiper|.

.. Warning::
   Do not take bias frames too soon after (within less than 20 minutes)
   powering on the CCDs to avoid higher than normal dark current. |makebias|
   include a plot option to check this. Make sure to look at this if the bias
   is taken not long after a power on.

Flat fields
===========

Flat fields with |hiper| are taken at twilight in areas of sky with as few
stars as possible. One should offset the telescope between images to allow
stars to be removed. At the GTC, |hiper|'s driving routine, ``hdriver`` can
drive the telescope as well as the instrument, making spiralling during sky
flats straightforward. One can normally aquire more than 100 frames in a
single run, but the different CCDs will have different count levels on any one
frame, and will come out of saturation at different times. The count levels
will also be falling or rising according to whether the flats were taken at
evening or morning twilight.

The task of making the flat fields is thus to combine a series of frames with
differing count levels while removing features that vary between images. In
order to do this, one must normalise the images by their mean levels, but
weight them appropriately in the final combination to avoid giving too much
weight to under-exposed images. This is tedious by hand, and therefore the
command |makeflat| was written to carry out all the necessary tasks.

As with the biases, it is strongly recommended that you inspect the frames to
be combined using |rtplot| to avoid including any disastrous ones. Saturated
frames are spotted using user-defined mean levels at which to reject
frames. The documentation of |makeflat| has details of how it works, and you
are referred to this for more information.

.. Warning::
   It is highly advisable to compute multiple versions of the flat field
   using different values of the parameter ``ngroup`` which can have a
   significant effect on the removal of stars from the final frame and
   then to compare the results against each other. See |makeflat|.

Fringing / darks
================

With |hiper| there is some fringing in the z-band and to a much
smaller extent in the i-band. Fringing does not flat-field away because it is
the result of illumination by the emission-line dominated night sky whereas
twilight flats come from broad-band illumination. It should be possible to
correct for it, but I have yet to do so. This is mostly here for warning.

.. Warning::
   I have yet to implement corrections for fringing or dark current in the
   |hiper| pipeline.

Aperture files
==============

The pipeline photometry provides straightforward aperture photometry. Many of
the details can be defined when setting the apertures using |setaper|. Not
only can you choose your targets, but you can mask nearby stars from the sky
aperture, and you can to a certain extent sculpt your target apertures which
can help with blended interlopers by including them in an over-sized aperure.

A key decision to be made at this stage is whether you think your target will
remain detectable on each frame throughout the run. Detectable means that it's
position can be measured and thus the relevant aperture re-positioned.  If
not, then |setaper| gives you the option of ``linking`` any target to another,
with the idea that a brighter target can define the position shifts which are
applied to the fainter target.

An example of a set of apertures showing all these features is shown in
|fig-1|.

.. _fig-1:

.. figure:: complex_mask.png
   :scale: 50 %
   :alt: Click on image to enlarge
   :align: center

   Example of a complex set of apertures. The target is at the centre of the
   circles in the lower-right. The comparison star is in the upper-left.
   [click on the image to enlarge it].

In this case the target star has two nearby companions which causes three
problems: (1) the sky annulus may include flux from the companions, (2) the
target aperture can include a variable contribution from the companions,
depending upon the seeing, and (3) it is hard to locate the object because the
position location can jump to the nearby objects. The set of apertures shown
in |fig-1| combats these problems as follows. First, there are pink/purplish
dashed circles connected to the centre of the apertures. These are *mask*
regions which exclude the circled regions from any consideration in sky
estimation. NB they *do not* exclude the pixels from inclusion in target
apertures; this is not possible without systematic bias without full-blown
profile fits. Second, are somewhat similar brown/green dashed circles. These
are *extra* apertures which indicate that the flux in the regions enclosed is
to be added to the flux in the target aperture. This offsets the problem of
variable amounts of nearby stars' flux being included in the
aperture. Finally, the thick pink arrow pointing from the lower-right (target)
aperture to the upper-left reference aperture (green circle) *links* the
target aperture. This means its position is calculated using a fixed offset
from the reference aperture. This is often useful for very faint targets, or
those, which like the one shown here, havd close-by objects that can confuse
the re-positioning code.

Reduction files
===============

Once you have defined the apertures for all CCDs that you can, you can create
a reduction file using |genred|. This reads which CCDs and apertures have been
set and tries to create a file that at will work with |reduce|, even if it may
well not be ideal. Once you have this file, then you should expect to go
through a sequence of running |reduce|, adjusting the reduce file, re-running
|reduce| until you have a useful plot going. As mentioned above, it is inside
the reduce file that you can set the name of your bias, and also the flat
field file. You should experiment with changing settings such as the
extraction lines and aperture re-position section to get a feel for the
different parameters.

The two keys elements to get right in the reduction files are the sections
concerning how the apertures are re-located and the extraction. The former is
the trickier one. The apertures are re-located in a two-step process. First a
search is made in a square box centred on the last position measured for a
given target. The size of this box (``search_half_width``) is important. If it
too small, your target may jump beyond it and the reduction is likely to fail
as that will set an error condition and the target may never be recovered. On
the other hand too large, and you could easily end up jumping to a nearby
object. (The search identifies the brightest target in the box after gaussian
smoothing to reduce the effects of cosmic rays.) This is where rerence targets
come in. A well chosen reference target (bright, in an isolated region), can
allow you to set a relatively large search box, robust against large
jumps. The shifts from the reference targets are used to place the profile fit
boxes and avoid the need to searches over non-reference targets. This can be a
great help on faint targets. The main decision on extraction are variable or
fixed apertures (usually recommend variable), optimal or normal extraction
(usually recommend normal), and the various aperture radii. It is recommended
to plot the images in reduce at least once zoomed in on your target to get a
feel for this. Depending on circumstances, significant improvments to the
photometry can be made with careful choices on these parameters; do not assume
that the file produced by |genred| is in any way the last word; there is no
automatic way to come up with the ideal choices which depend upon the nature
of the field and conditions.

Plotting results
================

|reduce| delivers a basic view of your data as it comes in, which is usually
enough at the telescope. If you want to look at particular features, then you
should investigate the command |plog|. This allows you to plot one parameter
versus another, including division by comparison stars. The |plog| code is a
good place to start from when analysing your data in more detail. In
particular it shows you how to load in the rather human-unreadable |hiper| log
files (huge numbers of columns and rows).

Customisation
=============

You may well find that your data has particular features that the current
pipeline does not allow for. An obvious one is with crowded fields, which
can only roughly be accomodated with judicious application of the options
within |setaper|. The pipeline does not aim to replicate packages designed to
handle crowded fields, and you are best advised to port the data over into
single frames using |grab|, remembering that the 'hcm' files are nothing more
than a particular form of FITS. If your data requires only a small level of
tweaking then there are a few simple aritematical commands such as |add| and
|cmul| that might help, but it is not the intention to provide a full suite of
tools that can deal with all cases. Instead, the recommended route is to code
Python scripts to manipulate your data, and the :doc:`api` is designed to make
it relatively easy to access |hiper| data. If you devise routines of generic
interest, you are encouraged to submit them for possible inclusion within the
main pipeline commands. The existing pipeline commands are a good place to
start when looking for examples.

The reduce log files
====================

|reduce| writes all results to an ASCII log file. This can be pretty enormous
with many columns. The log file is self-documenting with an extensive header
section which is worth a look through. In particular the columns are named and
given data types to aid ingestion into numpy recarrays. The pipeline command
|plog| provides a crude interface to plotting these files, and module
:mod:`hipercam.hlog` should allow you to develop scripts to access the data
and to make your own plots.

Trouble shooting reduction
==========================

There are several things you can do to avoid problems during reduction. The
main thing to avoid is that |reduce| simply loses your target or the
comparison stars.

Use reference apertures
-----------------------

If you identify a star (or stars) as "reference apertures", their position is
the first to be determined and then used to offset the location before
carrying out profile fits for other stars. If you choose a well-isolated
reference star, this can allow cope with large changes in position from frame
to frame, while maintaining a tight search on non-reference stars which may be
close to other objects and be difficult to locate using a more wide-ranging
search. Sensible use of this can avoid the need to link apertures in some
cases.

Defocussed images
-----------------

Defocussing is often used in exoplanet work. Defocussed images are not well
fit by either gaussian or moffat profiles. In this case, when measuring the
object position, you should hold the FWHM fixed and use a large FWHM,
comparable to the width of the image. Experiment for best results. You should
also raise the threshold for bad data rejection as well. The idea is simply to
get a sort of weighted centroid, and you will not get a good profile fit. For
very defocussed images, it is important to avoid too narrow a FWHM otherwise
you could end up zeroing in on randome peaks in the profile.

Linked apertures
----------------

Linked apertures can be very useful if your target is simply too faint or near
others to track well. However, they bring the danger, especially in long runs,
of atmospheric refraction shifting the relative positions of the linked
apertures, leading to loss of flux. This is particularly the case in the
u-band. I have not, but will at some point, implement a two-pass reduction to
get round this issue, but in the meantime, if you have to link an aperture,
try to do so with a nearby object to minimise the distortion. It does not need
to be super-bright (although preferably it should be brighter than your
target), or your main comparison; the key point should be that its position
can be securely tracked.

Experiment
----------

Try different settings. Especially the extraction settings can make a
significant difference to the results. Compare the results visually. There is
no one prescription that works for all cases. Faint targets normally require
different setting from bright ones for the best results.
