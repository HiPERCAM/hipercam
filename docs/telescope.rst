.. Observer's guide created on 25/10/2017

.. include:: globals.rst

.. highlightlang:: rest

At the telescope
****************

Here it is assumed that you are observing with |hiper|, ULTRACAM or
ULTRASPEC for the first time and want to look at your data. It is
written in tutorial style for someone starting afresh who just wants
to get going, so it is as short as possible, with details listed
elsewhere. The first thing to realise is that the pipeline is run
entirely from command lines entered in a terminal, so first of all
open a terminal on the "drpc" (data reduction PC).

Plotting the data coming in
===========================

Assuming that you are observing and data are being acquired to the
rack PC at the telescope, the first command to try is |rtplot| [#f0]_
. Try typing in a terminal as follows::

  rtplot

You should see something like [#f1]_ ::

  rtplot
  run - run name [run0068]:

which is prompting you to supply a name for the run. If this is your first
time, the suggested value 'run0068', which remains from the last time the
command was run, is likely to be of little use [#f2]_ . So you need to type the
name of your run, which might be 'run0002'. If you don't know it, either look
at the observing GUI to see the run number, or ctrl-C out of |rtplot| (a safe
option to quit any pipeline command) and type instead::

  hls

for |hiper|, or::

  uls

for ULTRACAM and ULTRASPEC (which works with a different server).
This returns a list of runs; you probably want the last one of the list, which
will be the most recent. Note that all |hiper| runs appear as 'run####', while
ULTRACAM and ULTRASPEC have just three digits. Returning to |rtplot| we can
enter more inputs::

  rtplot
  run - run name [run0068]: run0002
  first - first frame to plot [100]: 11
  ccd - CCD(s) to plot [0 for all] [3]: 2 3 4 5
  nx - number of panels in X [2]:
  bias - bias frame ['none' to ignore] [none]:
  defect - defect file ['none' to ignore] [defect]:
  setup - display current hdriver window settings [True]:
  msub - subtract median from each window? [True]:
  iset - set intensity a(utomatically), d(irectly) or with p(ercentiles)? [p]:
  plo - lower intensity limit percentile [10.0]:
  phi - upper intensity limit percentile [99.0]:
  xlo - left-hand X value [1360.0]: min
  xhi - right-hand X value [1400.0]: max
  ylo - lower Y value [280.0]: min
  yhi - upper Y value [310.0]: max

Each line consists of an input parameter name (e.g. `bias`, `xlo`), followed
by a dash, a prompt, a default value in square brackets then a colon. In some
lines the default has been accepted (e.g. `plo`) by hitting <enter>,
whereas in others, something new has been entered, e.g. `ccd` where the
default value '3' has been replaced with '2 3 4 5'. (for ULTRACAM you would
want some combination of '1', '2' and '3'; ULTRASPEC has a single CCD and
you will not be prompted).

If we wanted to change one parameter, say `first`, the first frame to plot,
we can type::

  rtplot
  run - run name [run0002]:
  first - first frame to plot [11]: 12
  ccd - CCD(s) to plot [0 for all] [3]: \

The backslash '\\\\' means adopt the default values for all remaining
parameters. An alternative with the same effect is::

  rtplot first=12 \\

where the '\\\\\\\\' (double backslash) is needed because the shell
will strip one off. A simple ``rtplot \\`` at the terminal will repeat
the command with no changes of any parameters. The ability to specify
parameters on the command line is particularly useful for the later
ones e.g.::

  rtplot xlo=100 xhi=1100 \\

will reduce the X-range without the need to wind through the initial
parameters. The one other thing to be aware of are "hidden" parameters. For
instance::

  rtplot source=hs

ensures that |rtplot| tries to get the data from the |hiper| server on
the rack (which must itself be running of course), while you would
want to substitute 'us' when using the ULTRACAM|SPEC server, or 'ul'
when accessing ULTRACAM|SPEC data locally. One very useful feature of
|rtplot| is for quick examination of the FWHM of images through a
parameter called ``profit`` that will only be prompted if you display
just a single CCD.  This is useful to see if your target has
reasonable count levels, and as a quick check on focus, although
focussing can only be done properly by looking at the FWHMs of all
CCDs during reduction.  The section on :ref:`command-calling` has many
further details on the parameter input system and should be read the first
time you use the software. From now it will be assumed that you know the
parameter input system.

.. Note::

  |rtplot| is one of the most useful commands during observing runs. Three
  useful features are (1) the ability to fit the profiles of selected targets
  when displaying a single CCD; (2) the option to plot a file of CCD defects
  (:download:`hipercam defect file <hipercam_defects_2018_05_21.dft>`;
  :download:`ultracam defect file <ultracam_defects_2018_06_02.dft>`). Use
  this to avoid bad regions of the CCDs; (3) the `setup` option which
  interrogates the `hdriver` GUI for the current set of windows. This is
  useful for optimising the CCD setup on a target.

Commands used in this section: |rtplot|, |hls|

Looking in more detail
======================

To look in more detail at particular images, you can save individual
exposures to disk using the command |grab| as e.g.::

  grab
  run - run name [run0002]:
  ndigit - number of digits in frame identifier [3]:
  first - first frame to grab [1]:
  last - last frame to grab [2]: 5
  bias - bias frame ['none' to ignore] [none]:
  Written frame 1 to run0002_001.hcm
  Written frame 2 to run0002_002.hcm
  Written frame 3 to run0002_003.hcm
  Written frame 4 to run0002_004.hcm
  Written frame 5 to run0002_005.hcm

In this case the first 5 exposures have been saved to "hcm" files. These
are in fact FITS files but the extension is to indicate that they have
a format specific to |hiper|. For more on the format of these files consult
:doc:`hcmformat`.

You can then look at these with the |hiper| command |hplot| which, if you use
the matplotlib option (``device=/mpl``), can allow you to zoom in on the same
part of all CCDs simultaneously and also enables some simple statistics. These
frames can also be displayed with other packages such as
`ds9 <http://ds9.si.edu/site/Home.html>`_ or
`gaia <http://star-www.dur.ac.uk/~pdraper/gaia/gaia.html>`_.

Commands mentioned in this section: |grab|, |hplot|.

Plotting defects during target acquisition
==========================================

It is **strongly** recommended that you plot CCD defects during
acquisition.  For ULTRACAM|SPEC we have created defect files ('.dft')
for this purpose from flat-field and bias frames. Flat-field defects
in particular can often be hard to see, during acquistion, especially
if the sky is dark. Plotting them can help you avoid the worst defects
which you only later discover during reduction, when it is too late.

Getting reduction going
=======================

Typically one of the first things you want to do during a run is to
start a reduction which allows you to monitor the conditions, above
all the transparency and seeing (and therefore the telescope
focus). For multi-CCD images it is important to measure the FWHM in
all CCDs because if the cameras are not par-focal, you may have good
images in one but not another, but even with just one CCD as with
ULTRASPEC, one of the first things you should do is check the focus.

The |hiper| reduction is essentially a three-step process:

 1. Define target apertures with the |setaper| command.

 2. Generate a text file ("reduce file") defining the reduction
    with the |genred| command.

 3. Run the |reduce| command to see lightcurves.

Setting the apertures
---------------------

Beginning with |setaper|, one selects the targets to be extracted using a
cursor over an image. For the image, you can simply extract a single frame
from the run using |grab|, but it is better practice to average a few frames
which is easily done with |averun|::

  averun run0076 1 10 none run0076

which averages (actually take the median in this case) the first 10 frames
without bias subtraction, and stores the result in :file:`run0076.hcm`. Next
run |setaper| on this image. This should display the image and allow you to
add apertures one by one. Standard |hiper| practice is to label the main
target '1', the main comparison star '2', and then number them from then on in
order of increasing faintness. See the docs on |setaper| for the many features
of the program. The end result of |setaper| is an "aperture file" with
extension `.ape`, which it is recommended you call by the same name as the run
itself, and thus in this case you would end up with a file called
:file:`run0076.ape`. This is a JSON-format text file, and can be edited
directly, although it is much safer to let |setaper| do the job for you.

Commands mentioned in this section: |averun|, |grab|, |setaper|.

Generating a reduce file
------------------------

Once you have defined the apertures, you are ready to create a "reduce file",
which is a text file with a series of directives specifying how the reduction
is to be carried out. The easiest way to do this is through the command
|genred| which generates the file given the aperture file you just created. In
this case the command would be::

  genred run0076 run0076 "A test" none none none false

which generates a reduce file called :file:`run0076.red` with no bias,
flat or dark frame ('none none none') and using a magnitude scale
(linear='false') at the end. It is likely that the resultant reduce
file will not be right for you, and you will want to adjust it, but
the aim is to provide a start to get you going quickly. The reduce
file is just a text file and is easy to edit.

.. Note::

   |genred| has a **huge** number of hidden parameters. The very first
   time you run it for a given observing run, specify "prompt" on the
   command line so that these hidden parameters are prompted. Once you
   have done this once, they can usually be left the same value from
   one run to the next. This way you only need modify a few obvious
   parameters for subsequent invocations of |genred|.

Commands mentioned in this section: |genred|.

Reducing your data
------------------

You are now ready to run the |reduce|. This should generate a plot in
which you can see light curves, the object X and Y position, the
transmission and the FWHM of the images. You will almost certainly
need to adjust the plots, especially the light curve plots, by adding
offsets to the ``plot`` lines in the ``[light]`` section of the reduce
file. Apart from this, the first thing to do is to look at the FWHM
plots (at the bottom of the light curve plot panel) and see whether
the telescope focus can be adjusted to improve them. It is also worth
displaying images at least once (parameter ``implot``) in order to see
that things are doing what you expected, but, once you are going, you
may want to switch off this option for speed.

If you have got this far, well done. You will see some steps, such as bias
subtraction have been skipped. This was in the spirit of getting going fast.
Look at :doc:`reduction` for a longer and slower look at how to get more
out of the reduction.

Commands mentioned in this section: |reduce|.

.. rubric:: Footnotes

.. [#f0] All |hiper| commands are implemented as entry points to functions,
         all of which are part of the :mod:`hipercam.scripts` module. This
         means that documentation on any |hiper| pipeline command can be
         looked up from a terminal with `pydoc`, e.g. for |rtplot| you would
         type ``pydoc hipercam.scripts.rtplot``. To invoke the same command,
         simply type ``rtplot`` in a terminal.


.. [#f1] This may fail if the data source is set incorrectly. If so, re-type
         as``rtplot prompt`` to see all the options, many of which are hidden
         by default. The one you might want to change is ``source`` which
         should be ``hs`` for data from the |hiper| server at the telescope,
         or ``hl`` for data from a disk file local to your computer (or ``us``
         and ``ul`` for ULTRACAM).

.. [#f2] All raw |hiper| file names take the form 'run####.fits' but
         the extension '.fits' can be omitted. ULTRACAM|SPEC files take the
         form 'run###.dat' / 'run###.xml'. Note that |hiper| files have
         4 digits compared to ULTRACAM|SPEC's 3.
