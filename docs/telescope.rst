.. Observer's guide created on 25/10/2017

.. include:: globals.rst

.. highlightlang:: rest

At the telescope
****************

Here it is assumed that you are observing with |hiper| for the first time and
want to look at your data. It is written in tutorial style for someone
starting afresh who just wants to get going, so it is as short as
possible. The first thing to realise is that the pipeline is run entirely from
command lines entered in a terminal, so first of all open a terminal on the
"drpc" (data reduction PC).

Plotting the data coming in
===========================

The first command to try is |rtplot| [#f0]_ . Try typing in a terminal as
follows::

  rtplot

You should see something like [#f1]_ ::

  rtplot
  run - run name [run0068]:

which is prompting you to supply a name for the run. If this is your first
time, the suggested value 'run0068', which remains from the last time the
command was run, is likely to be of little use [#f2]_ . So you need to type the
name of your run, which might be 'run0002'. If you don't know it, either look
at the 'hdriver' observing GUI to see the run number, or ctrl-C out of
|rtplot| (a safe option to quit any pipeline command) and type instead::

  hls

This returns a list of runs; you probably want the last one of the list, which
will be the most recent. Note that all |hiper| runs appear as 'run####'.
Returning to |rtplot| we can enter more inputs::

  rtplot
  run - run name [run0068]: run0002
  first - first frame to plot [100]: 11
  ccd - CCD(s) to plot [0 for all] [3]: 2 3 4 5
  nx - number of panels in X [2]:
  bias - bias frame ['none' to ignore] [none]:
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
default value '3' has been replaced with '2 3 4 5'.

If we wanted to change one parameter, say `first`, the first frame to plot,
we can type::

  rtplot
  run - run name [run0002]:
  first - first frame to plot [11]: 12
  ccd - CCD(s) to plot [0 for all] [3]: \

The backslash '\\\\' means adopt the default values for all remaining
parameters. An alternative with the same effect is::

  rtplot first=12 \\

where the '\\\\\\\\' (double backslash) is needed because the shell will
strip one off. A simple ``rtplot \\`` will repeat the command with
no changes of any parameters. The ability to specify parameters on the
command line is particularly useful for the later ones e.g.::

  rtplot xlo=100 xhi=1100 \\

will reduce the X-range without the need to type through the initial
parameters. See the :ref:`command-calling` section for more details on the
parameter input system. From now it will be assumed that you know the parameter
input system.

|rtplot| can also be used to examine the FWHM of images through a hidden
parameter called ``profit`` (see :ref:`command-calling` for more on hidden
parameters). This is useful also to see if your target has reasonable count
levels.

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

Commands used in this section: |grab|, |hplot|

Getting reduction going
=======================

Typically one of the first things you want to do during a run is to start
a reduction which allows you to monitor the conditions, above all the
transparency and seeing (and telescope focus). For multi-CCD images it is
important to measure the FWHM in all CCDs because if the cameras are not
par-focal, you may have good images in one but not another.

The |hiper| reduction is a three-step process:

 1. Define target apertures with the |setaper| command.

 2. Generate a text file ("reduce file") defining the reduction
    with the |genred| command.

 3. Run the |reduce| command to see lightcurves.

Setting the apertures
---------------------

Beginning with |setaper|, one selects the targets to be extracted using a
cursor over an image. You can simply select a single frame from the run using
|grab|, but it is better practice to average a few frames which is easily 
done with |averun|::

  averun run0076 1 10 none run0076

which averages the first 10 frames without bias subtraction, and stores the
result in :file:`run0076.hcm`. Next run |averun| on this image. This should
display the image and allow you to add apertures one by one. Standard |hiper|
practice is to label the main target '1', the main comparison star '2', and
then number them from then on in order of increasing faintness. See the docs
on |averun| for the many features of the program. The end result of |averun|
is an "aperture file", which it is recommended you call by the same name as
the run itself, which in this case would produce a file called 
:file:`run0076.ape`.

Generating a reduce file
------------------------

Once you have defined the apertures, you are ready to create a "reduce file",
which is a text file with a series of directives specifying how the reduction
is to be carried out. The easiest way to do this is through the command
|genred| which generate the file given the aperture file just created. In this
case the command would be::

  genred run0076 run0076 "A test" none none none false

which generates a reduce file called :file:`run0076.red` with no bias, flat or
dark frame ('none none none') and using a magnitude scale (linear='false')
at the end. It is very likely, almost certain in fact, that the resultant
reduce file will not be right, and you will want to adjust it. This is easy
since it is just a text file and can easily be edited.

Running |reduce|
----------------

Finally you are ready to run the |reduce|. This should generate a plot in
which you can see light curves, the object position, the transmission and the
FWHM of the images. You will almost certainly need to adjust the plots,
especially the light curve plots, by adding offsets to the ``plot`` lines in
the ``[light]`` section of the reduce file. Apart from this, the first thing
to do is to look at the FWHM plots and see whether the focus can be adjusted
to improve them. It is also worth displaying images at least once in order to
see that things are behaving as you expect.

If you have got this far, well done. You will see some steps, such as bias
subtraction have been skipped. This was in the spirit of getting going fast.
Look at :doc:`reduction` for a longer and slower look at how to get more
out of reduction.

.. rubric:: Footnotes

.. [#f0] All |hiper| commands are implemented as entry points to functions,
         all of which are part of the :mod:`hipercam.scripts` module, hence
         the slightly odd-looking command name. The same string, without the
         closing pair of braces, is needed when looking up help from the
         terminal using ``pydoc``.

.. [#f1] This may fail if the data source is set incorrectly. If so, re-type
         as "rtplot PROMPT" to see all the options, many of which are hidden
         by default. 

.. [#f2] All raw |hiper| file names take the form 'run####.fits' but
         the extension '.fits' can be omitted.
