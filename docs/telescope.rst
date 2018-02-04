.. Observer's guide created on 25/10/2017

At the telescope
****************

Here it is assumed that you are observing with HiPERCAM for the first time and
want to look at your data. It is written in tutorial style for someone
starting afresh who just wants to get going. The first thing to realise is
that the pipeline is run from the command line entered in a terminal. 

Plotting the data coming in
===========================

The first command to try is 'rtplot'. Try typing in a terminal as follows::

  rtplot

You should see something like [#f1]_ ::

  rtplot
  run - run name [run0068]:

which is prompting you to supply a name for the run. If this is your first
time, the suggested value 'run0068', which remains from the last time the
command was run, is likely to be of little use [#f2]_ . So you need to type the
name of your run, which might be 'run0002'. If you don't know it, either look
at the 'hdriver' observing GUI to see the run number, or ctrl-C out of
'rtplot' (a safe option to quit any pipeline command) and type instead::

  hls

This returns a list of runs; you probably want the last one of the list, which
will be the most recent. Note that all HiPERCAM runs appear as 'run####'.
Returning to 'rtplot' we can enter more inputs::

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

Each line consists of an input parameter name (e.g. 'bias', 'xlo'), followed
by a dash, a prompt, a default value in square brackets then a colon. In some
lines the default has been accepted (e.g. 'plo') by hitting the <enter> key,
whereas in others, something new has been entered, e.g. 'ccd' where the
default value '3' has been replaced with '2 3 4 5'.

If we wanted to change one parameter, say 'first', the first frame to plot,
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
parameters.

Looking in more detail
======================

You can grab individual images using `grab`::

  grab
  run -



.. rubric:: Footnotes

.. [#f1] This may fail if the data source is set incorrectly. If so, re-type
         as "rtplot PROMPT" to see all the options, many of which are hidden
         by default. 

.. [#f2] All raw HiPERCAM file names take the form 'run####.fits' but
         the extension '.fits' can be omitted.
