.. how to make movies doc created on 03/02/2018

.. include:: globals.rst

.. highlight:: rest

Making |hiper| movies
*********************

The |makemovie| command allows you to generate stills from a movie. It
is a little painful to run and even more so to actually compile the
stills into a movie, so this is just to provide an easy reminder.

Making the stills
=================

|makemovie| can plot CCD images only or images plus a light curve if
you have a |reduce| log file available. I would suggest that you
either run it within an empty directory or create such a directory for
directing the output since it will create one image per
exposure. Start by creating just a few images and look at them with
your favourite PNG viewer. You will need to spend a while tweaking
things. The one thing that will change when you switch to a full run
is the X-scale on the light curve plot since that is adjusted to fit
the times of the images shown. You can adjust the apparent font size
by altering the plot width and height, at the same time as the dpi to
control the resolution. Crude, but it works. I haven't got it working
without rather a stack of "whitespace"; something for the future.

Here is a full listing of parameters I used to make a movie of the
eclipsing double white dwarf ZTF J1539+5027. In this case I was in the
directory containing the reduce log within which I had created an
empty directory "tmp". One up from the working directory was where the
raw data, run0005.fits sat.

.. code-block:: console

  makemovie
  source = hl
  run = ../run0005
  first = 2
  last = 1000
  trim = False
  ccd = 2
  bias = bias3x3_slow.hcm
  dark = none
  flat = none
  fmap = none
  defect = none
  log = run0005.log
  targ = 1
  comp = 2
  ymin = 0.0
  ymax = 0.05
  ynorm = [1.0]
  yoffset = [0.0]
  location = e
  fraction = 0.55
  lpad = (0.15, 0.2, 0.1, 0.15)
  cmap = viridis
  width = 12.0
  height = 3.0
  dstore = tmp
  ndigit = 4
  fext = png
  msub = False
  iset = d
  ilo = [0.0]
  ihi = [100.0]
  xlo = 0.0
  xhi = 2048.0
  ylo = 0.0
  yhi = 1024.0
  dpi = 200
  style = dots
  ms = 2.0


Making the movie from the stills
================================

Once you have the stills, you still need to actually make the movie. I use the command 'ffmpeg' to do so like this:

  ``ffmpeg -start_number 123 -i run0005_%04d.png -c:v libx264 -pix_fmt yuv420p movie.mp4 -y``

which would find any files of the form 'run0005_0340.png' starting
from 'run0005_0123.png'. 'ffmpeg' has a million-and-one options and
there be ways to tweak it, but this does a reasonable job I find.
You can loop the movie with something like:

  ``ffmpeg -start_number 3 -i run0011_%04d.png -filter_complex loop=5:362 -c:v libx264 -pix_fmt yuv420p movie.mp4 -y``

which would add 6 loops.
