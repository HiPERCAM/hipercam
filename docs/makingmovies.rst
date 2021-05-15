.. how to make movies doc created on 03/02/2018

.. include:: globals.rst

.. highlightlang:: rest

Making movies from |hiper| data
*******************************

The |makemovie| command allows you to generate stills from a movie. It
is however painful to use and even more to actually compile, so this
is just to provide an easy reminder.

Making the stills
=================

TBD

Making the movie from the stills
================================

I use 'ffmpeg' to make the movie with a command like::

  ``ffmpeg -start_number 123 -i run0005_%04d.png -c:v libx264 -pix_fmt yuv420p movie.mp4 -y``

which would find any files of the form 'run0005_0340.png' starting from 'run0005_0123.png'.




