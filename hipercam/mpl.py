# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
matplotlib plotting functions.
"""

# Standard pre-amble from astropy
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.extern import six
from .window import *

def pwin(axes, win, col='k', lw=1, **kwargs):
    """
    Plots boundary of a :class:`Window` as a line. (matplotlib)
    Arguments::

      axes : (:class:`matplotlib.axes.Axes`)
           the Axes to plot to.

      win : (Window)
           the :class:`Window` to plot

      col : (matplotlib colour)
           the colour to use

      lw : (float)
           linewidth in points

      kwargs : (keyword arguments)
           other arguments to feed to :func:`matplotlib.pyplot.plot`
    """
    left,right,bottom,top = win.extent()
    axes.plot([left,right,right,left,left],[bottom,bottom,top,top,bottom],
              c=col,lw=lw,**kwargs)

def pwind(axes, wind, col='k', lw=1, aspect='equal', **kwargs):
    """Plots :class:`Windata` as an image with a line border. (matplotlib).
    Note that the keyword arguments are only passed to :func:`imshow` and
    you should plot the border separately if you want anything out of the
    ordinary. If no colour map ('cmap') is specified, it will be set to
    'Greys'

    Arguments::

      axes : (:class:`matplotlib.axes.Axes`)
           the Axes to plot to.

      wind : (Windata)
           the :class:`Windata` to plot

      col : (matplotlib colour)
           the colour to use for the border

      lw : (float)
           linewidth of the border in points (0 to skip)

      kwargs : (keyword arguments)
           other arguments to feed to :func:`matplotlib.pyplot.imshow`.
           Useful ones are 'vmin', 'vmax', 'cmap'

    """
    left,right,bottom,top = wind.extent()

    if 'cmap' in kwargs:
        axes.imshow(wind.data,extent=(left,right,bottom,top),
                    aspect=aspect,origin='lower',
                    interpolation='nearest',**kwargs)
    else:
        axes.imshow(wind.data,extent=(left,right,bottom,top),
                    aspect=aspect,origin='lower',cmap='Greys',
                    interpolation='nearest',**kwargs)
    pwin(axes, wind, col, lw)
