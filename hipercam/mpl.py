# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
matplotlib plotting functions.
"""

from .core import *
from .window import *
from .ccd import *

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
              color=col, lw=lw, **kwargs)

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

      aspect : (string)
           aspect ratio

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

def pccd(axes, ccd, iset='p', plo=5., phi=95., dlo=0., dhi=1000.,
         label=True, col='k', lw=1, aspect='equal', **kwargs):
    """Plots :class:`CCD` as a set of :class:`Windata` objects correctly
    positioned with respect to each other. The keyword arguments
    are passed to :func:`imshow`. If no colour map
    ('cmap') is specified, it will be set to 'Greys'

    Arguments::

      axes : (:class:`matplotlib.axes.Axes`)
           the Axes to plot to.

      ccd : (CCD)
           the :class:`CCD` to plot

      iset : (string)
           how to set the intensity scale to be used. 'p' for percentiles
           (set using plo and phi); 'r' for min/max range; 'd' for direct
           value set using dlo and dhi.

      plo : (float)
           lower percentile limit to use (if iset='p')

      phi : (float)
           upper percentile limit to use (if iset='p')

      dlo : (float)
           value to use for lower intensity limit (if iset='d')

      dhi : (float)
           value to use for upper intensity limit (if iset='d')

      label : (boolean)
           whether to attach a label to the Windows or not.

      col : (matplotlib colour)
           the colour to use for the borders

      lw : (float)
           linewidth of the borders in points (0 to skip)

      kwargs : (keyword arguments)
           other arguments to feed to :func:`matplotlib.pyplot.imshow`.
           'cmap' is useful. 'vmin' and 'vmax' are overridden by the
           iset method.

    Returns: (vmin,vmax) the intensity limits used.
    """
    if iset == 'p':
        # Set intensities from percentiles
        vmin,vmax = ccd.percentile((plo,phi))

    elif iset == 'r':
        # Set intensities from min/max range
        vmin, vmax = ccd.min(), ccd.max()

    elif iset == 'd':
        vmin, vmax = dlo, dhi

    else:
        raise HipercamError('mpl.pccd: did not recognise iset = "' +
                            iset + '"')
    kwargs['vmin'] = vmin
    kwargs['vmax'] = vmax

    for key, wind in ccd.items():
        pwind(axes, wind, col, lw, aspect, **kwargs)

    # plot outermost border of CCD
    axes.plot([0.5,ccd.nxtot+0.5,ccd.nxtot+0.5,0.5,0.5],
              [0.5,0.5,ccd.nytot+0.5,ccd.nytot+0.5,0.5],
              color=col, lw=lw)

    return (vmin,vmax)
