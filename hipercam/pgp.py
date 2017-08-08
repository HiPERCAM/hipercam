# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
PGPLOT plotting functions.
"""

from .core import *
from .window import *
from .ccd import *

def pwin(win):
    """
    Plots boundary of a :class:`Window` as a line. (PGPLOT)
    Plots to the current device with current line width, colour
    etc. This changes the fill-area style and the colour index.

    Arguments::

      win : (Window)
           the :class:`Window` to plot
    """
    left,right,bottom,top = win.extent()
    pgsfs(2)
    pgsci(6)
    pgrect(left,right,bottom,top)

def pwind(wind, vmin, vmax):
    """Plots :class:`Windata` as an image with a line border. (PGPLOT).

    Arguments::

      wind : (Windata)
           the :class:`Windata` to plot

      vmin : (float)
           value at minimum of scale

      vmax : (float)
           value at maximum of scale
    """
    tr = [wind.llx+(wind.xbin-1)/2,wind.xbin,0,
          wind.lly+(wind.ybin-1)/2,0,wind.ybin]
    pggray(wind.data, vmax, vmin, tr)
    pwin(wind)

def pccd(ccd, iset='p', plo=5., phi=95., dlo=0., dhi=1000.):
    """Plots :class:`CCD` as a set of :class:`Windata` objects correctly
    positioned with respect to each other.

    Arguments::

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

    Returns: (vmin,vmax), the intensity limits used.
    """
    if iset == 'p':
        # Set intensities from percentiles
        vmin, vmax = ccd.percentile((plo,phi))

    elif iset == 'r':
        # Set intensities from min/max range
        vmin, vmax = ccd.min(), ccd.max()

    elif iset == 'd':
        vmin, vmax = dlo, dhi

    else:
        raise ValueError(
            'pgp.pccd: did not recognise iset = "' + iset + '"'
        )

    for key, wind in ccd.items():
        pwind(wind, vmin, vmax)

    # plot outermost border of CCD
    pgsci(5)
    pgrect(0.5,ccd.nxtot+0.5,0.5,ccd.nytot+0.5])

    return (vmin,vmax)
