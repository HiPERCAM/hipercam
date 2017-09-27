# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
PGPLOT plotting functions.

PGPLOT is a venerable plotting package that is still good when it comes to
speed hence it has a place in the hipercam package although I anticipate
replacing it in the long term.
"""

from trm.pgplot import *
from .core import *
from .window import *
from .ccd import *

__all__ = (
    'Params', 'Device',
    'pWin', 'pWind', 'pCcd'
    )


# some look-and-feel globals.
Params = {

    # colour indices rgb values
    'cis' : (
        (1,1,1),      # 0
        (0,0,0),      # 1
        (0.5,0,0),    # 2
        (0,0.5,0),    # 3
        (0,0,0.5),    # 4
        (0.5,0.5,0),  # 5
        (0.5,0,0.5),  # 6
        (0,0.5,0.5),  # 7
    ),

    # axis label character height
    'axis.label.ch' : 1.5,

    # axis number character height
    'axis.number.ch' : 1.5,

    # axis colour index
    'axis.ci' : 4,

    # window box colour index
    'ccd.box.ci'   : 8,

    # window label character height
    'win.label.ch' : 1.2,

    # window label colour index
    'win.label.ci' : 6,

    # window box colour index
    'win.box.ci'   : 7,
}

class Device(PGdevice):
    """Sub-class of PGdevice that after opening the plot device, re-defines colour indices
    according to the values set in pgp.Params

    """

    def __init__(self, device):
        super(Device, self).__init__(device)

        for i, (r,g,b) in enumerate(Params['cis']):
            pgscr(i,r,g,b)

def pWin(win, label=None):
    """
    Plots boundary of a :class:`Window` as a line. (PGPLOT)
    Plots to the current device with current line width, colour
    etc. This changes the fill-area style and the colour index.

    Arguments::

      win : (Window)
           the :class:`Window` to plot

      label : (string / None)
           label to plot at lower-left corner of Window
    """
    left,right,bottom,top = win.extent()

    pgsfs(2)
    pgsci(Params['win.box.ci'])
    pgrect(left,right,bottom,top)

    if label is not None:
        pgsci(Params['win.label.ci'])
        pgsch(Params['win.label.ch'])
        pgptxt(left,bottom,0,1.3,label)

def pWind(wind, vmin, vmax, label=None):
    """Plots :class:`Windata` as an image with a line border. (PGPLOT).

    Arguments::

      wind : (Windata)
           the :class:`Windata` to plot

      vmin : (float)
           value at minimum of scale

      vmax : (float)
           value at maximum of scale

      label : (string / None)
           label to plot at lower-left corner of Window
    """
    tr = [wind.llx+(wind.xbin-1)/2,wind.xbin,0,
          wind.lly+(wind.ybin-1)/2,0,wind.ybin]
    pggray(wind.data, vmax, vmin, tr)
    pWin(wind, label)

def pCcd(ccd, iset='p', plo=5., phi=95., dlo=0., dhi=1000.):
    """Plots :class:`CCD` as a set of :class:`Windata` objects correctly
    positioned with respect to each other.

    Arguments::

      ccd : (CCD)
           the :class:`CCD` to plot

      iset : (string)
           how to set the intensity scale to be used. 'p' for percentiles (set
           using plo and phi); 'a' automatic for min/max range; 'd' for direct
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

    elif iset == 'a':
        # Set intensities from min/max range
        vmin, vmax = ccd.min(), ccd.max()

    elif iset == 'd':
        vmin, vmax = dlo, dhi

    else:
        raise ValueError(
            'did not recognise iset = "' + iset + '"'
        )

    for key, wind in ccd.items():
        pWind(wind, vmin, vmax, '{!s}'.format(key))

    # plot outermost border of CCD
    pgsci(Params['ccd.box.ci'])
    pgrect(0.5,ccd.nxtot+0.5,0.5,ccd.nytot+0.5)

    return (vmin,vmax)
