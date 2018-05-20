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
from . import utils
from . import defect

__all__ = (
    'Params', 'Device',
    'pWin', 'pWind', 'pCcd'
    )


# some look-and-feel globals.
Params = {

    # axis label character height
    'axis.label.ch' : 1.2,

    # axis label colour
    'axis.label.ci' : 2,

    # axis number character height
    'axis.number.ch' : 1.,

    # axis colour index
    'axis.ci' : 4,

    # window box colour index
    'ccd.box.ci'   : 7,

    # window label character height
    'win.label.ch' : 1.2,

    # window label colour index
    'win.label.ci' : 5,

    # window box colour index
    'win.box.ci'   : 7,

    # aperture target colour
    'aper.target.ci' : 15,

    # aperture reference target colour
    'aper.reference.ci' : 3,

    # aperture sky colour
    'aper.sky.ci' : 2,

    # aperture label colour
    'aper.label.ci' : 6,

    # aperture link colour
    'aper.link.ci' : 5,

    # aperture mask colour
    'aper.mask.ci' : 5,

    # aperture extra colour
    'aper.extra.ci' : 14,

    # moderate defect colour
    'defect.moderate.ci' : 15,

    # severe defect colour
    'defect.severe.ci' : 2,

}

class Device(PGdevice):
    """Sub-class of PGdevice that after opening the plot device, re-defines colour
    indices according to a standardised set in core.CIS

    """

    def __init__(self, device):
        super(Device, self).__init__(device)

        for i, (r,g,b) in enumerate(CIS):
            pgscr(i,r,g,b)

def pWin(win, label=''):
    """
    Plots boundary of a :class:`Winhead` as a line. (PGPLOT)
    Plots to the current device with current line width, colour
    etc. This changes the fill-area style and the colour index.

    Arguments::

      win : :class:`Winhead`
           the :class:`Winhead` to plot

      label : string
           label to plot at lower-left corner of Winhead
    """
    left,right,bottom,top = win.extent()

    pgsfs(2)
    pgsci(Params['win.box.ci'])
    pgrect(left,right,bottom,top)

    if label != '':
        pgsci(Params['win.label.ci'])
        pgsch(Params['win.label.ch'])
        pgptxt(left,bottom,0,1.3,label)

def pWind(wind, vmin, vmax, label=''):
    """Plots :class:`Window` as an image with a line border. (PGPLOT).

    Arguments::

      wind : Window
           the :class:`Window` to plot

      vmin : float
           value at minimum of scale

      vmax : float
           value at maximum of scale

      label : string
           label to plot at lower-left corner of Winhead
    """
    tr = [wind.llx+(wind.xbin-1)/2,wind.xbin,0,
          wind.lly+(wind.ybin-1)/2,0,wind.ybin]
    pggray(wind.data, vmax, vmin, tr)
    pWin(wind, label)

def pCcd(ccd, iset='p', plo=5., phi=95., dlo=0., dhi=1000., tlabel=''):
    """Plots :class:`CCD` as a set of :class:`Window` objects correctly
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

      tlabel : (string)
           label to use for top of plot; 'X' and 'Y' will also be added if
           this is not blank.

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

    if tlabel != '':
        pgsch(Params['axis.label.ch'])
        pgsci(Params['axis.label.ci'])
        pglab('X','Y',tlabel)

    for key, wind in ccd.items():
        pWind(wind, vmin, vmax, '{!s}'.format(key))

    # plot outermost border of CCD
    pgsci(Params['ccd.box.ci'])
    pgrect(0.5,ccd.nxtot+0.5,0.5,ccd.nytot+0.5)

    return (vmin,vmax)

def pAper(aper, label='', ccdAper=None):
    """
    Plots an :class:`Aperture` object to the current PGPLOT device

    Arguments::

      aper    : (Aperture)
           the :class:`Aperture` to plot

      label   : (string)
           a string label to add

      ccdAper : (CcdAper)
           needed if plotting multiple apertures with links

    """

    # draw circles to represent the aperture. 'objs' is a list of the
    # objects that we keep to return for possible later deletion.
    pgsfs(2)
    pgsls(2)
    pgsls(1)
    if aper.ref:
        pgsci(Params['aper.reference.ci'])
    else:
        pgsci(Params['aper.target.ci'])
    pgcirc(aper.x,aper.y,aper.rtarg)

    pgsci(Params['aper.sky.ci'])
    pgcirc(aper.x,aper.y,aper.rsky1)
    pgcirc(aper.x,aper.y,aper.rsky2)

    if aper.link != '':
        # indicate a link with an arrow
        if ccdAper is None:
            raise ValueError(
                'to plot a linked aperture, need to pass through an CcdAper')
        else:
            laper = ccdAper[aper.link]

            # draw arrow starting at target aperture of one
            # aperture to the other.
            p1 = utils.Vec2D(aper.x, aper.y)
            p2 = utils.Vec2D(laper.x, laper.y)
            v  = p2-p1
            uv = v.unit()
            r1 = aper.rtarg*uv
            r2 = laper.rtarg*uv
            v -= r1+r2
            p1 += r1
            pgsci(Params['aper.link.ci'])
            pgarro(p1.x, p1.y, p1.x+v.x, p1.y+v.y)

    # draw dashed lines connecting the aperture to the centres of mask
    # indicated with circles.
    for xoff,yoff,r in aper.mask:
        # draw the line
        pgsci(Params['aper.mask.ci'])
        pgsls(2)
        pgline([aper.x,aper.x+xoff],[aper.y,aper.y+yoff])

        # and now the circle
        pgcirc(aper.x+xoff,aper.y+yoff,r)

    # draw dashed lines connecting the aperture to the centres of extra
    # indicated with circles.
    for xoff,yoff in aper.extra:
        # draw the line
        pgsci(Params['aper.extra.ci'])
        pgsls(2)
        pgline([aper.x,aper.x+xoff],[aper.y,aper.y+yoff])

        # and now the circle
        pgcirc(aper.x+xoff,aper.y+yoff,aper.rtarg)

    if label != '':
        pgsci(Params['aper.label.ci'])
        pgptxt(
            aper.x-aper.rsky2, aper.y-aper.rsky2,
            0., 1., label
        )

def pCcdAper(ccdaper):
    """
    Plots a :class:`CcdAper` object

    Arguments::

      ccdaper : (CcdAper)
           the :class:`CcdAper` to plot

    """
    for key, aper in ccdaper.items():
        pAper(aper, key, ccdaper)


def pDefect(dfct):
    """y
    Plots a :class:`Defect` o the current PGPLOT device

    Arguments::

      dfct    : Defect
           the :class:`Defect` to plot

    """

    if isinstance(dfct, defect.Point):

        if dfct.severity == defect.Severity.MODERATE:
            pgsci(Params['defect.moderate.ci'])
            pgsch(1.3)
            pgpt1(dfct.x, dfct.y, 17)

        elif dfct.severity == defect.Severity.SEVERE:
            pgsci(Params['defect.severe.ci'])
            pgsch(1.7)
            pgpt1(dfct.x, dfct.y, 17)

def pCcdDefect(ccd_dfct):
    """
    Plots a :class:`CcdDefect` object

    Arguments::

      ccd_dfct  : CcdDefect
           the :class:`CcdDefect` to plot

    """
    for key, dfct in ccd_dfct.items():
        pDefect(dfct)
