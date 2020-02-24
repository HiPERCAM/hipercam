# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
matplotlib plotting functions. 
"""

import os
from matplotlib.patches import Circle

from .core import *
from .window import *
from .group import *
from .ccd import *
from . import utils
from . import defect


# Font scale factor
SCALE_FACTOR = float(os.environ['HIPERCAM_MPL_FSCALE']) \
               if 'HIPERCAM_MPL_FSCALE' in os.environ else 1.0

# some look-and-feel globals.
Params = {

    # axis label fontsize
    'axis.label.fs' : 14*SCALE_FACTOR,

    # axis label colour
    'axis.label.col' : CIS[2],

    # axis number fontsize
    'axis.number.fs' : 14*SCALE_FACTOR,

    # axis colour
    'axis.col' : CIS[4],

    # window box colour
    'ccd.box.col' : CIS[1],

    # window label font size
    'win.label.fs' : 12*SCALE_FACTOR,

    # window label colour
    'win.label.col' : CIS[5],

    # window box colour
    'win.box.col' : CIS[7],

    # aperture target colour
    'aper.target.col' : CIS[1],

    # aperture reference target colour
    'aper.reference.col' : CIS[3],

    # aperture sky colour
    'aper.sky.col' : CIS[2],

    # aperture label colour
    'aper.label.col' : CIS[4],

    # aperture link colour
    'aper.link.col' : CIS[5],

    # aperture mask colour
    'aper.mask.col' : CIS[5],

    # aperture extra colour
    'aper.extra.col' : CIS[14],

    # moderate defect colour
    'defect.moderate.col' : CIS[15],

    # severe defect colour
    'defect.severe.col' : CIS[2],
}

def pWin(axes, win, label=''):
    """
    Plots boundary of a :class:`Winhead` as a line. (matplotlib)
    Arguments::

      axes   : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

      win    : :class:`Winhead`
           the :class:`Winhead` to plot

      label  : string
           label to attach to the Winhead

      kwargs : keyword arguments
           other arguments to feed to :func:`matplotlib.pyplot.plot`
    """
    left,right,bottom,top = win.extent()
    axes.plot([left,right,right,left,left],[bottom,bottom,top,top,bottom],
              color=Params['win.box.col'])
    if label != '':
        axes.text(left-3,bottom-3,label,fontsize=Params['win.label.fs'],
                  color=Params['win.label.col'], ha='right', va='top')

def pWind(axes, wind,  vmin, vmax, label=''):
    """Plots :class:`Window` as an image with a line border. (matplotlib).
    Note that the keyword arguments are only passed to :func:`imshow` and
    you should plot the border separately if you want anything out of the
    ordinary. If no colour map ('cmap') is specified, it will be set to
    'Greys'

    Arguments::

      axes : (:class:`matplotlib.axes.Axes`)
           the Axes to plot to.

      wind : (Window)
           the :class:`Window` to plot

      vmin : (float)
           image value at minimum intensity

      vmax : (float)
           image value at maximum intensity
    """
    left,right,bottom,top = wind.extent()

    # plot image
    axes.imshow(wind.data,extent=(left,right,bottom,top),
                aspect='equal',origin='lower',cmap='Greys',
                interpolation='nearest',vmin=vmin,vmax=vmax)

    # plot window border
    pWin(axes, wind, label)

def pCcd(axes, ccd, iset='p', plo=5., phi=95., dlo=0., dhi=1000., tlabel='',
         xlo=None, xhi=None, ylo=None, yhi=None):
    """Plots :class:`CCD` as a set of :class:`Window` objects correctly
    positioned with respect to each other.

    Arguments::

      axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

      ccd : CCD
           the :class:`CCD` to plot

      iset : string
           how to set the intensity scale to be used. 'p' for percentiles (set
           using plo and phi); 'a' for automatic min/max range; 'd' for direct
           value set using dlo and dhi.

      plo : float
           lower percentile limit to use (if iset='p')

      phi : float
           upper percentile limit to use (if iset='p')

      dlo : float
           value to use for lower intensity limit (if iset='d')

      dhi : float
           value to use for upper intensity limit (if iset='d')

      tlabel : string
           label to use for top of plot; 'X' and 'Y' will also be
           added if this is not blank.

      xlo : int | None
           left-hand limit to define region for computing percentile

      xhi : int | None
           right-hand limit to define region for computing percentile

      ylo : int | None
           lower limit to define region for computing percentile

      yhi : int | None
           upper limit to define region for computing percentile

    Returns: (vmin,vmax), the intensity limits used.

    """
    if iset == 'p':
        # Set intensities from percentiles
        vmin,vmax = ccd.percentile((plo,phi),xlo,xhi,ylo,yhi)

    elif iset == 'a':
        # Set intensities from min/max range
        vmin, vmax = ccd.min(), ccd.max()

    elif iset == 'd':
        vmin, vmax = dlo, dhi

    else:
        raise ValueError('did not recognise iset = "' + iset + '"')

    for key, wind in ccd.items():
        pWind(axes, wind, vmin, vmax, key)

    # plot outermost border of CCD
    axes.plot([0.5,ccd.nxtot+0.5,ccd.nxtot+0.5,0.5,0.5],
              [0.5,0.5,ccd.nytot+0.5,ccd.nytot+0.5,0.5],
              color=Params['ccd.box.col'])
    if tlabel != '':
        axes.set_title(tlabel,color=Params['axis.label.col'],fontsize=Params['axis.label.fs'])
        axes.set_xlabel('X',color=Params['axis.label.col'],fontsize=Params['axis.label.fs'])
        axes.set_ylabel('Y',color=Params['axis.label.col'],fontsize=Params['axis.label.fs'])

    return (vmin,vmax)

def pAper(axes, aper, label='', ccdAper=None):
    """
    Plots an :class:`Aperture` object, returning references to the
    plot objects.

    Arguments::

      axes    : (:class:`matplotlib.axes.Axes`)
           the Axes to plot to.

      aper    : (Aperture)
           the :class:`Aperture` to plot

      label   : (string)
           a string label to add

      ccdAper : (CcdAper)
           needed if plotting multiple apertures with links

    Returns a tuple containing references to all the objects created
    when plotting the Aperture in case of a later need to delete them
    """

    # draw circles to represent the aperture. 'objs' is a list of the
    # objects that we keep to return for possible later deletion.
    objs = []

    # target
    objs.append(
        axes.add_patch(Circle(
            (aper.x,aper.y),aper.rtarg,fill=False,
            color=Params['aper.reference.col'] if aper.ref else
            Params['aper.target.col']
    )))

    # sky
    objs.append(
        axes.add_patch(Circle(
            (aper.x,aper.y),aper.rsky1,fill=False,
            color=Params['aper.reference.col'] if aper.ref else
            Params['aper.sky.col']
    )))

    objs.append(
        axes.add_patch(Circle(
            (aper.x,aper.y),aper.rsky2,fill=False,
            color=Params['aper.reference.col'] if aper.ref else
            Params['aper.sky.col']
    )))

    if aper.link != '':
        # indicate a link with an arrow
        if ccdAper is None:
            raise ValueError(
                'to plot a linked aperture, need to pass through an CcdAper'
            )
        else:
            laper = ccdAper[aper.link]

            # draw arrow starting at target aperture of one
            # aperture to the other.
            p1 = utils.Vec2D(aper.x, aper.y)
            p2 = utils.Vec2D(laper.x, laper.y)
            v = p2-p1
            uv = v.unit()
            r1 = aper.rtarg*uv
            r2 = laper.rtarg*uv
            v -= (aper.rtarg+laper.rtarg)*uv
            p1 += r1
            objs.append(
                axes.arrow(
                    p1.x, p1.y, v.x, v.y, width=0.5, length_includes_head=True,
                    overhang=0.8, lw=2, color=Params['aper.link.col']
                )
            )

    # draw dashed lines connecting the aperture to the centres of mask
    # indicated with circles. NB plot returns a list of the lines added, only
    # one in this case, so we extract it rather than storing a list
    for xoff,yoff,r in aper.mask:
        # draw the line
        objs.append(
            axes.plot([aper.x, aper.x+xoff], [aper.y, aper.y+yoff], '--',
                      color=Params['aper.mask.col'])[0]
            )

        # and now the circle
        objs.append(
            axes.add_patch(Circle((aper.x+xoff,aper.y+yoff),r,fill=False,ls='dashed',
                                  color=Params['aper.mask.col']))
            )

    # draw dashed lines connecting the aperture to the centres of mask
    # indicated with circles. NB plot returns a list of the lines added, only
    # one in this case, so we extract it rather than storing a list
    for xoff,yoff in aper.extra:
        # draw the line
        objs.append(
            axes.plot(
                [aper.x, aper.x+xoff], [aper.y, aper.y+yoff], '--',
                color=Params['aper.extra.col'])[0]
            )

        # and now the circle
        objs.append(
            axes.add_patch(
                Circle((aper.x+xoff,aper.y+yoff),aper.rtarg,fill=False,ls='dashed',
                       color=Params['aper.extra.col']))
            )

    if label != '':
        objs.append(
            axes.text(
                aper.x-aper.rsky2,aper.y-aper.rsky2,label,
                color=Params['aper.label.col'],ha='right',va='top',
                bbox=dict(ec='0.8',fc='0.8',alpha=0.5)
            )
        )

    return tuple(objs)

def pCcdAper(axes, ccdAper):
    """
    Plots a :class:`CcdAper` object, returning references to the plot
    objects.

    Arguments::

      axes    : (:class:`matplotlib.axes.Axes`)
           the Axes to plot to.

      ccdAper : (CcdAper)
           the :class:`CcdAper` to plot

    Returns a Group keyed on the same keys as ccdAper but containing
    tuples of the plot objects used to plot each Aperture. This can be
    used to delete them if need be.
    """
    g = Group(tuple)
    for key, aper in ccdAper.items():
        objs = pAper(axes, aper, key, ccdAper)
        g[key] = objs

    return g

def pDefect(axes, dfct):
    """Plots a :class:`Defect` object, returning references to the
    plot objects.

    Arguments::

      axes    : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

      dfct    : Defect
           the :class:`Defect` to plot

    Returns a reference to the plotting object created in case when plotting
    the Defect in case of a later need to delete it. For forwards
    compatibility, this is returned in a tuple

    """

    if isinstance(dfct, defect.Point):
        if dfct.severity == defect.Severity.MODERATE:
            obj = axes.plot(
                dfct.x, dfct.y, 'o', color=Params['defect.moderate.col'],
                ms=2.5
            )
        elif dfct.severity == defect.Severity.SEVERE:
            obj = axes.plot(
                dfct.x, dfct.y, 'o', color=Params['defect.severe.col'],
                ms=2.5
            )

    elif isinstance(dfct, defect.Line):
        if dfct.severity == defect.Severity.MODERATE:
            obj = axes.plot(
                [dfct.x1,dfct.x2], [dfct.y1,dfct.y2],
                '--', color=Params['defect.moderate.col']
            )
        elif dfct.severity == defect.Severity.SEVERE:
            obj = axes.plot(
                [dfct.x1,dfct.x2], [dfct.y1,dfct.y2],
                color=Params['defect.severe.col']
            )

    elif isinstance(dfct, defect.Hot):
        if dfct.severity == defect.Severity.MODERATE:
            obj = axes.plot(
                dfct.x, dfct.y, '*', color=Params['defect.moderate.col'],
                ms=2.5
            )
        elif dfct.severity == defect.Severity.SEVERE:
            obj = axes.plot(
                dfct.x, dfct.y, '*', color=Params['defect.severe.col'],
                ms=2.5
            )

    else:
        raise HipercamError('Did not recognise Defect')

    return tuple(obj)

def pCcdDefect(axes, ccdDefect):
    """Plots a :class:`CcdDefect` object, returning references to the plot
    objects.

    Arguments::

      axes      : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

      ccdDefect : CcdDefect
           the :class:`CcdDefect` to plot

    Returns a Group keyed on the same keys as ccdDefect but containing tuples
    of the plot objects used to plot each Defect. This can be used to
    delete them if need be.

    """
    g = Group(tuple)
    for key, dfct in ccdDefect.items():
        g[key] = pDefect(axes, dfct)

    return g
