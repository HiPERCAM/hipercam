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
SCALE_FACTOR = float(os.environ.get("HIPERCAM_MPL_FSCALE",1.))

# some look-and-feel globals.
Params = {
    # axis label fontsize
    "axis.label.fs": 14 * SCALE_FACTOR,
    # axis label colour
    "axis.label.col": CIS[2],
    # axis number fontsize
    "axis.number.fs": 14 * SCALE_FACTOR,
    # axis colour
    "axis.col": CIS[4],
    # window box colour
    "ccd.box.col": CIS[1],
    # window label font size
    "win.label.fs": 12 * SCALE_FACTOR,
    # window label colour
    "win.label.col": CIS[5],
    # window box colour
    "win.box.col": CIS[7],
    # aperture target colour
    "aper.target.col": CIS[1],
    # aperture reference target colour
    "aper.reference.col": CIS[3],
    # aperture sky colour
    "aper.sky.col": CIS[2],
    # aperture label colour
    "aper.label.col": CIS[4],
    # aperture link colour
    "aper.link.col": CIS[5],
    # aperture mask colour
    "aper.mask.col": CIS[5],
    # aperture extra colour
    "aper.extra.col": CIS[1],
    # moderate defect colour
    "defect.moderate.col": CIS[15],
    # severe defect colour
    "defect.severe.col": CIS[2],
}


def pWin(axes, win, label="", animated=False, artists=None):
    """
    Plots boundary of a :class:`Winhead` as a line and optionally
    labels it.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       win : :class:`Winhead`
           the :class:`Winhead` to plot

       label : str
           label to attach to the Winhead.

       animated : bool
           whether any artists created or updated are animated (had
           better be in the update case)

       artists : None | list(Artist)
           If not None, should be a saved version of a list of
           returned from a previous call to this routine with
           animated=True. In this case they are updated rather
           than re-created.

    Returns with a list of created / updated artists.

    """

    left, right, bottom, top = win.extent()

    if artists is None:
        artists = []
        artists.append(
            axes.plot(
                [left, right, right, left, left],
                [bottom, bottom, top, top, bottom],
                color=Params["win.box.col"],
                animated=animated
            )[0]
        )

        if label != "":
            artists.append(
                axes.text(
                    left - 3, bottom - 3,
                    label,
                    fontsize=Params["win.label.fs"],
                    color=Params["win.label.col"],
                    ha="right",
                    va="top",
                    clip_on=True,
                    animated=animated
                )
            )
    else:
        # just update. not sure why one would ever want to do this ...
        line = artists[0]
        line.set_data(
            [left, right, right, left, left], [bottom, bottom, top, top, bottom]
        )

        if label != "":
            text = artists[1]
            text.set_position(
                (left - 3, bottom - 3)
            )

    return artists

def pWind(axes, wind, vmin, vmax, label="", cmap="Greys", animated=False, artists=None):
    """Plots :class:`Window` as an image with a line border.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       wind : Window
           the :class:`Window` to plot

       vmin : float
           image value at minimum intensity

       vmax : float
           image value at maximum intensity

       label : str
           label to attach to the Window

       animated : bool
           whether any artists created or updated are animated (had
           better be in the update case)

       artists : None | list(Artist)
           If not None, should be a saved version of a list of
           returned from a previous call to this routine with
           animated=True. In this case they are updated rather
           than re-created. NB the border is assumed fixed and
           not updated, i.e. only the image is updated for new
           data and scaling.

    Returns with a list of created / updated artists.
    """
    left, right, bottom, top = wind.extent()

    if artists is None:

        artists = []

        # plot image
        artists.append(
            axes.imshow(
                wind.data,
                extent=(left, right, bottom, top),
                aspect="equal",
                origin="lower",
                cmap=cmap,
                interpolation="nearest",
                vmin=vmin,
                vmax=vmax,
                animated=animated
            )
        )

        # plot window border
        alist = pWin(axes, wind, label, animated)
        artists += alist

    else:
        # update only

        img = artists[0]
        img.set_data(wind.data)
        img.set_clim(vmin,vmax)

    return artists

def pCcd(
        axes, ccd,
        iset="p", plo=5.0, phi=95.0, dlo=0.0, dhi=1000.0,
        tlabel="", xlabel="X", ylabel="Y",
        xlo=None, xhi=None, ylo=None, yhi=None,
        cmap="Greys",
        animated=False, artists=None
):
    """Plots :class:`CCD` as a set of :class:`Window` objects correctly
    positioned with respect to each other.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       ccd : CCD
           the :class:`CCD` to plot

       iset : str
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

       tlabel : str
           label to use for top of plot

       xlabel : str
           label for X axis

       ylabel : str
           label for Y axis

       xlo : int | None
           left-hand limit to define region for computing percentile

       xhi : int | None
           right-hand limit to define region for computing percentile

       ylo : int | None
           lower limit to define region for computing percentile

       yhi : int | None
           upper limit to define region for computing percentile

       animated : bool
           whether any artists created or updated are animated (had
           better be in the update case)

       artists : None | dict
           If not None, should be a saved version of the dict of
           returned from a previous call to this routine with
           animated=True. In this case they are updated rather
           than re-created.

    Returns: (vmin,vmax,artists), the intensity limits used and
    a dictionary, the values of which are all lists of artists.

    """
    if iset == "p":
        # Set intensities from percentiles
        if xlo is not None and xhi is not None:
            xlo = min(xlo, xhi)
            xhi = max(xlo, xhi)
        if ylo is not None and yhi is not None:
            ylo = min(ylo, yhi)
            yhi = max(ylo, yhi)
        vmin, vmax = ccd.percentile(
            (plo, phi), xlo, xhi, ylo, yhi
        )

    elif iset == "a":
        # Set intensities from min/max range
        vmin, vmax = ccd.min(), ccd.max()

    elif iset == "d":
        vmin, vmax = dlo, dhi

    else:
        raise ValueError('did not recognise iset = "' + iset + '"')

    if artists is None:
        # create the artists
        artists = {}
        for key, wind in ccd.items():
            artists[key] = pWind(axes, wind, vmin, vmax, key, cmap, animated)

        # plot outermost border of CCD
        artists['border'] = \
            axes.plot(
                [0.5, ccd.nxtot + 0.5, ccd.nxtot + 0.5, 0.5, 0.5],
                [0.5, 0.5, ccd.nytot + 0.5, ccd.nytot + 0.5, 0.5],
                color=Params["ccd.box.col"], animated=animated
            )

        artists['labels'] = []
        if tlabel != "":
            artists['labels'].append(
                axes.set_title(
                    tlabel, color=Params["axis.label.col"],
                    fontsize=Params["axis.label.fs"], animated=animated
                )
            )

        if xlabel != "":
            artists['labels'].append(
                axes.set_xlabel(
                    xlabel, color=Params["axis.label.col"],
                    fontsize=Params["axis.label.fs"], animated=animated
                )
            )

        if ylabel != "":
            artists['labels'].append(
                axes.set_ylabel(
                    ylabel, color=Params["axis.label.col"],
                    fontsize=Params["axis.label.fs"], animated=animated
                )
            )

    else:
        # only need to update the images
        for key, wind in ccd.items():
            artists[key] = pWind(axes, wind, vmin, vmax, key, cmap, animated, artists[key])

    return (vmin, vmax, artists)

def pAper(axes, aper, label="", ccdAper=None, animated=False, artists=None):
    """Plots an :class:`Aperture` object, returning a list of the artists
    created or updated.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       aper : Aperture
           the :class:`Aperture` to plot

       label : str
           a string label to add

       ccdAper : CcdAper
           needed if plotting multiple apertures with links

       animated : bool
           whether to animate the artist (when creating). This allows later updating
           on the fly using the next argument.

       artists : None | list(artist)
           If not None, this should be the artists returned from a
           previous call to this routine with animated=True. Then the
           objects will be updated rather than created from scratch.
           NB This assumes that the basic structure, i.e. masks, extra
           apertures and links is unchanged, and only updates sizes and
           positions

    Returns a list containing references to all the artists created
    when plotting the Aperture in case of a later need to delete them,
    or to update them by passing back into this routine. Since an Aperture
    can feature masks, extras and links, the list can be quite a long one.

    """

    if artists is None:
        # Draw circles & lines to represent the aperture. Retain list of the
        # artists created
        artists = []

        # target
        artists.append(
            axes.add_patch(
                Circle(
                    (aper.x, aper.y), aper.rtarg, fill=False,
                    color=Params["aper.reference.col"]
                    if aper.ref else Params["aper.target.col"],
                    animated=animated
                )
            )
        )

        # sky
        artists.append(
            axes.add_patch(
                Circle(
                    (aper.x, aper.y), aper.rsky1, fill=False,
                    color=Params["aper.reference.col"]
                    if aper.ref else Params["aper.sky.col"],
                    animated=animated
                )
            )
        )

        artists.append(
            axes.add_patch(
                Circle(
                    (aper.x, aper.y), aper.rsky2, fill=False,
                    color=Params["aper.reference.col"]
                    if aper.ref else Params["aper.sky.col"],
                    animated=animated
                )
            )
        )

        if aper.link != "":
            # indicate a link with an arrow
            if ccdAper is None:
                raise ValueError(
                    "to plot a linked aperture, need to pass through an CcdAper"
                )
            else:
                laper = ccdAper[aper.link]

                # draw arrow starting at target aperture of one
                # aperture to the other.
                p1 = utils.Vec2D(aper.x, aper.y)
                p2 = utils.Vec2D(laper.x, laper.y)
                v = p2 - p1
                uv = v.unit()
                r1 = aper.rtarg * uv
                r2 = laper.rtarg * uv
                v -= (aper.rtarg + laper.rtarg) * uv
                p1 += r1
                artists.append(
                    axes.arrow(
                        p1.x, p1.y, v.x, v.y,
                        width=3.0, length_includes_head=True,
                        overhang=0.5, head_width=18,
                        head_length=30, lw=2,
                        color=Params["aper.link.col"],
                        animated=animated
                    )
                )

        # draw dashed lines connecting the aperture to the centres of mask
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff, r in aper.mask:
            # draw the line
            artists.append(
                axes.plot(
                    [aper.x, aper.x + xoff],
                    [aper.y, aper.y + yoff],
                    "--",
                    color=Params["aper.mask.col"],
                    animated=animated
                )[0]
            )

            # and now the circle
            artists.append(
                axes.add_patch(
                    Circle(
                        (aper.x + xoff, aper.y + yoff), r,
                        fill=False, ls="dashed",
                        color=Params["aper.mask.col"],
                        animated=animated
                    )
                )
            )

        # draw dashed lines connecting the aperture to the centres of extras
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff in aper.extra:
            # draw the line
            artists.append(
                axes.plot(
                    [aper.x, aper.x + xoff],
                    [aper.y, aper.y + yoff],
                    "--",
                    color=Params["aper.extra.col"],
                    animated=animated
                )[0]
            )

            # and now the circle
            artists.append(
                axes.add_patch(
                    Circle(
                        (aper.x + xoff, aper.y + yoff), aper.rtarg,
                        fill=False, ls="dashed",
                        color=Params["aper.extra.col"],
                        animated=animated
                    )
                )
            )

        if label != "":
            artists.append(
                axes.text(
                    aper.x - aper.rsky2,
                    aper.y - aper.rsky2,
                    label,
                    color=Params["aper.label.col"],
                    ha="right",
                    va="top",
                    bbox=dict(ec="0.8", fc="0.8", alpha=0.5),
                    animated=True
                )
            )
    else:

        # Update rather than create the artists. The update assumes
        # only that positions and sizes change. One exception is the
        # arrow used for links which is first removed then re-created.

        # counter within the list of artists
        n = 0

        # target
        circ = artists[n]
        circ.set_center((aper.x, aper.y))
        circ.set_radius(aper.rtarg)
        n += 1

        # sky
        circ = artists[n]
        circ.set_center((aper.x, aper.y))
        circ.set_radius(aper.rsky1)
        n += 1

        circ = artists[n]
        circ.set_center((aper.x, aper.y))
        circ.set_radius(aper.rsky2)
        n += 1

        if aper.link != "":
            # indicate a link with an arrow
            if ccdAper is None:
                raise ValueError(
                    "to plot a linked aperture, need to pass through a CcdAper"
                )
            else:
                laper = ccdAper[aper.link]

                # draw arrow starting at target aperture of one
                # aperture to the other.
                p1 = utils.Vec2D(aper.x, aper.y)
                p2 = utils.Vec2D(laper.x, laper.y)
                v = p2 - p1
                uv = v.unit()
                r1 = aper.rtarg * uv
                r2 = laper.rtarg * uv
                v -= (aper.rtarg + laper.rtarg) * uv
                p1 += r1

                # Seems to be no easy way to update a FancyArrow
                # just remove then recreate. Hope it works.
                artists[n].remove()
                artists[n] = axes.arrow(
                    p1.x, p1.y, v.x, v.y,
                    width=3.0, length_includes_head=True,
                    overhang=0.5, head_width=18,
                    head_length=30, lw=2,
                    color=Params["aper.link.col"],
                    animated=animated
                )
                n += 1

        # draw dashed lines connecting the aperture to the centres of mask
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff, r in aper.mask:
            # draw the line
            line = artists[n]
            line.set_data(
                [aper.x, aper.x + xoff],
                [aper.y, aper.y + yoff]
            )
            n += 1

            # and now the circle
            circ = artists[n]
            circ.set_center(
                (aper.x + xoff, aper.y + yoff)
            )
            circ.set_radius(r)
            n += 1

        # draw dashed lines connecting the aperture to the centres of extras
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff in aper.extra:
            # draw the line
            line = artists[n]
            line.set_data(
                [aper.x, aper.x + xoff], [aper.y, aper.y + yoff]
            )
            n += 1

            # and now the circle
            circ = artists[n]
            circ.set_center(
                (aper.x + xoff, aper.y + yoff)
            )
            circ.set_radius(aper.rtarg)
            n += 1

        if label != "":
            text = artists[n]
            text.set_position(
                (aper.x - aper.rsky2, aper.y - aper.rsky2)
            )
            n += 1

    return artists

def pCcdAper(axes, ccdAper, animated=False, artists=None):
    """Plots a :class:`CcdAper` object, returning references to the plot
    objects.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       ccdAper : CcdAper
           the :class:`CcdAper` to plot

       animated : bool
           whether any artists created or updated are animated (had
           better be in the update case)

       artists : None | dict
           If not None, should be a saved version of the dict of
           returned from a previous call to this routine with
           animated=True. In this case they are updated rather
           than re-created.

    Returns a dict keyed on the same keys as ccdAper but containing
    the artists used to plot each Aperture. This can be used to delete
    or update them if need be.
    """

    pobjs = {}
    if ccdAper is not None:
        for key, aper in ccdAper.items():
            if artists is None:
                pobjs[key] = pAper(axes, aper, key, ccdAper, animated)
            else:
                pobjs[key] = pAper(axes, aper, key, ccdAper, animated, artists[key])

    return pobjs

def pDefect(axes, dfct, animated=False, artists=None):
    """Plots a :class:`Defect` object, returning a list of
    the artists created.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       dfct : Defect
           the :class:`Defect` to plot

       animated : bool
           whether to animate the artist (when creating). This allows later updating
           on the fly using the next argument.

       artists : None | list(artist)
           If not None, this should be the artists returned from a
           previous call to this routine with animated=True. Then the
           objects will be updated rather than created from scratch.

    Returns a reference to the list of artists created in case of a
    later need to delete or update them.

    """

    if artists is None:
        # create plot objects

        if isinstance(dfct, defect.Point):
            if dfct.severity == defect.Severity.MODERATE:
                artists = axes.plot(
                    dfct.x, dfct.y, "o", color=Params["defect.moderate.col"], ms=2.5,
                    animated=animated
                )
            elif dfct.severity == defect.Severity.SEVERE:
                artists = axes.plot(
                    dfct.x, dfct.y, "o", color=Params["defect.severe.col"], ms=2.5,
                    animated=animated
                )

        elif isinstance(dfct, defect.Line):
            if dfct.severity == defect.Severity.MODERATE:
                artists = axes.plot(
                    [dfct.x1, dfct.x2],
                    [dfct.y1, dfct.y2],
                    "--",
                    color=Params["defect.moderate.col"],
                    animated=animated
                )
            elif dfct.severity == defect.Severity.SEVERE:
                artists = axes.plot(
                    [dfct.x1, dfct.x2],
                    [dfct.y1, dfct.y2],
                    color=Params["defect.severe.col"],
                    animated=animated
                )

        elif isinstance(dfct, defect.Hot):
            if dfct.severity == defect.Severity.MODERATE:
                artists = axes.plot(
                    dfct.x, dfct.y, "*", color=Params["defect.moderate.col"], ms=2.5,
                    animated=animated
                )
            elif dfct.severity == defect.Severity.SEVERE:
                artists = axes.plot(
                    dfct.x, dfct.y, "*", color=Params["defect.severe.col"], ms=2.5,
                    animated=animated
                )

        else:
            raise HipercamError("Did not recognise Defect")

    else:
        # update previously created artists
        if isinstance(dfct, defect.Point) or isinstance(dfct, defect.Hot):
            artists[0].set_data(dfct.x, dfct.y)

        elif isinstance(dfct, defect.Line):
            artists[0].set_data([dfct.x1, dfct.x2], [dfct.y1, dfct.y2])

        else:
            raise HipercamError("Did not recognise Defect")

    return artists

def pCcdDefect(axes, ccdDefect, animated=False, artists=None):
    """Plots a :class:`CcdDefect` object, returning references to the
    artists.

    Arguments::

       axes      : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       ccdDefect : CcdDefect
           the :class:`CcdDefect` to plot

       animated : bool
           whether any artists created or updated are animated (had
           better be in the update case)

       artists : None | dict
           If not None, should be a saved version of the dict of
           returned from a previous call to this routine with
           animated=True. In this case they are updated rather
           than re-created.

    Returns a dict keyed on the same keys as ccdDefect with values
    storing the artists for each Defect. This can be used to delete
    or update them (if animated) if need be.

    """
    pobjs = {}
    if ccdDefect is not None:
        for key, dfct in ccdDefect.items():
            if artists is None:
                pobjs[key] = pDefect(axes, dfct, animated)
            else:
                pobjs[key] = pDefect(axes, dfct, animated, artists[key])

    return pobjs

def pFringePair(axes, fpair):
    """Plots a :class:`FringePair` object, returning
    a list of artists created.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       fpair : FringePair
           the :class:`FringePair` to plot

    """

    return axes.plot(
        [fpair.x1,fpair.x2], [fpair.y1,fpair.y2],
        mfc='b', mec='b', color='r', marker='o', lw=4
    )

def pCcdFringePair(axes, ccdFringePair):
    """Plots a :class:`CcdFringePair` object, returning references to the plot
    objects.

    Arguments::

       axes : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       ccdFringePair : CcdFringePair
           the :class:`CcdFringePair` to plot

    Returns a Group keyed on the same keys as ccdFringePair but
    containing tuples of the plot objects used to plot each
    Defect. This can be used to delete them if need be.

    """
    pobjs = {}
    for key, fpair in ccdFringePair.items():
        pobjs[key] = pFringePair(axes, fpair)

    return pobjs
