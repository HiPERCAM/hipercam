# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
matplotlib plotting functions.
"""

from __future__ import annotations

import os
import re
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from . import defect, utils
from .cline import Cline
from .core import CIS, HipercamError

if TYPE_CHECKING:
    from .ccd import CCD, MCCD
    from .defect import CcdDefect


# Font scale factor
SCALE_FACTOR = float(os.environ.get("HIPERCAM_MPL_FSCALE", "1.0"))

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
    "defect.moderate.col": CIS[8],
    # severe defect colour
    "defect.severe.col": CIS[2],
    # setup window box colour
    "setup.box.col": (0.8, 0.0, 0.0),
}


def pWin(ax, win, label="", animated=False, artists=None):
    """
    Plots boundary of a :class:`Winhead` as a line and optionally
    labels it.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
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
            ax.plot(
                [left, right, right, left, left],
                [bottom, bottom, top, top, bottom],
                color=Params["win.box.col"],
                animated=animated,
            )[0]
        )

        if label != "":
            artists.append(
                ax.text(
                    left - 3,
                    bottom - 3,
                    label,
                    fontsize=Params["win.label.fs"],
                    color=Params["win.label.col"],
                    ha="right",
                    va="top",
                    clip_on=True,
                    animated=animated,
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
            text.set_position((left - 3, bottom - 3))

    return artists


def pWind(ax, wind, vmin, vmax, label="", cmap="Greys", animated=False, artists=None):
    """Plots :class:`Window` as an image with a line border.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
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
            ax.imshow(
                wind.data,
                extent=(left, right, bottom, top),
                aspect="equal",
                origin="lower",
                cmap=cmap,
                interpolation="nearest",
                vmin=vmin,
                vmax=vmax,
                animated=animated,
            )
        )

        # plot window border
        alist = pWin(ax, wind, label, animated)
        artists += alist

    else:
        # update only

        img = artists[0]
        img.set_data(wind.data)
        img.set_clim(vmin, vmax)

    return artists


def pCcd(
    ax,
    ccd,
    vmin=None,
    vmax=None,
    cmap="Greys",
    tlabel="",
    xlabel="X",
    ylabel="Y",
    animated=False,
    artists=None,
):
    """Plots :class:`CCD` as a set of :class:`Window` objects correctly
    positioned with respect to each other.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       ccd : CCD
           the :class:`CCD` to plot

       vmin : float | None
           image value at minimum intensity

       vmax : float | None
           image value at maximum intensity

       cmap : str
           the name of the matplotlib colormap to use

       tlabel : str
           label to use for top of plot

       xlabel : str
           label for X axis

       ylabel : str
           label for Y axis

       animated : bool
           whether any artists created or updated are animated (had
           better be in the update case)

       artists : None | dict
           If not None, should be a saved version of the dict of
           returned from a previous call to this routine with
           animated=True. In this case they are updated rather
           than re-created.

    Returns: a dictionary, the values of which are all artists or lists of artists:
    - <wnam>:  the artists for Windows with name <wnam>, returned from pWind
    - "border" : the artist for the CCD border
    - "labels" : list of artists for the title and axis labels

    """
    if artists is None:
        # create the artists
        artists = {}

        # plot each window in the CCD
        for key, wind in ccd.items():
            artists[key] = pWind(ax, wind, vmin, vmax, key, cmap, animated)

        # plot outermost border of CCD
        artists["border"] = ax.plot(
            [0.5, ccd.nxtot + 0.5, ccd.nxtot + 0.5, 0.5, 0.5],
            [0.5, 0.5, ccd.nytot + 0.5, ccd.nytot + 0.5, 0.5],
            color=Params["ccd.box.col"],
            animated=animated,
        )

        ax.tick_params(labelsize=Params["axis.number.fs"])

        # plot labels
        artists["labels"] = []
        if tlabel != "":
            artists["labels"].append(
                ax.set_title(
                    tlabel,
                    color=Params["axis.label.col"],
                    fontsize=Params["axis.label.fs"],
                    animated=animated,
                )
            )
        if xlabel != "":
            artists["labels"].append(
                ax.set_xlabel(
                    xlabel,
                    color=Params["axis.label.col"],
                    fontsize=Params["axis.label.fs"],
                    animated=animated,
                )
            )
        if ylabel != "":
            artists["labels"].append(
                ax.set_ylabel(
                    ylabel,
                    color=Params["axis.label.col"],
                    fontsize=Params["axis.label.fs"],
                    animated=animated,
                )
            )

    else:
        # only need to update the images
        for key, wind in ccd.items():
            artists[key] = pWind(
                ax, wind, vmin, vmax, key, cmap, animated, artists[key]
            )

    return artists


def pAper(ax, aper, label="", ccdAper=None, animated=False, artists=None):
    """Plots an :class:`Aperture` object, returning a list of the artists
    created or updated.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
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
            ax.add_patch(
                Circle(
                    (aper.x, aper.y),
                    aper.rtarg,
                    fill=False,
                    color=(
                        Params["aper.reference.col"]
                        if aper.ref
                        else Params["aper.target.col"]
                    ),
                    ls="dashdot" if aper.compo else "solid",
                    animated=animated,
                )
            )
        )

        # sky
        artists.append(
            ax.add_patch(
                Circle(
                    (aper.x, aper.y),
                    aper.rsky1,
                    fill=False,
                    color=(
                        Params["aper.reference.col"]
                        if aper.ref
                        else Params["aper.sky.col"]
                    ),
                    ls="dashdot" if aper.compo else "solid",
                    animated=animated,
                )
            )
        )

        artists.append(
            ax.add_patch(
                Circle(
                    (aper.x, aper.y),
                    aper.rsky2,
                    fill=False,
                    color=(
                        Params["aper.reference.col"]
                        if aper.ref
                        else Params["aper.sky.col"]
                    ),
                    ls="dashdot" if aper.compo else "solid",
                    animated=animated,
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
                    ax.arrow(
                        p1.x,
                        p1.y,
                        v.x,
                        v.y,
                        width=3.0,
                        length_includes_head=True,
                        overhang=0.5,
                        head_width=18,
                        head_length=30,
                        lw=2,
                        color=Params["aper.link.col"],
                        animated=animated,
                    )
                )

        # draw dashed lines connecting the aperture to the centres of mask
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff, r in aper.mask:
            # draw the line
            artists.append(
                ax.plot(
                    [aper.x, aper.x + xoff],
                    [aper.y, aper.y + yoff],
                    "--",
                    color=Params["aper.mask.col"],
                    animated=animated,
                )[0]
            )

            # and now the circle
            artists.append(
                ax.add_patch(
                    Circle(
                        (aper.x + xoff, aper.y + yoff),
                        r,
                        fill=False,
                        ls="dashed",
                        color=Params["aper.mask.col"],
                        animated=animated,
                    )
                )
            )

        # draw dashed lines connecting the aperture to the centres of extras
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff in aper.extra:
            # draw the line
            artists.append(
                ax.plot(
                    [aper.x, aper.x + xoff],
                    [aper.y, aper.y + yoff],
                    "--",
                    color=Params["aper.extra.col"],
                    animated=animated,
                )[0]
            )

            # and now the circle
            artists.append(
                ax.add_patch(
                    Circle(
                        (aper.x + xoff, aper.y + yoff),
                        aper.rtarg,
                        fill=False,
                        ls="dashed",
                        color=Params["aper.extra.col"],
                        animated=animated,
                    )
                )
            )

        if label != "":
            artists.append(
                ax.text(
                    aper.x - aper.rsky2,
                    aper.y - aper.rsky2,
                    label,
                    color=Params["aper.label.col"],
                    ha="right",
                    va="top",
                    bbox=dict(ec="0.8", fc="0.8", alpha=0.7, boxstyle="round, pad=0"),
                    clip_on=True,
                    animated=animated,
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
                artists[n] = ax.arrow(
                    p1.x,
                    p1.y,
                    v.x,
                    v.y,
                    width=3.0,
                    length_includes_head=True,
                    overhang=0.5,
                    head_width=18,
                    head_length=30,
                    lw=2,
                    color=Params["aper.link.col"],
                    animated=animated,
                )
                n += 1

        # draw dashed lines connecting the aperture to the centres of mask
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff, r in aper.mask:
            # draw the line
            line = artists[n]
            line.set_data([aper.x, aper.x + xoff], [aper.y, aper.y + yoff])
            n += 1

            # and now the circle
            circ = artists[n]
            circ.set_center((aper.x + xoff, aper.y + yoff))
            circ.set_radius(r)
            n += 1

        # draw dashed lines connecting the aperture to the centres of extras
        # indicated with circles. NB plot returns a list of the lines added, only
        # one in this case, so we extract it rather than storing a list
        for xoff, yoff in aper.extra:
            # draw the line
            line = artists[n]
            line.set_data([aper.x, aper.x + xoff], [aper.y, aper.y + yoff])
            n += 1

            # and now the circle
            circ = artists[n]
            circ.set_center((aper.x + xoff, aper.y + yoff))
            circ.set_radius(aper.rtarg)
            n += 1

        if label != "":
            text = artists[n]
            text.set_position((aper.x - aper.rsky2, aper.y - aper.rsky2))
            n += 1

    return artists


def pCcdAper(ax, ccdAper, animated=False, artists=None):
    """Plots a :class:`CcdAper` object, returning references to the plot
    objects.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
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
                pobjs[key] = pAper(ax, aper, key, ccdAper, animated)
            else:
                pobjs[key] = pAper(ax, aper, key, ccdAper, animated, artists[key])

    return pobjs


def pDefect(ax, dfct, animated=False, artists=None):
    """Plots a :class:`Defect` object, returning a list of
    the artists created.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
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
                artists = ax.plot(
                    dfct.x,
                    dfct.y,
                    "o",
                    color=Params["defect.moderate.col"],
                    ms=3,
                    animated=animated,
                )
            elif dfct.severity == defect.Severity.SEVERE:
                artists = ax.plot(
                    dfct.x,
                    dfct.y,
                    "o",
                    color=Params["defect.severe.col"],
                    ms=3,
                    animated=animated,
                )

        elif isinstance(dfct, defect.Line):
            if dfct.severity == defect.Severity.MODERATE:
                artists = ax.plot(
                    [dfct.x1, dfct.x2],
                    [dfct.y1, dfct.y2],
                    "--",
                    color=Params["defect.moderate.col"],
                    animated=animated,
                )
            elif dfct.severity == defect.Severity.SEVERE:
                artists = ax.plot(
                    [dfct.x1, dfct.x2],
                    [dfct.y1, dfct.y2],
                    color=Params["defect.severe.col"],
                    animated=animated,
                )

        elif isinstance(dfct, defect.Hot):
            if dfct.severity == defect.Severity.MODERATE:
                col = Params["defect.moderate.col"]
            elif dfct.severity == defect.Severity.SEVERE:
                col = Params["defect.severe.col"]

            artists = [
                ax.text(
                    dfct.x,
                    dfct.y,
                    f"{int(round(dfct.rate))}",
                    color=col,
                    ha="center",
                    va="center",
                    bbox=dict(ec="0.8", fc="0.8", alpha=0.5, boxstyle="round, pad=0"),
                    clip_on=True,
                    animated=animated,
                )
            ]

        else:
            raise HipercamError(f"Did not recognise Defect: {type(dfct)}")

    else:
        # update previously created artists
        if isinstance(dfct, defect.Point):
            artists[0].set_data((dfct.x,), (dfct.y,))

        elif isinstance(dfct, defect.Line):
            artists[0].set_data([dfct.x1, dfct.x2], [dfct.y1, dfct.y2])

        elif isinstance(dfct, defect.Hot):
            artists[0].set_position((dfct.x, dfct.y))

        else:
            raise HipercamError("Did not recognise Defect")

    return artists


def pCcdDefect(ax, ccdDefect, animated=False, artists=None):
    """Plots a :class:`CcdDefect` object, returning references to the
    artists.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
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
                pobjs[key] = pDefect(ax, dfct, animated)
            else:
                pobjs[key] = pDefect(ax, dfct, animated, artists[key])

    return pobjs


def pFringePair(ax, fpair):
    """Plots a :class:`FringePair` object, returning
    a list of artists created.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       fpair : FringePair
           the :class:`FringePair` to plot

    """

    return ax.plot(
        [fpair.x1, fpair.x2],
        [fpair.y1, fpair.y2],
        mfc="b",
        mec="b",
        color="r",
        marker="o",
        lw=4,
    )


def pCcdFringePair(ax, ccdFringePair):
    """Plots a :class:`CcdFringePair` object, returning references to the plot
    objects.

    Arguments::

       ax : :class:`matplotlib.axes.Axes`
           the Axes to plot to.

       ccdFringePair : CcdFringePair
           the :class:`CcdFringePair` to plot

    Returns a Group keyed on the same keys as ccdFringePair but
    containing tuples of the plot objects used to plot each
    Defect. This can be used to delete them if need be.

    """
    pobjs = {}
    for key, fpair in ccdFringePair.items():
        pobjs[key] = pFringePair(ax, fpair)

    return pobjs


def get_scale_limits(ccd, mode, limits, region=None):
    """Returns the intensity limits for a CCD based on the given mode.

    Arguments::
        ccd : :class:`CCD`
            The CCD to compute intensity limits for.

        mode : str
            How to set the intensity scale for the images:
                'a' for automatic min/max range
                'd' for direct values (set using `limits`)
                'p' for percentiles (set using `limits`)

        limits : tuple[float, float] | None
            Lower and upper limits for the intensity scaling:
                If mode='d' these are the lower and upper intensity limits
                If mode='p' these are the lower and upper percentiles
                If mode='a' this argument is ignored
            If None, default values are used:
                (0.0, 1000.0) if mode='d'
                (5.0, 95.0) if mode='p'

        region : tuple[int, int, int, int] | None
            The region to consider for computing intensity limits,
            values are (xmin, xmax, ymin, ymax).
            Only used if mode is 'p'. If None, the entire CCD is considered.

    Returns::
        vmin, vmax : tuple[float, float]
            The computed minimum and maximum intensity limits.
    """
    if mode not in ("a", "p", "d"):
        raise ValueError(f"Invalid mode '{mode}'. Must be 'a', 'p', or 'd'.")
    if limits is None:
        if mode == "p":
            limits = (5.0, 95.0)
        elif mode == "d":
            limits = (0.0, 1000.0)
        else:  # mode == "a"
            limits = (None, None)

    if mode == "a":
        # Set intensities from min/max range
        vmin, vmax = ccd.min(), ccd.max()

    elif mode == "d":
        # Set intensities from direct values
        if limits[0] > limits[1]:
            raise ValueError("Lower limit must be less than or equal to upper limit.")
        vmin, vmax = limits[0], limits[1]

    elif mode == "p":
        # Set intensities from percentiles for the given region
        if limits[0] > limits[1]:
            raise ValueError("Lower limit must be less than or equal to upper limit.")
        if limits[0] < 0.0 or limits[1] > 100.0:
            raise ValueError("Percentile limits must be in the range [0, 100].")
        if region is not None:
            xlo, xhi, ylo, yhi = region
        else:
            xlo, xhi, ylo, yhi = None, None, None, None
        if xlo is not None and xhi is not None:
            pxlo = min(xlo, xhi)
            pxhi = max(xlo, xhi)
        else:
            pxlo, pxhi = None, None
        if ylo is not None and yhi is not None:
            pylo = min(ylo, yhi)
            pyhi = max(ylo, yhi)
        else:
            pylo, pyhi = None, None
        vmin, vmax = ccd.percentile((limits[0], limits[1]), pxlo, pxhi, pylo, pyhi)

    return vmin, vmax


def get_window_position(fig):
    """Extracts the position of an open figure window."""
    manager = fig.canvas.manager
    backend = plt.get_backend().lower()

    if "qt" in backend:
        # Qt uses QRect: x(), y()
        rect = manager.window.geometry()
        xpos, ypos = rect.x(), rect.y()
        return [xpos, ypos]

    elif "tk" in backend:
        # Tk returns a string formatted as 'WxH+X+Y', or '-X-Y' for negative positions
        geom = manager.window.geometry()
        match = re.match(r"(\d+)x(\d+)([+-]\d+)([+-]\d+)", geom)
        _width, _height, xpos, ypos = match.groups()
        return [int(xpos), int(ypos)]

    elif "wx" in backend:
        # wxPython uses tuples for position and size
        xpos, ypos = manager.frame.GetPosition()
        return [xpos, ypos]

    else:
        raise NotImplementedError(
            f"Window position extraction not implemented for backend '{backend}', "
            "set memory=False."
        )


def set_window_position(fig, xpos, ypos):
    """Sets the position of the figure window."""
    manager = fig.canvas.manager
    backend = plt.get_backend().lower()
    xpos = int(xpos)
    ypos = int(ypos)

    if "qt" in backend:
        rect = manager.window.geometry()
        manager.window.setGeometry(xpos, ypos, rect.width(), rect.height())
    elif "tk" in backend:
        xpos_str = f"-{-xpos}" if xpos < 0 else f"+{xpos}"
        ypos_str = f"-{-ypos}" if ypos < 0 else f"+{ypos}"
        manager.window.geometry(f"{xpos_str}{ypos_str}")
    elif "wx" in backend:
        width, height = manager.frame.GetSize()
        manager.frame.SetSize(xpos, ypos, width, height)


def store_window(fig, command, args, prefix="im"):
    """Store the figure size and position in the Cline file for a script.

    This is to replicate the PGPLOT behaviour of keeping the windows open
    after a script has finished. We can't do that without a long-lasting server,
    but this allows you to position the window and have it reopen in the same
    position when the script is restarted.

    Arguments::

        fig : matplotlib.figure.Figure
            The figure to store the size and position for.

        command : str
            The command to store the window size and position for.
            This is sent directly to Cline, so should be the name of the script
            that is calling this function.

        args : list[str]
            The arguments to the script command, which will be stored in the Cline file.
            See `trm.cline.script_args()` for details.

        prefix : str
            The prefix to use for the parameter names in the Cline file.
            Parameters will be stored as:
                {prefix}width, {prefix}height, {prefix}xpos, {prefix}ypos
            Default is "im", so the parameters will be:
                "imwidth", "imheight", "imxpos", "imypos"

    """
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:
        # get config param names
        width_param = f"{prefix}width"
        height_param = f"{prefix}height"
        xpos_param = f"{prefix}xpos"
        ypos_param = f"{prefix}ypos"

        # get and store the current figure size (in inches)
        width = fig.get_figwidth()
        height = fig.get_figheight()
        cl.register(width_param, Cline.LOCAL, Cline.HIDE)
        cl.register(height_param, Cline.LOCAL, Cline.HIDE)
        cl.set_default(width_param, width)
        cl.set_default(height_param, height)

        # get and store the current screen position of the figure (in pixels)
        xpos, ypos = get_window_position(fig)
        cl.register(xpos_param, Cline.LOCAL, Cline.HIDE)
        cl.register(ypos_param, Cline.LOCAL, Cline.HIDE)
        cl.set_default(xpos_param, xpos)
        cl.set_default(ypos_param, ypos)

    return width, height, xpos, ypos


def create_figure(
    cnams: list[str],
    ncols: int,
    region: tuple[int, int, int, int] | None = None,
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
    position: tuple[int, int] | None = None,
):
    """
    Sets up the figure and axes for plotting an MCCD.

    This is separate from plot_frame because for an animation,
    the figure and axes are set up once and then reused.

    Arguments::

        cnams : list[str]
            The names of the CCDs to create axes for.

        ncols : int
            The number of columns in the grid of subplots

        region : tuple[int, int, int, int] | None
            The region to plot, values are (xmin, xmax, ymin, ymax)
            If None, the entire CCD is displayed.

        title : str
            The title of the figure, used for the window title.

        figsize : tuple[float, float] | None
            The width and height of the figure, in inches.
            If None, the default size is used.

        position : tuple[int, int] | None
            The position of the figure window on the screen, in pixels.
            If None, the default position is used.

    Returns::
        fig: matplotlib.figure.Figure
            The figure object.

        axs: dict[str, matplotlib.axes.Axes]
            A dictionary of the Axes objects, one per CCD, labeled by their names.
    """
    if region is None:
        region = [None, None, None, None]
    xlim = region[0:2]
    ylim = region[2:4]

    # Create the figure for the image plot
    fig = plt.figure()
    if title is not None:
        fig.canvas.manager.set_window_title(title)
    if figsize is not None and figsize[0] > 0 and figsize[1] > 0:
        fig.set_size_inches(figsize)
    if position is not None:
        set_window_position(fig, position[0], position[1])

    # Define the panel grid
    ncols = min(len(cnams), ncols) if len(cnams) > 1 else 1
    nrows = len(cnams) // ncols if len(cnams) % ncols == 0 else len(cnams) // ncols + 1

    # Create the axes
    axs = {}
    for n, cnam in enumerate(cnams):
        if n == 0:
            ax = first_ax = fig.add_subplot(nrows, ncols, n + 1)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
        else:
            # subsequent axes share the same x and y as the first
            ax = fig.add_subplot(
                nrows,
                ncols,
                n + 1,
                sharex=first_ax,
                sharey=first_ax,
            )
            ax.set_aspect("equal")

        # store the axes for each CCD in a dictionary for later reference
        axs[cnam] = ax

    return fig, axs


def plot_frame(
    mccd: MCCD,
    cnams: str | list[str] = "0",
    ncols: int = 3,
    defects: dict[str, CcdDefect] | None = None,
    msub: bool = False,
    cmap: str = "Greys",
    iset: str = "p",
    ilim: tuple[float, float] | None = None,
    region: tuple[int, int, int, int] | None = None,
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
    position: tuple[int, int] | None = None,
) -> tuple[plt.Figure, dict[str, plt.Axes]]:
    """Creates a plot of an :class:`Mccd` object.

    Each CCD is added to a subplot in a grid.

    Arguments::

        mccd : :class:`MCCD`
            The MCCD to plot

        cnams : str | list[str]
            The name(s) of the CCDs to plot.
            If "0" then all CCDs in the MCCD are plotted.

        ncols : int
            The number of columns in the grid of subplots

        defects: dict[str, :class:`CcdDefect`] or None
            groups of defects to plot on each CCD
            keys should be cnams

        msub : bool
            Whether to subtract the median from each window before plotting

        cmap : str
            The name of the colormap to use for the images

        iset : str
            How to set the intensity scale for the images:
                'a' for automatic min/max range
                'd' for direct values (set using `ilim`)
                'p' for percentiles (set using `ilim`)
            See :func:`get_scale_limits` for details.

        ilim : tuple[float, float] | None
            Lower and upper limits for the intensity scaling:
                If iset='d' these are the lower and upper intensity limits
                If iset='p' these are the lower and upper percentiles
                If iset='a' this argument is ignored
            See :func:`get_scale_limits` for details.

        region : tuple[int, int, int, int] | None
            The region to plot, values are (xmin, xmax, ymin, ymax)
            If None, the entire CCD is displayed.
            If iset='p', only this region is used to compute the percentiles for
            intensity scaling.

        title : str | None
            The title of the figure, used for the window title.

        figsize : tuple[float, float] | None
            The width and height of the figure, in inches.
            If None, the default size is used.

        position : tuple[int, int] | None
            The position of the figure window on the screen, in pixels.
            If None, the default position is used.

    Returns:
        fig : :class:`matplotlib.figure.Figure`
            The Figure object.

        axs : dict[str, :class:`matplotlib.axes.Axes`]
            A dictionary mapping CCD names to their corresponding Axes objects.

    """
    if cnams == "0" or len(mccd) == 1:
        cnams = list(mccd.keys())
    if isinstance(cnams, str):
        cnams = [cnams]

    # Create the figure and axes
    fig, axs = create_figure(
        cnams=cnams,
        ncols=ncols,
        region=region,
        title=title,
        figsize=figsize,
        position=position,
    )

    # Plot each CCD on its axes
    for cnam, ax in axs.items():
        # Subtract median from each window, if requested
        if msub:
            for wind in mccd[cnam].values():
                wind -= wind.median()

        # Get scale limits
        vmin, vmax = get_scale_limits(mccd[cnam], iset, ilim, region)

        # Plot the ccd
        pCcd(
            ax,
            mccd[cnam],
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            tlabel=f"CCD {cnam}",
        )

        # Plot the defects, if any
        if defects is not None:
            pCcdDefect(ax, defects[cnam] if cnam in defects else None)

    # Adjust the layout to prevent overlap
    fig.tight_layout()

    return fig, axs


class AnimationManager(ABC):
    """Abstract class to manage an animated matplotlib plot.

    Basic design comes from the example in the matplotlib documentation:
    https://matplotlib.org/stable/users/explain/animations/blitting.html
    """

    def __init__(self, fig):
        """Initialises the animated plot using the given figure."""
        self.fig = fig  # ensure we have a reference to the figure

        # default holding dict for animated artists
        self.artists = {}

        # register a callback to save the background when the figure is drawn
        self._bg = None
        self.fig.canvas.mpl_connect("draw_event", self.on_draw)

    @abstractmethod
    def draw_artists(self):
        """Draw all of the animated axes and artists.

        Must be implemented by subclasses to draw the specific artists
        using `self.fig.draw_artist()`.

        Note that includes the axes themselves, which are not automatically updated
        when created with `animated=True`.
        So you need to call `self.fig.draw_artist(ax)` for each axes you want to draw,
        in addition to any other artists.
        """
        pass

    def on_draw(self, event):
        """Callback to register with 'draw_event'."""
        if event is not None and event.canvas != self.fig.canvas:
            raise RuntimeError
        # save the background
        self._bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        # redraw the animated artists on top of the background
        self.draw_artists()

    def refresh(self, redraw=False):
        """Redraws the animated artists.

        Safe to call repeatedly to keep the window responsive even if the data has
        not changed (with `redraw=False`).
        If the data has changed (i.e. in `update()`), then call with `redraw=True`.
        """
        if self._bg is None:
            # no background saved yet, so just draw everything
            self.on_draw(None)
        elif redraw:
            # restore the background
            self.fig.canvas.restore_region(self._bg)
            # draw the animated artists on top of the background
            self.draw_artists()
            # redraw the canvas with the updated artists
            self.fig.canvas.blit(self.fig.bbox)

        # flush events to keep the window responsive
        self.fig.canvas.flush_events()

    @abstractmethod
    def update(self, *args, **kwargs):
        """Create or update the animated artists.

        Must be implemented by subclasses to update the data of the artists.

        Normal usage is to create the artists if they do not exist yet,
        and store references to them in the instance (e.g. in a dict or list called
        `self.artists`), and then on subsequent calls update the data of the artists
        using `set_data()` or similar methods.

        After updating the data, you should call `self.refresh(redraw=True)`
        to redraw the image artists.
        """
        pass

    def show(self):
        """Show the figure."""
        plt.show(block=False)
        plt.pause(0.1)

    @property
    def is_open(self):
        """Check if the figure is still open."""
        return plt.fignum_exists(self.fig.number)


def setup_plot(fig, disable_shortcuts=True):
    """Set up the plot for interactive use.

    These are various tweaks to make interactive plotting more user-friendly.
    """
    backend = mpl.get_backend().lower()

    # Start the interactive plot in pan mode, if not using the Tk backend
    # (commit 11a9b65 added the skip "to help target selection when using TkAgg",
    #  not sure why or if it's still necessary?)
    if backend != "tkagg":
        fig.canvas.manager.toolbar.pan()

    # Disable the default key shortcuts as we often use them for other things
    if disable_shortcuts:
        for key in [k for k in mpl.rcParams.keys() if k.startswith("keymap.")]:
            mpl.rcParams[key] = ""

    # Replace the default matplotlib cursors with a more precise pointer.
    from matplotlib.backend_bases import cursors

    if backend.startswith("qt"):
        from matplotlib.backends.backend_qt import cursord
    elif backend.startswith("gtk3"):
        from matplotlib.backends.backend_gtk3 import cursord
    else:
        return

    try:
        # Pan mode
        cursord[cursors.MOVE] = cursord[cursors.POINTER]
        # Zoom mode
        cursord[cursors.SELECT_REGION] = cursord[cursors.POINTER]
    except KeyError:
        pass
