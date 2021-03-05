import sys
import os

import numpy as np
import matplotlib as mpl

# re-configure the cursors: backend specific.
# aim to get rid of irritating 'hand' icon in
# favour of something pointier.

backend = mpl.get_backend()

if backend == "Qt4Agg" or "Qt5Agg":
    from matplotlib.backends.backend_qt5 import cursord as curs
elif backend == "GTK3agg":
    from matplotlib.backends.backend_gtk3 import cursord as curs
else:
    curs = None

if curs is not None:
    from matplotlib.backend_bases import cursors

    try:
        curs[cursors.HAND] = curs[cursors.POINTER]
        curs[cursors.WAIT] = curs[cursors.POINTER]
        curs[cursors.SELECT_REGION] = curs[cursors.POINTER]
        curs[cursors.MOVE] = curs[cursors.POINTER]
    except AttributeError:
        pass

import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = [
    "flagcloud",
]

COL_AP1 = 'g'
COL_AP2 = 'b'
COL_CLOUD = 'darkorange'
COL_JUNK = 'r'
COL_RAT = 'k'

########################################################################
#
# flagcloud -- to flag cloud-affected/ bad data in a hipercam reduce log
#
########################################################################


def flagcloud(args=None):
    """``flagcloud hlog aper1 aper2 ccd delta output``

    Interactive flagging of cloud-affected or otherwise bad points in a
    |hiper| log file. You either mark a range of times as cloudy,
    or individual points as junk. If you mark a time range, then *all*
    apertures of *all* CCDs will be flagged with the bitmask value
    CLOUDS. Individual points will be flagged as JUNK. Note that nothing
    is done to the data apart from changing the bitmask flags, so it is
    then up to you to test for these later on. It is also possible to
    flag individual points as CLOUDS but these will not propagate across
    CCDs so it is probably not advisable to do this.

    Junk points are marked red, cloudy points orange. OK aperture 1
    points are plotted green, aperture 2 blue, their ratio black. What
    is meant by 'junk' as opposed to 'cloud' is really down to the
    user. I tend to reserve junk for one-off points affected by bad
    cosmic rays and satellites, but in general it is probably sensible
    to think of junk as points you never want to see again versus
    clouds meaning data that you might want to mask or use down the
    line or perhaps just grey out in plots. Some genuine "cloudy" data
    will be so bad that it will be better flagged as junk however. You can
    flag the same point as both "cloud" and "junk", but "junk" is the
    stronger condition.

    Bitmasks propagate when data are combined so a point flagged junk in
    aperture 2 but not aperture 1 will be flagged junk in the ratio of 1
    divided by 2. You can also recover points in this routine; doing so will
    clear both their junk and/or cloud status. At the moment this is only
    possible on a point-by-point basis.

    Interaction is via the cursor and hitting specific keys. Common options are
    lower case; less common ones upper case. Hitting the X on the plot will abort
    without saving the results. 'q' to quite saves the results. 'h' will give some
    help on the options.

    Parameters:

      hlog : str
         ASCII log file, as produced by |reduce|.

      aper1 : str
         the name of first aperture to look at

      aper2 : str
         the name of second aperture to look at. The ratio aper1 / aper2
         will be plotted along with the two separately, scaled by their
         maximum, all in the same panel.

      ccd : str
         CCD(s) to plot, '0' for all. If not '0' then '1', '2', or even '1 2 3'
         are possible inputs (without the quotes). Note the space separation when multiple
         CCDs are specified. If you want to plot more than one CCD, then you will get multiple
         panels in the Y direction, but their X-axes are kept in lock step when panning or zooming.

      delta : float
         separation to use to space the plots in a given panel, each of which is normalised
         to 1. A value of 1 is recommended because then the second aperture should end with
         a typical level of 0, and any dips below this show the extent of the cloud. e.g. -0.9
         would suggest an approximate 90% loss of flux due to cloud, ignoring extinction.

      output : str
         name of modified version of the Hlog for output. Can overwrite the original if you dare.

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("hlog", Cline.LOCAL, Cline.PROMPT)
        cl.register("aper1", Cline.LOCAL, Cline.PROMPT)
        cl.register("aper2", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("delta", Cline.LOCAL, Cline.PROMPT)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        hlog = cl.get_value(
            "hlog", "hipercam ASCII log file", cline.Fname("hcam", hcam.LOG)
        )
        hlog = hcam.hlog.Hlog.read(hlog)
        cnams = sorted(hlog.keys())
        apnames = set()
        for cnam in cnams:
            apnames |= set(hlog.apnames[cnam])
        apnames = sorted(apnames)

        aper1 = cl.get_value("aper1", "first aperture", apnames[0], lvals=apnames)
        if len(apnames) > 1:
            aper2 = cl.get_value("aper2", "second aperture", apnames[-1], lvals=apnames)
        else:
            aper2 = None

        max_ccd = len(hlog)
        if max_ccd > 1:
            ccd = cl.get_value("ccd", "CCD(s) to plot [0 for all]", "0")
            if ccd == "0":
                ccds = cnams
            else:
                ccds = ccd.split()
                if set(ccds) <= set(cnams):
                    print(f'   selected CCDs ({ccds}) not amongst those in the reduce log ({cnams})')
                    return
        else:
            ccds = cnams

        delta = cl.get_value(
            "delta", "vertical separation between plots", 1, 0.
        )

        output = cl.get_value(
            "output",
            "name for output log file",
            cline.Fname("run", hcam.LOG, cline.Fname.NEW),
        )

    # Inputs obtained.
    print()

    # re-configure keyboard shortcuts to avoid otherwise confusing behaviour
    # quit_all does not seem to be universal, hence the try/except
    try:
        mpl.rcParams["keymap.back"] = ""
        mpl.rcParams["keymap.forward"] = ""
        mpl.rcParams["keymap.fullscreen"] = ""
        mpl.rcParams["keymap.grid"] = ""
        mpl.rcParams["keymap.home"] = ""
        mpl.rcParams["keymap.pan"] = ""
        mpl.rcParams["keymap.quit"] = ""
        mpl.rcParams["keymap.save"] = ""
        mpl.rcParams["keymap.pan"] = ""
        mpl.rcParams["keymap.save"] = ""
        mpl.rcParams["keymap.xscale"] = ""
        mpl.rcParams["keymap.yscale"] = ""
        mpl.rcParams["keymap.zoom"] = ""
    except KeyError:
        pass

    cnams, anams, plots = {}, {}, {}
    ax = None

    fig,axs = plt.subplots(len(ccds),1,sharex=True)
    if len(ccds) == 1:
        axs = [axs]

    # get the navigation toolbar. Go straight into pan mode where we
    # want to stay.
    toolbar = fig.canvas.manager.toolbar
    if backend != "TkAgg":
        toolbar.pan()

    ny = len(ccds)
    T0 = None
    for cnam, ax in zip(ccds,axs):

        # prep data
        a1 = hlog.tseries(cnam, aper1)
        a1.normalise()
        if T0 is None:
            T0 = a1.t[0]

        if aper2:
            a2 = hlog.tseries(cnam, aper2)
            a2.normalise()
            rat = a1 / a2
            a2.t -= T0
            rat.t -= T0
            rat += delta
            a2 -= delta

            rat.mplot(ax, COL_RAT, bitmask=hcam.JUNK | hcam.CLOUDS)
            rat.mplot(ax, COL_CLOUD, bitmask=hcam.CLOUDS, flagged=True)
            rat.mplot(ax, COL_JUNK, bitmask=hcam.JUNK, flagged=True)

            a2.mplot(ax, COL_AP2, bitmask=hcam.JUNK | hcam.CLOUDS)
            a2.mplot(ax, COL_CLOUD, bitmask=hcam.CLOUDS, flagged=True)
            a2.mplot(ax, COL_JUNK, bitmask=hcam.JUNK, flagged=True)

        a1.t -= T0
        a1.mplot(ax, COL_AP1, bitmask=hcam.JUNK | hcam.CLOUDS)
        a1.mplot(ax, COL_CLOUD, bitmask=hcam.CLOUDS, flagged=True)
        a1.mplot(ax, COL_JUNK, bitmask=hcam.JUNK, flagged=True)

        # store the plots needed to identify which point has been selected
        # along with names of apertures needed to flag points
        plots[cnam] = {
            "aper1" : aper1,
            "a1": a1,
            "aper2" : aper2,
            "a2": a2 if aper2 is not None else None,
            "rat": rat if aper2 is not None else None
        }

        # keep track of the CCD associated with each axes
        cnams[ax] = cnam

        # and the axes associated with each CCD
        anams[cnam] = ax
        ax.set_ylabel("CCD {:s}".format(cnam))

    ax.set_xlabel(f'Time [MJD - {T0}]')

    # create the picker
    picker = PickPoints(fig, hlog, cnams, anams, plots, T0, output)

    try:
        plt.tight_layout()
    except:
        pass

    PickPoints.action_prompt(False)

    # squeeze space a bit
    plt.subplots_adjust(hspace=0.1)

    # finally show stuff ....
    plt.xlabel("Time [MJD - {:.7f}]".format(T0))
    plt.show()


def nearest(x, y, lc, ax, fig, dmax=0.02, rmin=1.5):
    """Given an x, y location in an Axes ax, plotting a Tseries lc, this
    comes back with the index of the nearest point in the light curve,
    and the distance from it as a fraction of xwidth of the plot as a
    two-element tuple. It returns (None,None) if no point is found
    less than a fraction dmax of the x width from x,y. The nearest
    point must also be clearly the nearest point such that it is at
    least rmin times closer than the next nearest.

    Arguments::

      x : float
        X-position near point

      y : float
        Y-position near point

      lc : Tseries
        the light-curve to test against

      ax : Axes
        Axes instance

      fig : Figure
        containing Figure

      dmax : float
        maximum distance. If no point is closer than this as a fraction
        of the figure width, then nothing will be selected.

      rmin : float
        minimum ratio between second closest and closest distances. The best point
        must be clearly better.

    """

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()

    # select only points in range
    ok = (((lc.t > x1) & (lc.t < x2)) | ((lc.t < x1) & (lc.t > x2))) & (
        ((lc.y > y1) & (lc.y < y2)) | ((lc.y < y1) & (lc.y > y2))
    )

    if len(lc.t[ok]):
        indices = np.mgrid[0 : len(lc)]
        dsq = (width * (lc.t[ok] - x) / (x2 - x1)) ** 2 + (
            height * (lc.y[ok] - y) / (y2 - y1)
        ) ** 2
        imin = dsq.argmin()
        dmin = np.sqrt(dsq[imin])
        if len(dsq) > 1:
            # check that second nearest point is distinctly further away
            rest = np.mgrid[0 : len(dsq)] != imin
            dnext = np.sqrt(dsq[rest].min())
            if dnext < rmin * dmin:
                return (None, None)

        if dmin < width * dmax:
            return (indices[ok][imin], dmin)
        else:
            return (None, None)
    else:
        return (None, None)


class PickPoints:
    """This is where all the action occurs. A rather complicated matter
    of handling events. Note that standard terminal input with 'input'
    becomes impossible, explaining some of the weirdness. Effectively
    the class is used here to define a scope for variables that would
    otherwise be treated as globals
    """

    def __init__(self, fig, hlog, cnams, anams, plots, T0, oname):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect("key_press_event", self._keyPressEvent)
        self.hlog = hlog
        self.cnams = cnams
        self.anams = anams
        self.plots = plots
        self.T0 = T0
        self.oname = oname

        # flag to indicate the action we are in which requires extra input
        self._range_mode = False

    @staticmethod
    def action_prompt(cr):
        """Prompts user for an action. cr controls whether there is an intial
        carriage return or not.  It leaves the cursor at the end of the
        line. This appears in multiple places hence why it is a method.

        """
        if cr:
            print()

        print(
            "c(range), C(point), j(point), J(range), r(point), R(range), h(elp), q(uit): ",
            end="",
            flush=True,
        )

    def _keyPressEvent(self, event):
        """This is where we do the hard work. Every key press event is
        diverted to this method. It either takes an action based on
        the input, such as flagging a point as cloudy, or sometimes it
        causes a state change such that input is diverted to and
        accumulated in a buffer until 'enter' is hit.  The latter
        stage comes first.

        """

        if self._range_mode:

            if event.key == "q":
                # trap 'q' for quit
                self._range_mode = False
                print("no time range was defined")
                PickPoints.action_prompt(True)

            elif event.key == "c" or event.key == "J" or event.key == "R":
                # range mode.
                self._cnam = self.cnams[event.inaxes]
                self._axes = event.inaxes
                self._x = event.xdata
                self._y = event.ydata
                self._range(event.key)

        else:

            # standard mode action
            self._standard(event.key, event.xdata, event.ydata, event.inaxes)

    def _standard(self, key, x, y, axes):
        """Carries out the work needed when we are in the standard mode. Just
        pass through the value of the key pressed, the associated x, y
        position and the axes instance (all from the event) and this
        will handle the rest.

        """

        if axes is not None:

            # store information in attributes accessible to all methods, for
            # later access: name, the key hit, the axes instance, x, y
            self._cnam = self.cnams[axes]
            self._key = key
            self._axes = axes
            self._x = x
            self._y = y

            if key == "h":
                # help text
                print(key)
                print("""

Help on the actions available in 'flagcloud':

  c   : mark a time range over which all apertures of all CCDs are cloudy
  C   : mark a single point of one aperture of one CCD as cloudy [unusual]
  j   : mark a point of one aperture of a CCD as junk
  J   : mark a time range over which all apertures of all CCDs are junk
  r   : mark a point to clear any cloud and/or junk flag
  R   : mark a time range over which all apertures of all CCDs are cleared of junk and clouds
  h   : print this help text
  q   : quit 'flagcloud', and save the Hlog to disk.

Hitting X-on the gui or ctrl-C will abort without saving anything. 'c'. 'J' and 'R' all
require setting 2 points which should be from the same window.
"""
                      )

            elif key == "C" or key == "j" or key == "r":

                # selects a single point to flag or unflag as junk / cloud
                print(key)

                # short-hand
                cnam = self.cnams[axes]
                plot = self.plots[cnam]
                p1, p2, rat = plot['a1'], plot['a2'], plot['rat']

                i1, d1 = nearest(x, y, p1, axes, self.fig)
                if p2 is not None:
                    i2, d2 = nearest(x, y, p2, axes, self.fig)
                else:
                    i2 = None

                ap = None
                if i1 is not None and i2 is not None:
                    if 1.5*min(d1,d2) > max(d1,d2):
                        print('ambiguous choice: must be >1.5x closer a point in one aperture than to any other')

                    elif d1 < d2:
                        ap, ind = plot['aper1'], i1
                        col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_AP1
                        axes.errorbar(p1.t[i1],p1.y[i1],p1.ye[i1],fmt='.',color=col,zorder=100)
                        col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_RAT
                        axes.errorbar(rat.t[i1],rat.y[i1],rat.ye[i1],fmt='.',color=col,zorder=100)

                    elif d2 < d1:
                        ap, ind = plot['aper2'], i2
                        col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_AP2
                        axes.errorbar(p2.t[i2],p2.y[i2],p2.ye[i2],fmt='.',color=col,zorder=100)
                        col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_RAT
                        axes.errorbar(rat.t[i2],rat.y[i2],rat.ye[i2],fmt='.',color=col,zorder=100)

                elif i1 is not None:
                    ap, ind = plot['aper1'], i1
                    col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_AP1
                    axes.errorbar(p1.t[i1],p1.y[i1],p1.ye[i1],fmt='.',color=col,zorder=100)
                    if rat is not None:
                        col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_RAT
                        axes.errorbar(rat.t[i1],rat.y[i1],rat.ye[i1],fmt='.',color=col,zorder=100)

                elif i2 is not None:
                    ap, ind = plot['aper2'], i2
                    col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_AP2
                    axes.errorbar(p2.t[i2],p2.y[i2],p2.ye[i2],fmt='.',color=col,zorder=100)
                    col = COL_JUNK if key == 'j' else COL_CLOUD if key == 'C' else COL_RAT
                    axes.errorbar(rat.t[i2],rat.y[i2],rat.ye[i2],fmt='.',color=col,zorder=100)

                if ap is not None:
                    # propagate change into the hlog
                    flag = self.hlog[cnam][f'flag_{ap}'][ind]
                    self.hlog[cnam][f'flag_{ap}'][ind] = \
                        (flag | hcam.JUNK) if key == 'j' else \
                        (flag | hcam.CLOUDS) if key == 'c' else \
                        (flag & ~hcam.JUNK & ~hcam.CLOUDS)

                plt.draw()
                PickPoints.action_prompt(False)

            elif key == "c" or key == "J" or key == "R":
                # define a time range
                print(key)

                self._range_stage = 0
                self._range_mode = True
                self._range(key)

            elif key == "q":
                print(key)

                # quit and clear up
                plt.close()

                # old files are over-written at this point
                self.hlog.write(self.oname)
                print(f"\nHlog saved to {self.oname}.\n")

            elif (
                key == "shift"
                or key == "alt"
                or key == "control"
                or key == "pagedown"
                or key == "pageup"
            ):
                # trap some special keys to avoid irritating messages
                pass

            else:
                print(f'\nNo action is defined for key = "{key}"')
                PickPoints.action_prompt(False)

    def _range(self, key):
        """Defines time range
        """

        self._range_stage += 1

        if self._range_stage == 1:
            # first time through, store the CCD of the first point selected
            # (insist that this matches the second point later)
            self._range_cnam = self._cnam
            print(
                f"'{key}' to select the other end of the time range or 'q' to quit: ",
                end="", flush=True
            )
            self._t1 = self._x
            self._key = key

        else:

            # second time through. Check we are in the same CCD
            print(key)

            if self._cnam != self._range_cnam:
                # wrong CCD / axes: quit because this suggests not
                # doing right thing at all
                self._range_mode = False
                print("  *** must define time ranges on the same CCD / plot axis; nothing done")

            elif key != self._key:
                # hit different key; probably a mistake. Get another chance
                print(
                    f"  invalid response ['{key}']. Press '{self._key}' to select the other end of the time range or 'q' to quit: ",
                    end="", flush=True
                )
                self._range_stage -= 1

            else:
                self._range_mode = False
                # order the times
                t1, t2 = self._t1, self._x
                t1, t2 = (t1, t2) if t1 < t2 else (t2, t1)

                # mark points across all CCDs
                for cnam, plot in self.plots.items():
                    axes = self.anams[cnam]
                    p1, p2, rat = plot['a1'], plot['a2'], plot['rat']

                    col = COL_CLOUD if key == 'c' else \
                        COL_JUNK if key == 'J' else COL_AP1
                    mask = (p1.t > t1) & (p1.t < t2)
                    p1.mplot(axes,col,zorder=100,mask=mask)

                    if p2 is not None:
                        col = COL_CLOUD if key == 'c' else \
                            COL_JUNK if key == 'J' else COL_AP2
                        p2.mplot(axes,col,zorder=100,mask=mask)
                        col = COL_CLOUD if key == 'c' else \
                            COL_JUNK if key == 'J' else COL_RAT
                        rat.mplot(axes,col,zorder=100,mask=mask)

                # Propagate the changes into the hlog
                t1 += self.T0
                t2 += self.T0
                for cnam in self.hlog:
                    data = self.hlog[cnam]
                    apnames = self.hlog.apnames[cnam]
                    ts = self.hlog[cnam]['MJD']
                    mask = (ts > t1) & (ts < t2)
                    for ap in apnames:
                        if key == 'c':
                            # clouds
                            data[f'flag_{ap}'][mask] |= hcam.CLOUDS
                        elif key == 'J':
                            # it's junk
                            data[f'flag_{ap}'][mask] |= hcam.JUNK
                        elif key == 'R':
                            # clear junk and clouds
                            data[f'flag_{ap}'][mask] &= (~hcam.JUNK & ~hcam.CLOUDS)

                plt.draw()
                PickPoints.action_prompt(True)

