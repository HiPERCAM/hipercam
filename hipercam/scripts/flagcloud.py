import sys
import os

import numpy as np
import matplotlib as mpl

# re-configure the cursors: backend specific.
# aim to get rid of irritating 'hand' icon in
# favour of something pointier.

backend = mpl.get_backend()

if backend == 'Qt4Agg' or 'Qt5Agg':
    from matplotlib.backends.backend_qt5 import cursord as curs
elif backend == 'GTK3agg':
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

__all__ = ['flagcloud',]

#############################################
#
# setaper -- defines apertures given an image
#
#############################################

def flagcloud(args=None):
    """``flagcloud hlog aper1 aper2 ccd output``

    Interactive flagging of cloud-affected or bad points in a |hipercam| log file.

    Parameters:

      hlog : string
         ASCII log file.

      aper1 : string
         the name of first aperture to look at

      aper2 : string
         the name of second aperture to look at. The ratio aper1 / aper2
         will be plotted along with the two separately, scaled by their
         maximum, all in the same panel.

      ccd : string
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction.

      output : string
         name of modified version of the Hlog for output.
    """

    command, args = utils.script_args(args)

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('hlog', Cline.LOCAL, Cline.PROMPT)
        cl.register('aper1', Cline.LOCAL, Cline.PROMPT)
        cl.register('aper2', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('output', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        hlog = cl.get_value(
            'hlog', 'hipercam ASCII log file',
            cline.Fname('hcam', hcam.LOG)
        )
        hlog = hcam.hlog.Hlog.read(hlog)

        aper1 = cl.get_value('aper1', 'first aperture', '2')
        aper2 = cl.get_value('aper2', 'second aperture', '3')

        max_ccd = len(hlog)
        if max_ccd > 1:
            ccd = cl.get_value('ccd', 'CCD(s) to plot [0 for all]', '0')
            if ccd == '0':
                ccds = list(hlog.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(hlog.keys())

        output = cl.get_value(
            'output', 'name for output log file',
            cline.Fname('run', hcam.LOG, cline.Fname.NEW)
        )

    # Inputs obtained.

    # re-configure keyboard shortcuts to avoid otherwise confusing behaviour
    # quit_all does not seem to be universal, hence the try/except
    try:
        mpl.rcParams['keymap.back'] = ''
        mpl.rcParams['keymap.forward'] = ''
        mpl.rcParams['keymap.fullscreen'] = ''
        mpl.rcParams['keymap.grid'] = ''
        mpl.rcParams['keymap.home'] = ''
        mpl.rcParams['keymap.pan'] = ''
        mpl.rcParams['keymap.quit'] = ''
        mpl.rcParams['keymap.save'] = ''
        mpl.rcParams['keymap.pan'] = ''
        mpl.rcParams['keymap.save'] = ''
        mpl.rcParams['keymap.xscale'] = ''
        mpl.rcParams['keymap.yscale'] = ''
        mpl.rcParams['keymap.zoom'] = ''
    except KeyError:
        pass

    cnams, anams = {}, {}
    plots = {}
    ax = None

    fig = plt.figure()

    ny = len(ccds)
    T0 = None
    for n, cnam in enumerate(ccds):

        if ax is None:
            axes = ax = fig.add_subplot(ny, 1, n+1)
        else:
            axes = fig.add_subplot(ny, 1, n+1, sharex=ax)

        # prep data
        a1 = hlog.tseries(cnam,aper1)
        a2 = hlog.tseries(cnam,aper1)
        rat = (a1/a2).normalise()
        a1 /= np.percentile(a1.y,99)
        a2 /= np.percentile(a2.y,99)

        if T0 is None:
            T0 = a1.t[0]

        a1.t -= T0
        a2.t -= T0
        rat.t -= T0

        # three vector plots
        (rat+0.1).mplot(plt,'k')
        a1.mplot(plt,'g')
        (a2-0.1).mplot(plt,'b')

        # store the plots needed to identify which point has been selected
        plots[cnam] = {'a1' : a1, 'a2' : a2-0.1, 'rat' : rat+0.1}

        # keep track of the CCD associated with each axes
        cnams[axes] = cnam

        # and the axes associated with each CCD
        anams[cnam] = axes
        plt.ylabel('CCD {:s}'.format(cnam))

    # create the picker
    picker = PickPoints(fig, hlog, cnams, anams, plots, output)

    try:
        plt.tight_layout()
    except:
        pass

    PickPoints.action_prompt(False)

    # squeeze space a bit
    plt.subplots_adjust(hspace=0.1)

    # finally show stuff ....
    plt.xlabel('Time [MJD - {:.7f}]'.format(T0))
    plt.show()

def nearest(x, y, lc, ax, fig, dmax=0.02, rmin=1.5):
    """Given an x, y location in an Axes ax, plotting a Tseries lc, this
    comes back with the index of the nearest point in the light curve,
    and the distance from it as a fraction of xwidth of the plot as a
    two-element tuple. It returns (None,None) if no point is found
    less than a fraction dlim of the x width from x,y. 

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
        must be clearly better

    """

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()

    # select only points in range
    ok = (((lc.t > x1) & (lc.t < x2)) | ((lc.t < x1) & (lc.t > x2))) & \
         (((lc.y > y1) & (lc.y < y2)) | ((lc.y < y1) & (lc.y > y2)))

    if len(lc.t[ok]):
        indices = np.mgrid[0:len(lc)]
        dsq = (width*(lc.t[ok]-x)/(x2-x1))**2 + (height*(lc.y[ok]-y)/(y2-y1))**2
        imin = dsq.argmin()
        dmin = np.sqrt(dsq[imin])
        if len(dsq) > 1:
            rest = (np.mgrid[0:len(dsq)] != imin)
            dnext = np.sqrt(dsq[rest].min())
            if dnext < rmin*dmin:
                return (None,None)

        if dmin < width*dmax:
            print(indices.shape, ok.shape, imin, dmin)
            return (indices[ok][imin], dmin)
        else:
            return (None,None)
    else:
        return (None,None)

class PickPoints:
    """This is where all the action occurs. A rather complicated matter
    of handling events. Note that standard terminal input with 'input'
    becomes impossible, explaining some of the weirdness. Effectively
    the class is used here to define a scope for variables that would
    otherwise be treated as globals

    """

    ADD_PROMPT = "enter a label for the aperture, '!' to abort: "

    def __init__(self, fig, hlog, cnams, anams, plots, oname):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect('key_press_event', self._keyPressEvent)
        self.hlog = hlog
        self.cnams = cnams
        self.anams = anams
        self.plots = plots
        self.oname = oname

        # then mutually exclusive flags to indicate the action we are in
        # for actions that require extra input. We are not in these at the
        # start so we set them False
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
            'c(loud), j(unk), J(unk range), h(elp), r(estore), q(uit): ',
            end='', flush=True
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

            if event.key == 'q':
                # trap 'q' for quit
                self._link_mode = False
                print('no extra aperture added')
                PickPoints.action_prompt(True)

            elif event.key == 's':

                # range mode.
                self._cnam = self.cnams[event.inaxes]
                self._axes = event.inaxes
                self._x = event.xdata
                self._y = event.ydata
                self._range()

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

            if key == 'h':
                # help text
                print(key)
                print("""

Help on the actions available in 'flagcloud':

  c(loud)    : mark time range on all CCDs as cloud
  h(elp)     : print this help text
  j(unk)     : mark a point on one aperture of CCD as junk
  J(unk)     : mark a range on one CCD as junk
  r(estore)  : restore a point
  q(uit)     : quit 'flagcloud' and save the Hlog to disk
""")

            elif key == 'j':

                # need to select a point that we define as junk
                print(key)

                cnam = self.cnams[axes]
                plot = self.plots[cnam]
                i1, d1 = nearest(x, y, plot['a1'], axes, self.fig)
                i2, d2 = nearest(x, y, plot['a2'], axes, self.fig)
                print(i1, d1, i2, d2)
                if i1 is not None and i2 is not None:
                    if d1 < d2:
                        print('i1 =',i1,plot['a2'].t[i1],plot['a2'].y[i1])
                        plot['a1'].mask[i1] |= hcam.JUNK
                        axes.plot(plot['a1'].t[i1], plot['a1'].y[i1],'.r',ms=20)
                        axes.errorbar(
                            plot['a1'].t[i1], plot['a1'].y[i1], 
                            plot['a1'].ye[i1],
                            fmt='.r',zorder=100
                        )
                    else:
                        print('i2 =',i2,plot['a2'].t[i2],plot['a2'].y[i2])
                        plot['a2'].mask[i2] |= hcam.JUNK
                        axes.plot(plot['a2'].t[i2], plot['a2'].y[i2],'.r',ms=20)
                        axes.errorbar(
                            plot['a2'].t[i2], plot['a2'].y[i2], 
                            plot['a2'].ye[i2],
                            fmt='.r',zorder=100
                        )

                elif i1 is not None:
                    print('i1 =',i1,plot['a2'].t[i1],plot['a2'].y[i1])
                    plot['a1'].mask[i1] |= hcam.JUNK
                    axes.plot(plot['a1'].t[i1], plot['a1'].y[i1],'.r',ms=20)
                    axes.errorbar(
                        plot['a2'].t[i2], plot['a2'].y[i2], plot['a2'].ye[i2],
                        fmt='.r',zorder=100
                    )

                elif i2 is not None:
                    print('i2 =',i2,plot['a2'].t[i2],plot['a2'].y[i2])
                    plot['a2'].mask[i2] |= hcam.JUNK
                    axes.errorbar(
                        plot['a2'].t[i2], plot['a2'].y[i2], plot['a2'].ye[i2],
                        fmt='.r',zorder=100
                    )

            elif key == 'c':
                # define range of cloudy data
                print(key)
                print('not implemented')

            elif key == 'e':
                # add extra target pixels to an aperture
                print(key)
                print('not implemented')

                # switch to extra mode
#                self._extra_mode = True
#                self._extra_stage = 0
#                self._extra()

            elif key == 'q':
                print(key)

                # quit and clear up
                plt.close()

                # old files are over-written at this point
                self.hlog.write(self.oname)
                print('\nHlog saved to {:s}.\nBye'.format(self.oname))

            elif key == 'r':
                print(key)
                print('not implemented')

            elif key == 'C':
                print(key)
                print('not implemented')

            elif key == 'shift' or key == 'alt' or key == 'control' or \
                 key == 'pagedown' or key == 'pageup':
                # trap some special keys to avoid irritating messages
                pass

            else:
                print('\nNo action is defined for key = "{:s}"'.format(key))
                PickPoints.action_prompt(False)

    def _select(self):
        """Selects a range
        """


        self._select_stage += 1

        if self._select_stage == 1:

            # first end of range
            self._select_cnam = self._cnam
            print(" 's' to select other end of range ['q' to quit]")

        else:
            # second time through. Check we are in the same CCD
            # but on a different aperture first change the link status
            self._select_mode = False

            # then see if we can make a valid link.
            if self._cnam != self._select_cnam:
                print('  *** cannot select across CCDs; no range set')

            else:
                # add link to the first aperture
                self._link_aper.set_link(apnam)

                # re-plot new version, over-writing plot objects
                self.pobjs[self._cnam][self._link_apnam] = hcam.mpl.pAper(
                    self._axes, self._link_aper, self._link_apnam,
                    self.mccdaper[self._link_cnam])
                plt.draw()

                print('  linked aperture {:s} to aperture {:s}'
                      ' in CCD {:s}'.format(
                          self._link_apnam, apnam, self._link_cnam))
                PickPoints.action_prompt(True)


