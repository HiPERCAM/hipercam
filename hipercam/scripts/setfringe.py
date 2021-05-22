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
from hipercam import cline, utils, fringe
from hipercam.cline import Cline

__all__ = [
    "setfringe",
]

#############################################
#
# setfringe -- defines pairs of points for fringe measurement
#
#############################################


def setfringe(args=None):
    """``setfringe fmap fringe ccd [width height] nx [cmap] iset (ilo ihi
    | plo phi)``

    Interactive definition of CCD fringe pairs. The idea is to place a large
    number of such pairs at the peaks and valleys of a fringe map file. To
    defringe a file, the differences of intensity for each pair will be evaluated
    in the data and the fringe map and the ratio of these differences taken. The
    median of these ratios will be used to scale the fringe map which will then be
    subtracted. The use of a median allows the odd ratio to be ruined by stars,
    as long as a large number of pairs are defined.

    Parameters:

      fmap : str
         name of an MCCD file showing fringes with objects removed, e.g.
         as produced by 'makefringe'

      fringe : str
         the name of a file of fringe pairs. If it exists it will be
         read so that more fringe pairs can be added to it. If it does
         not exist, it will be created on exiting the routine. The
         files are in a fairly readable / editable text format

      ccd : str
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      width : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

      nx : int
         number of panels across to display, prompted if more than one CCD is
         to be plotted.

      cmap : str [hidden]
         The colour map to use. "Greys" is the usual; "Greys_r" reverses it.
         There are many others; typing an incorrect one will give a list. "none"
         for matplotlib default.

      hsbox  : int
         half-width in binned pixels of stats box as offset from central pixel
         hsbox = 1 gives a 3x3 box; hsbox = 2 gives 5x5 etc. This is used by
         the "show" option when setting FringePair

      iset : str [single character]
         determines how the intensities are determined. There are three
         options: 'a' for automatic simply scales from the minimum to the
         maximum value found on a per CCD basis. 'd' for direct just takes two
         numbers from the user. 'p' for percentile dtermines levels based upon
         percentiles determined from the entire CCD on a per CCD bais.

      ilo : float [if iset=='d']
         lower intensity level

      ihi : float [if iset=='d']
         upper intensity level

      plo : float [if iset=='p']
         lower percentile level

      phi : float [if iset=='p']
         upper percentile level

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("fmap", Cline.LOCAL, Cline.PROMPT)
        cl.register("fringe", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
        cl.register("hsbox", Cline.GLOBAL, Cline.HIDE)
        cl.register("iset", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ilo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ihi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("plo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("phi", Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        mccd = cl.get_value("fmap", "fringe map frame", cline.Fname("hcam", hcam.HCAM))
        mccd = hcam.MCCD.read(mccd)

        fpair = cl.get_value(
            "fringe",
            "name of fringe pair file",
            cline.Fname("fringe", hcam.FRNG, exist=False),
        )

        if os.path.exists(fpair):
            # read in old fringe pairs
            mccd_fpair = fringe.MccdFringePair.read(fpair)
            print(f"Loaded existing file = {fpair}")
        else:
            # create empty container
            mccd_fpair = fringe.MccdFringePair()
            print(
                "No file called {:s} exists; " "will create from scratch".format(fpair)
            )

        max_ccd = len(mccd)
        if max_ccd > 1:
            ccd = cl.get_value("ccd", "CCD(s) to plot [0 for all]", "0")
            if ccd == "0":
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(mccd.keys())

        width = cl.get_value("width", "plot width (inches)", 0.0)
        height = cl.get_value("height", "plot height (inches)", 0.0)

        nxdef = cl.get_default("nx", 3)

        # number of panels in X
        if len(ccds) > 1:
            nxdef = min(len(ccds), nxdef)
            cl.set_default("nx", nxdef)
            nx = cl.get_value("nx", "number of panels in X", 3, 1)
        else:
            nx = 1

        cmap = cl.get_value("cmap", "colour map to use ['none' for mpl default]", "Greys")
        cmap = None if cmap == "none" else cmap

        hsbox = cl.get_value("hsbox", "half-width of stats box (binned pixels)", 2, 1)
        iset = cl.get_value(
            "iset",
            "set intensity a(utomatically)," " d(irectly) or with p(ercentiles)?",
            "a",
            lvals=["a", "A", "d", "D", "p", "P"],
        )
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == "d":
            ilo = cl.get_value("ilo", "lower intensity limit", 0.0)
            ihi = cl.get_value("ihi", "upper intensity limit", 1000.0)
        elif iset == "p":
            plo = cl.get_value(
                "plo", "lower intensity limit percentile", 5.0, 0.0, 100.0
            )
            phi = cl.get_value(
                "phi", "upper intensity limit percentile", 95.0, 0.0, 100.0
            )

        nxmax, nymax = 0, 0
        for cnam in ccds:
            nxmax = max(nxmax, mccd[cnam].nxtot)
            nymax = max(nymax, mccd[cnam].nytot)

        # might be worth trying to improve this at some point
        xlo, xhi, ylo, yhi = 0, nxmax + 1, 0, nymax + 1

    # Inputs obtained.

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

    # start plot
    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width, height))
    else:
        fig = plt.figure()

    # get the navigation toolbar. Go straight into pan mode
    # where we want to stay.
    toolbar = fig.canvas.manager.toolbar
    toolbar.pan()

    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    # we need to store some stuff
    ax = None
    cnams = {}
    anams = {}

    # this is a container for all the objects used to plot FringePairs to
    # allow deletion.
    pobjs = {}

    for n, cnam in enumerate(ccds):
        if ax is None:
            axes = ax = fig.add_subplot(ny, nx, n + 1)
            axes.set_aspect("equal", adjustable="box")
            axes.set_xlim(xlo, xhi)
            axes.set_ylim(ylo, yhi)
        else:
            axes = fig.add_subplot(ny, nx, n + 1, sharex=ax, sharey=ax)
            axes.set_aspect("equal", adjustable="box")

        hcam.mpl.pCcd(
            axes, mccd[cnam], iset, plo, phi, ilo, ihi, f"CCD {cnam}", cmap=cmap
        )

        # keep track of the CCDs associated with each axes
        cnams[axes] = cnam

        # and axes associated with each CCD
        anams[cnam] = axes

        if cnam in mccd_fpair:
            # plot any pre-existing FringePairs, keeping track of
            # the plot objects
            pobjs[cnam] = hcam.mpl.pCcdFringePair(axes, mccd_fpair[cnam])

        else:
            # add in an empty CcdFringePair for any CCD not already present
            mccd_fpair[cnam] = fringe.CcdFringePair()

            # and an empty container for any new plot objects
            pobjs[cnam] = {}

    # create the FringePair picker (see below for class def)
    picker = PickFringePair(
        mccd, cnams, anams, toolbar, fig, mccd_fpair, fpair, hsbox, pobjs
    )

    picker.action_prompt(False)

    # squeeze space a bit
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    # finally show stuff ....
    plt.show()


# the next class is where all the action occurs. A rather complicated matter
# of handling events. note that standard terminal input with 'input' becomes
# impossible, explaining some of the weirdness. Effectively the class is used
# here to define a scope for variables that would otherwise be treated as globals


class PickFringePair:
    """Class to pick fringe pairs
    """

    def __init__(
        self, mccd, cnams, anams, toolbar, fig, mccd_fpair, fringenam, hsbox, pobjs
    ):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect("key_press_event", self._keyPressEvent)
        self.mccd = mccd
        self.cnams = cnams
        self.anams = anams
        self.toolbar = toolbar
        self.mccd_fpair = mccd_fpair
        self.fringenam = fringenam
        self.hsbox = hsbox
        self.pobjs = pobjs

        # then mutually exclusive flags to indicate the action we are in
        # for actions that require extra input. We are not in these at the
        # start so we set them False
        self._mid_pair = False

    def action_prompt(self, cr):
        """Prompts user for an action. cr controls whether there is an initial
        carriage return or not.  It leaves the cursor at the end of the
        line. This appears in multiple places hence why it is a method.

        """
        if cr:
            print()

        print(
            "a(dd), d(elete), h(elp), s(how), q(uit): ",
            end="",
            flush=True,
        )

    def _keyPressEvent(self, event):
        """This is where we do the hard work. Every key press event is
        diverted to this method. It either takes an action based on
        the input, such as removing a fringe pair, or sometimes it
        causes a state change to get other input.

        """

        if self._mid_pair:

            if event.key == "q":
                self._mid_pair = False
                print("pair not added")
                self.action_prompt(True)
                return

            axes = event.inaxes
            if axes is not None:
                self._cnam = self.cnams[axes]
                self._key = event.key
                self._axes = axes
                self._x = event.xdata
                self._y = event.ydata
                self._pair()

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

            # store information in attributes accessible to all
            # methods, for later access: name, the key hit, the axes
            # instance, x, y
            self._cnam = self.cnams[axes]
            self._key = key
            self._axes = axes
            self._x = x
            self._y = y

            if key == "h":
                # help text
                print(key)
                print(
                    """

Help on the actions available in 'setfringe':

  a(dd)      : add a fringe pair
  d(elete)   : delete a fringe pair
  h(elp)     : print this help text
  s(how)     : show image values
  q(uit)     : quit 'setfringe' and save the fringe pairs to disk

Hitting 'd' will delete the fringe pair nearest to the cursor, as long
as it is close enough (< 10 pixels)
"""
                )

            elif key == "a":

                # add a fringe pair
                print(key)

                # Try to calculate the largest number, label the new
                # FringePair with one more
                high = 0
                for frng in self.mccd_fpair[self._cnam]:
                    try:
                        high = max(high, int(frng))
                    except ValueError:
                        pass

                self._buffer = str(high + 1)
                self._mid_pair = True
                self._pair_stage = 0
                self._pair()

            elif key == "d":
                # delete an fringe pair
                print(key)
                self._delete()

            elif key == "s":
                # show some values
                print(key)
                self._show()

            elif key == "q":
                print(key)
                # quit and clear up
                plt.close()

                # old files are over-written at this point
                self.mccd_fpair.write(self.fringenam)
                print(f"\nFringePairs saved to {self.fringenam}.\nBye")

            elif key == "enter":
                self.action_prompt(True)

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
                self.action_prompt(False)

    def _pair(self):
        """Once all set to add a pair, this routine actually carries out the
        necessary operations

        """

        self._pair_stage += 1

        if self._pair_stage == 1:

            wnam = self.mccd[self._cnam].inside(self._x, self._y, 0)
            if wnam is None:
                self._line_mode = False
                print("  cannot set pairs outside windows")
                self.action_prompt(True)

            else:

                # store the CCD, window, and the first x,y position
                self._first_cnam = self._cnam
                self._first_wnam = wnam
                self._first_x = self._x
                self._first_y = self._y

                # prompt stage 2
                print(" second point: a(dd) or q(uit)")

        elif self._pair_stage == 2:

            wnam = self.mccd[self._cnam].inside(self._x, self._y, 0)
            self._mid_pair = False
            if wnam is None:
                print("  cannot set pairs outside windows")
                self.action_prompt(True)

            elif wnam != self._first_wnam:
                print("  cannot set pairs across different windows")
                self.action_prompt(True)

            elif self._cnam != self._first_cnam:
                print("  cannot set pairs across different CCDs")
                self.action_prompt(True)

            else:
                # add new pair
                frng = fringe.FringePair(
                    self._first_x,
                    self._first_y,
                    self._x,
                    self._y
                )
                self.mccd_fpair[self._cnam][self._buffer] = frng

                # add fringe pair to the plot, store plot objects
                self.pobjs[self._cnam][self._buffer] = \
                    hcam.mpl.pFringePair(self._axes, frng)

                # make sure it appears
                plt.draw()

                # let user know what has happened
                print(
                    (
                        f"added fringe pair to CCD {self._cnam} "
                        f"from x1,y1 = {self._first_x:.2f},{self._first_y:.2f} "
                        f"to x2,y2 = {self._x:.2f},{self._y:.2f}"
                    )
                )
                self.action_prompt(True)

    def _show(self):
        """
        Prints stats on pixels around selected place
        """

        # search for enclosing window, print stats
        wnam, wind = utils.print_stats(
            self.mccd[self._cnam], self._cnam,
            self._x, self._y, self.hsbox, False
        )
        if wnam is None:
            print('  must hit "s" inside a window')

        self.action_prompt(True)

    def _delete(self):
        """This deletes the nearest fringe pair to the currently selected
        position, if it is near enough (within 10 pixels)

        """

        # first see if there is a FringePair near enough the selected position
        frng, frngnam, dmin = self._find_fringe()

        if dmin is not None and dmin < 10:

            # delete plot objects
            for pobj in self.pobjs[self._cnam][frngnam]:
                pobj.remove()

            # and containers
            del self.pobjs[self._cnam][frngnam]
            del self.mccd_fpair[self._cnam][frngnam]

            # update plot
            plt.draw()
            print(f'  deleted fringe pair "{frngnam}"')

        else:
            print("  no fringe pair near enough the cursor position for deletion")

        self.action_prompt(True)

    def _find_fringe(self):
        """Finds the nearest FringePair to the currently selected position,

        It returns (fringe, fringenam, dmin) where fringe is the
        FringePair, fringenam its label, and dmin is the minimum
        distance. These are all returned as None if no suitable
        FringePair is found.

        """

        dmin = None
        fringemin = None
        dnmin = None
        for fringenam, frng in self.mccd_fpair[self._cnam].items():
            dist = frng.dist(self._x, self._y)
            if dmin is None or dist < dmin:
                dmin = dist
                fringemin = frng
                dnmin = fringenam

        return (fringemin, dnmin, dmin)
