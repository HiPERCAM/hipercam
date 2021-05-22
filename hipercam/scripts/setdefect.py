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
from hipercam import cline, utils, defect
from hipercam.cline import Cline

__all__ = [
    "setdefect",
]

#############################################
#
# setdefect -- defines CCD defects given an image
#
#############################################


def setdefect(args=None):
    """``setdefect mccd defect ccd [width height] nx msub [cmap] ffield
    hsbox iset (ilo ihi | plo phi)``

    Interactive definition of CCD defects. This is a matplotlib-based routine
    allowing you to define defects using the cursor.

    Parameters:

      mccd   : string
         name of an MCCD file, as produced by e.g. 'grab'

      defect : string
         the name of a defect file. If it exists it will be read so that
         defects can be added to it. If it does not exist, it will be
         created on exiting the routine. The defect files are in a fairly
         readable / editiable text format

      ccd : string
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

      msub   : bool
         True/False to subtract median from each window before scaling

      cmap : str [hidden]
         The colour map to use. "Greys" is the usual; "Greys_r" reverses it.
         There are many others; typing an incorrect one will give a list. "none"
         for matplotlib default.

      ffield : bool
         If True, all defects will be assumed to be flat-field or poor
         charge transfer defects as opposed to hot pixels. The latter are
         best set from dark frames, and have a different impact than the
         first two types in that they are worst for faint targets. Hot pixels
         and flat-field defects are shown with the same colours for moderate
         and severe, but different symbols (filled circles for flat-field
         defects, stars for hot pixels). If you say no to add hot pixels,
         the line defect option is not available.

      hsbox : int
         half-width in binned pixels of stats box as offset from central pixel
         hsbox = 1 gives a 3x3 box; hsbox = 2 gives 5x5 etc. This is used by
         the "show" option when setting defects.

      iset   : string [single character]
         determines how the intensities are determined. There are three
         options: 'a' for automatic simply scales from the minimum to the
         maximum value found on a per CCD basis. 'd' for direct just takes two
         numbers from the user. 'p' for percentile dtermines levels based upon
         percentiles determined from the entire CCD on a per CCD bais.

      ilo    : float [if iset=='d']
         lower intensity level

      ihi    : float [if iset=='d']
         upper intensity level

      plo    : float [if iset=='p']
         lower percentile level

      phi    : float [if iset=='p']
         upper percentile level

    There are a few conveniences to make setdefect easier:

      1. The plot is initialised in pan mode whereby you can move around and
         scale using the left and right mouse buttons.

      2. All input is accomplished with the keyboard; the mouse buttons are
         only for navigating the image.

      3. The label input can be switched between sequential numerical,
         single- and multi-character input ('linput').

    Various standard keyboard shortcuts (e.g. 's' to save) are disabled as
    they just confuse things and are of limited use in setdefect in any case.

    Some aspects of the usage of matplotlib in setdefect are tricky. It is
    possible that particular 'backends' will cause problems. I have tested
    this with Qt4Agg, Qt5agg and GTK3Agg. One aspect is the cursor icon in pan
    mode is a rather indistinct hand where one can't tell what is being
    pointed at. I have therefore suppressed this, but only for the tested
    backends. Others would need require further investigation.

    NB At the end of this routine, it re-orders the defects so that the severe
    ones follows the moderates. This helps emphasize the severe ones over the
    moderates when running rtplot.

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("mccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("defect", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("msub", Cline.GLOBAL, Cline.PROMPT)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
        cl.register("ffield", Cline.GLOBAL, Cline.PROMPT)
        cl.register("hsbox", Cline.GLOBAL, Cline.HIDE)
        cl.register("iset", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ilo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ihi", Cline.GLOBAL, Cline.PROMPT)
        cl.register("plo", Cline.GLOBAL, Cline.PROMPT)
        cl.register("phi", Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        mccd = cl.get_value("mccd", "frame to plot", cline.Fname("hcam", hcam.HCAM))
        mccd = hcam.MCCD.read(mccd)

        dfct = cl.get_value(
            "defect",
            "name of defect file",
            cline.Fname("defect", hcam.DFCT, exist=False),
        )

        if os.path.exists(dfct):
            # read in old defects
            mccd_dfct = defect.MccdDefect.read(dfct)
            print("Loaded existing file = {:s}".format(dfct))
        else:
            # create empty container
            mccd_dfct = defect.MccdDefect()
            print(
                "No file called {:s} exists; " "will create from scratch".format(dfct)
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

        # define the display intensities
        msub = cl.get_value("msub", "subtract median from each window?", True)

        cmap = cl.get_value("cmap", "colour map to use ['none' for mpl default]", "Greys")
        cmap = None if cmap == "none" else cmap

        ffield = cl.get_value("ffield", "flat field defects? [else hot pixels]", True)

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

    # this is a container for all the objects used to plot Defects to allow
    # deletion. 
    pobjs = {}

    for n, cnam in enumerate(ccds):
        if ax is None:
            axes = ax = fig.add_subplot(ny, nx, n + 1)
            axes.set_aspect("equal", adjustable="box")
            axes.set_xlim(xlo, xhi)
            axes.set_ylim(ylo, yhi)
        else:
            axes = fig.add_subplot(ny, nx, n + 1, sharex=ax, sharey=ax)
#            axes.set_aspect("equal", adjustable="datalim")
            axes.set_aspect("equal", adjustable="box")

        if msub:
            # subtract median from each window
            for wind in mccd[cnam].values():
                wind -= wind.median()

        hcam.mpl.pCcd(
            axes, mccd[cnam], iset, plo, phi, ilo, ihi, f"CCD {cnam}", cmap=cmap
        )

        # keep track of the CCDs associated with each axes
        cnams[axes] = cnam

        # and axes associated with each CCD
        anams[cnam] = axes

        if cnam in mccd_dfct:
            # plot any pre-existing Defects, keeping track of
            # the plot objects
            pobjs[cnam] = hcam.mpl.pCcdDefect(axes, mccd_dfct[cnam])

        else:
            # add in an empty CcdDefect for any CCD not already present
            mccd_dfct[cnam] = defect.CcdDefect()

            # and an empty container for any new plot objects
            pobjs[cnam] = {}

    # create the Defect picker (see below for class def)
    picker = PickDefect(
        mccd, cnams, anams, toolbar, fig, mccd_dfct, dfct, ffield, hsbox, pobjs
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


class PickDefect:
    """Class to pick defects
    """

    def __init__(
        self, mccd, cnams, anams, toolbar, fig, mccd_dfct, dfctnam, ffield, hsbox, pobjs
    ):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect("key_press_event", self._keyPressEvent)
        self.mccd = mccd
        self.cnams = cnams
        self.anams = anams
        self.toolbar = toolbar
        self.mccd_dfct = mccd_dfct
        self.dfctnam = dfctnam
        self.ffield = ffield
        self.hsbox = hsbox
        self.pobjs = pobjs

        # then mutually exclusive flags to indicate the action we are in
        # for actions that require extra input. We are not in these at the
        # start so we set them False
        self._line_mode = False

    def action_prompt(self, cr):
        """Prompts user for an action. cr controls whether there is an intial
        carriage return or not.  It leaves the cursor at the end of the
        line. This appears in multiple places hence why it is a method.

        """
        if cr:
            print()

        if self.ffield:
            print(
                "d(elete), h(elp), l(ine), m(odest), n(asty), s(how), q(uit): ",
                end="",
                flush=True,
            )
        else:
            print(
                "d(elete), h(elp), m(odest), n(asty), s(how), q(uit): ",
                end="",
                flush=True,
            )

    def _keyPressEvent(self, event):
        """
        This is where we do the hard work. Every key press event is diverted
        to this method. It either takes an action based on the input, such as
        removing a defect, or sometimes it causes a state change to get other
        input.
        """

        if self._line_mode:

            if event.key == "q":
                self._line_mode = False
                print("no line defect added")
                self.action_prompt(True)
                return

            elif event.key == "m":
                self._severity = defect.Severity.MODERATE

            elif event.key == "n":
                self._severity = defect.Severity.SEVERE

            axes = event.inaxes
            if axes is not None:
                self._cnam = self.cnams[axes]
                self._key = event.key
                self._axes = axes
                self._x = event.xdata
                self._y = event.ydata
                self._line()

        else:
            # standard mode action
            self._standard(event.key, event.xdata, event.ydata, event.inaxes)

    def _standard(self, key, x, y, axes):
        """Carries out the work needed when we are in the standard mode. Just pass
        through the value of the key pressed, the associated x, y position and
        the axes instance (all from the event) and this will handle the rest.

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
                print(
                    """

Help on the actions available in 'setdefect':

  d(elete)   : delete a defect
  h(elp)     : print this help text
  l(ine)     : add a line defect (flat-field defects, not hot pixels)
  m(odest)   : add a moderate-level point defect
  n(asty)    : add a severe-level point defect
  s(how)     : show image values
  q(uit)     : quit 'setdefect' and save the defects to disk

Hitting 'd' will delete the defect nearest to the cursor, as long as it is
close enough (< 10 pixels)
"""
                )

            elif key == "m" or key == "n":

                # add a point defect
                print(key)

                # Try to calculate the largest number, label the new Defect
                # with one more
                high = 0
                for dfct in self.mccd_dfct[self._cnam]:
                    try:
                        high = max(high, int(dfct))
                    except ValueError:
                        pass

                self._buffer = str(high + 1)

                if key == "m":
                    self._severity = defect.Severity.MODERATE
                elif key == "n":
                    self._severity = defect.Severity.SEVERE
                self._point()

            elif self.ffield and key == "l":
                # add a line defect
                print(key)
                self._line_stage = 0
                self._line_mode = True

                # Try to calculate the largest number, label the new Defect
                # with one more
                high = 0
                for dfct in self.mccd_dfct[self._cnam]:
                    try:
                        high = max(high, int(dfct))
                    except ValueError:
                        pass

                self._buffer = str(high + 1)

                self._line()

            elif key == "d":
                # delete an defect
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

                # Re-order to have severe after moderate defects
                # create empty output container
                out_dfct = defect.MccdDefect()
                for cnam, ccd_dfct in self.mccd_dfct.items():

                    # initialise output object for this CCD
                    ccd_dfct_out = out_dfct[cnam] = defect.CcdDefect()

                    # first add in moderates
                    n = 0
                    for dnam, dfct in ccd_dfct.items():
                        if dfct.severity == defect.Severity.MODERATE:
                            ccd_dfct_out[str(n + 1)] = dfct
                            n += 1

                    # then the severes
                    for dnam, dfct in ccd_dfct.items():
                        if dfct.severity == defect.Severity.SEVERE:
                            ccd_dfct_out[str(n + 1)] = dfct
                            n += 1

                # old files are over-written at this point
                out_dfct.write(self.dfctnam)
                print("\nDefects saved to {:s}.\nBye".format(self.dfctnam))

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

            elif key == "l" and not self.ffield:
                print("\nCannot add lines of hot pixels, only flat-field defects")

            else:
                print('\nNo action is defined for key = "{:s}"'.format(key))
                self.action_prompt(False)

    def _point(self):
        """Once all set to add a Point defect, this routine actually carries out the
        necessary operations

        """

        # search for enclosing window, print stats
        wnam, wind = utils.print_stats(
            self.mccd[self._cnam], self._cnam, self._x, self._y, self.hsbox, False
        )

        if wnam is None:
            print("  cannot set defects outside windows")

        else:

            if self.ffield:
                # flat-field defect
                dfct = defect.Point(self._severity, self._x, self._y)

            else:
                # hot pixel
                dfct = defect.Hot(self._severity, self._x, self._y)

            self.mccd_dfct[self._cnam][self._buffer] = dfct

            # add defect to the plot, store plot objects
            self.pobjs[self._cnam][self._buffer] = hcam.mpl.pDefect(self._axes, dfct)

            # make sure it appears
            plt.draw()

            # let user know what has happened
            level = (
                "moderate" if self._severity == defect.Severity.MODERATE else "severe"
            )

            print(
                "added {:s} level defect {:s} to CCD {:s} at x,y = {:.2f},{:.2f}".format(
                    level, self._buffer, self._cnam, self._x, self._y
                )
            )

        self.action_prompt(True)

    def _line(self):
        """Once all set to add a Line defect, this routine actually carries
        out the necessary operations

        """

        self._line_stage += 1

        if self._line_stage == 1:

            wnam = self.mccd[self._cnam].inside(self._x, self._y, 0)
            if wnam is None:
                self._line_mode = False
                print("  cannot set defects outside windows")
                self.action_prompt(True)

            else:

                # store the CCD, and the first x,y position
                self._line_cnam = self._cnam
                self._line_x1 = self._x
                self._line_y1 = self._y

                # prompt stage 2
                print(" second point: s(elect) or q(uit)")

        elif self._line_stage == 2:

            wnam = self.mccd[self._cnam].inside(self._x, self._y, 0)
            if wnam is None:
                self._line_mode = False
                print("  cannot set defects outside windows")
                self.action_prompt(True)

            else:

                # now the second x,y position
                if self._cnam != self._line_cnam:
                    self._line_mode = False
                    print("  cannot set defects across different CCDs")
                    self.action_prompt(True)

                self._line_x2 = self._x
                self._line_y2 = self._y

                # prompt stage 2
                print(" Defect level: m(odest), n(asty) [else q(uit)]")

        elif self._line_stage == 3:

            self._line_mode = False

            dfct = defect.Line(
                self._severity,
                self._line_x1,
                self._line_y1,
                self._line_x2,
                self._line_y2,
            )
            self.mccd_dfct[self._cnam][self._buffer] = dfct

            # add defect to the plot, store plot objects
            self.pobjs[self._cnam][self._buffer] = \
                hcam.mpl.pDefect(self._axes, dfct)

            # make sure it appears
            plt.draw()

            # let user know what has happened
            level = (
                "moderate" if \
                self._severity == defect.Severity.MODERATE else "severe"
            )

            print(
                (
                    "added {:s} level line defect {:s}"
                    " to CCD {:s} from x1,y1 = {:.2f},{:.2f}"
                    " to x2,y2 = {:.2f},{:.2f}"
                ).format(
                    level,
                    self._buffer,
                    self._cnam,
                    self._line_x1,
                    self._line_y1,
                    self._line_x2,
                    self._line_y2,
                )
            )
            self.action_prompt(True)

    def _show(self):
        """
        Prints stats on pixels around selected place
        """

        # search for enclosing window, print stats
        wnam, wind = utils.print_stats(
            self.mccd[self._cnam], self._cnam, self._x, self._y, self.hsbox, False
        )
        if wnam is None:
            print('  must hit "s" inside a window')

        self.action_prompt(True)

    def _delete(self):
        """This deletes the nearest defect to the currently selected
        position, if it is near enough (within 10 pixels)

        """

        # first see if there is a Defect near enough the selected position
        dfct, dfnam, dmin = self._find_defect()

        if dmin is not None and dmin < 10:

            # near enough for deletion
            self.pobjs[self._cnam][dfnam].remove()

            # delete Defect from containers
            del self.pobjs[self._cnam][dfnam]
            del self.mccd_dfct[self._cnam][dfnam]

            # update plot
            plt.draw()
            print('  deleted defect "{:s}"'.format(dfnam))

        else:
            print("  no defect near enough the cursor position for deletion")

        self.action_prompt(True)

    def _find_defect(self):
        """Finds the nearest Defect to the currently selected position,

        It returns (dfct, dfctnam, dmin) where dfct is the Defect, dfctnam its
        label, and dmin is the minimum distance. These are all returned as
        None if no suitable Defect is found.

        """

        dmin = None
        dfctmin = None
        dnmin = None
        for dfctnam, dfct in self.mccd_dfct[self._cnam].items():
            dist = dfct.dist(self._x, self._y)
            if dmin is None or dist < dmin:
                dmin = dist
                dfctmin = dfct
                dnmin = dfctnam

        return (dfctmin, dnmin, dmin)
