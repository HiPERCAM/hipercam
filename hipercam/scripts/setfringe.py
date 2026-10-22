import os

import matplotlib.pyplot as plt
from trm import cline
from trm.cline import Cline

import hipercam as hcam
from hipercam import fringe, utils

__all__ = [
    "setfringe",
]

#############################################
#
# setfringe -- defines pairs of points for fringe measurement
#
#############################################


def setfringe(args=None):
    """``setfringe fmap fringe ccd [width height] nx msub [cmap nhalf] iset (ilo ihi
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

      msub : bool
         True/False to subtract median from each window before scaling

      cmap : str [hidden]
         The colour map to use. "Greys" is the usual; "Greys_r" reverses it.
         There are many others; typing an incorrect one will give a list. "none"
         for matplotlib default.

      nhalf : int [hidden]
         half-width in binned pixels of stats box as offset from central pixel
         nhalf = 1 gives a 3x3 box; hsbox = 2 gives 5x5 etc. This is used by
         the "show" option when setting FringePair and when reporting FringePair
         differences.

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

    command, args = cline.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:
        # register parameters
        cl.register("fmap", Cline.LOCAL, Cline.PROMPT)
        cl.register("fringe", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)
        cl.register("nx", Cline.LOCAL, Cline.PROMPT)
        cl.register("msub", Cline.GLOBAL, Cline.PROMPT)
        cl.register("cmap", Cline.LOCAL, Cline.HIDE)
        cl.register("nhalf", Cline.GLOBAL, Cline.HIDE)
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
            print("No file called {:s} exists; will create from scratch".format(fpair))

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

        cmap = cl.get_value(
            "cmap", "colour map to use ['none' for mpl default]", "Greys"
        )
        cmap = None if cmap == "none" else cmap

        nhalf = cl.get_value("nhalf", "half-width of stats box (binned pixels)", 2, 0)

        iset = cl.get_value(
            "iset",
            "set intensity a(utomatically), d(irectly) or with p(ercentiles)?",
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

    # create the plot
    fig, caxes = hcam.mpl.plot_frame(
        mccd=mccd,
        cnams=ccds,
        ncols=nx,
        msub=msub,
        cmap=cmap,
        iset=iset,
        ilim=(plo, phi) if iset == "p" else (ilo, ihi),
        region=(xlo, xhi, ylo, yhi),
        title='setfringe',
        figsize=(width, height),
    )

    # add any existing fringe pairs to the CCDs, storing a reference to the plot objects
    pobjs = {}
    for cnam, axes in caxes.items():
        if cnam in mccd_fpair:
            pobjs[cnam] = hcam.mpl.pCcdFringePair(axes, mccd_fpair[cnam])
        else:
            # add in an empty group for any CCD not already present
            mccd_fpair[cnam] = fringe.CcdFringePair()
            pobjs[cnam] = {}

    # create the picker class
    picker = PickFringePair(
        mccd, caxes, fig, mccd_fpair, fpair, nhalf, pobjs
    )
    picker.action_prompt(False)

    # finally set up the plot for interactive mode and display
    hcam.mpl.setup_plot(fig)
    plt.show()


class PickFringePair:
    """Class to pick fringe pairs"""

    def __init__(
        self, mccd, caxes, fig, mccd_fpair, fringenam, nhalf, pobjs
    ):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect("key_press_event", self._keyPressEvent)
        self.mccd = mccd
        self.caxes = caxes
        self.axnams = {axes: cnam for cnam, axes in caxes.items()}  # reverse for events
        self.mccd_fpair = mccd_fpair
        self.fringenam = fringenam
        self.nhalf = nhalf
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

            if event.inaxes is not None:
                self._cnam = self.axnams[event.inaxes]
                self._key = event.key
                self._axes = event.inaxes
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
            self._cnam = self.axnams[axes]
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

                self._label = self.mccd_fpair[self._cnam].next_label()
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

            elif wnam != self._first_wnam:
                print("  cannot set pairs across different windows")

            elif self._cnam != self._first_cnam:
                print("  cannot set pairs across different CCDs")

            else:
                # add new pair
                frng = fringe.FringePair(self._first_x, self._first_y, self._x, self._y)
                self.mccd_fpair[self._cnam][self._label] = frng

                # add fringe pair to the plot, store plot objects
                self.pobjs[self._cnam][self._label] = hcam.mpl.pFringePair(
                    self._axes, frng
                )

                # make sure it appears
                plt.draw()

                # let user know what has happened
                diff = frng.diff(self.mccd[self._cnam], self.nhalf)
                print(
                    (
                        f"added fringe pair to CCD {self._cnam} "
                        f"from x1,y1 = {self._first_x:.2f},{self._first_y:.2f} "
                        f"to x2,y2 = {self._x:.2f},{self._y:.2f}. "
                        f"Median difference = {diff:.2f}"
                    )
                )

            self.action_prompt(True)

    def _show(self):
        """
        Prints stats on pixels around selected place
        """

        # search for enclosing window, print stats
        wnam, wind = utils.print_stats(
            self.mccd[self._cnam], self._cnam, self._x, self._y, self.nhalf, False
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
