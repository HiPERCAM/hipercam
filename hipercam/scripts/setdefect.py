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
from hipercam import cline, utils, defect
from hipercam.cline import Cline

__all__ = ['setdefect',]

#############################################
#
# setdefect -- defines CCD defects given an image
#
#############################################

def setdefect(args=None):
    """``setdefect mccd defect ccd [linput width height] rtarg rsky1 rsky2 nx
    msub iset (ilo ihi | plo phi) [profit method beta fwmin fwhm fwfix
    shbox smooth splot fhbox read gain thresh]``

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

      ccd    : string
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      width  : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

      nx     : int
         number of panels across to display, prompted if more than one CCD is
         to be plotted.

      msub   : bool
         True/False to subtract median from each window before scaling

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
    they just confuse things and are of limited use in setaper in any case.

    Some aspects of the usage of matplotlib in setaper are tricky. It is
    possible that particular 'backends' will cause problems. I have tested
    this with Qt4Agg, Qt5agg and GTK3Agg. One aspect is the cursor icon in pan
    mode is a rather indistinct hand where one can't tell what is being
    pointed at. I have therefore suppressed this, but only for the tested
    backends. Others would need require further investigation.

    When in `setaper`, help is always available by hitting 'h'. Most of the options
    are self-evident. A few which may not be are:

      | e(xtra)    : add extra pixels to the target aperture. This allows you to sculpt
      |              the shape of the extraction aperure to perhaps include the flux of
      |              a blended star.
      |
      | l(ink)     : link one aperture to another in the same CCD for re-positioning. If
      |              a target is very faint or will disappear, this will make sure that
      |              is position is defined relative to another, brighter star. Best to
      |              choose one close by because of potential refractive distortion.
      |
      | r(eference): toggle whether an aperture is a reference aperture. The re-positioning
      |              can work in two steps: first position reference stars, then the other stars,
      |              using the reference stars to provide a first cut at the position shift of
      |              the rest.
      |
      | C(opy)     : copy apertures of the CCD the cursor is in to all others. This basically
      |              clones apertures across the CCDs. You will need to re-centre each on afterwards.


    .. Warning::

       The aperture positionsimmediately after a copy reflect the origin CCD, and may be somewhat off
       if there are significant offsets between CCDs.


    """

    command, args = utils.script_args(args)

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('mccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('aper', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)
        cl.register('nx', Cline.LOCAL, Cline.PROMPT)
        cl.register('msub', Cline.GLOBAL, Cline.PROMPT)
        cl.register('iset', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ilo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ihi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('plo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('phi', Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        mccd = cl.get_value('mccd', 'frame to plot',
                            cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.read(mccd)

        dfct = cl.get_value('defect', 'name of defect file',
                            cline.Fname('hcam', hcam.DFCT, exist=False))

        if os.path.exists(aper):
            # read in old defects
            mccd_dfct = hcam.MccdDefect.read(dfct)
            print('Loaded existing file = {:s}'.format(dfct))
        else:
            # create empty container
            mccd_dfct = hcam.MccdDefect()
            print(
                'No file called {:s} exists; '
                'will create from scratch'.format(dfct)
            )

        # define the panel grid
        try:
            nxdef = cl.get_default('nx')
        except:
            nxdef = 3

        max_ccd = len(mccd)
        if max_ccd > 1:
            ccd = cl.get_value('ccd', 'CCD(s) to plot [0 for all]', '0')
            if ccd == '0':
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(mccd.keys())

        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

        # number of panels in X
        if len(ccds) > 1:
            nxdef = min(len(ccds), nxdef)
            cl.set_default('nx', nxdef)
            nx = cl.get_value('nx', 'number of panels in X', 3, 1)
        else:
            nx = 1

        # define the display intensities
        msub = cl.get_value('msub', 'subtract median from each window?', True)

        iset = cl.get_value(
            'iset', 'set intensity a(utomatically),'
            ' d(irectly) or with p(ercentiles)?',
            'a', lvals=['a','A','d','D','p','P'])
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == 'd':
            ilo = cl.get_value('ilo', 'lower intensity limit', 0.)
            ihi = cl.get_value('ihi', 'upper intensity limit', 1000.)
        elif iset == 'p':
            plo = cl.get_value(
                'plo', 'lower intensity limit percentile', 5., 0., 100.)
            phi = cl.get_value(
                'phi', 'upper intensity limit percentile', 95., 0., 100.)

        nxmax, nymax = 0, 0
        for cnam in ccds:
            nxmax = max(nxmax, mccd[cnam].nxtot)
            nymax = max(nymax, mccd[cnam].nytot)

        # might be worth trying to improve this at some point
        xlo, xhi, ylo, yhi = 0, nxmax+1, 0, nymax+1

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

    # start plot
    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width,height))
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

    # this is a container for all the objects used to plot apertures to allow
    # deletion. This is Group of Group objects supporting tuple storage. The idea
    # is that pobjs[cnam][anam] returns the objects used to plot aperture anam of
    # CCD cnam. It is initially empty,
    pobjs = hcam.Group(hcam.Group)

    for n, cnam in enumerate(ccds):
        if ax is None:
            axes = ax = fig.add_subplot(ny, nx, n+1)
            axes.set_aspect('equal', adjustable='box')
            axes.set_xlim(xlo,xhi)
            axes.set_ylim(ylo,yhi)
        else:
            axes = fig.add_subplot(ny, nx, n+1, sharex=ax, sharey=ax)
            axes.set_aspect('equal', adjustable='datalim')

        if msub:
            # subtract median from each window
            for wind in mccd[cnam].values():
                wind -= wind.median()

        hcam.mpl.pCcd(
            axes,mccd[cnam],iset,plo,phi,ilo,ihi,'CCD {:s}'.format(cnam)
            )

        # keep track of the CCDs associated with each axes
        cnams[axes] = cnam

        # and axes associated with each CCD
        anams[cnam] = axes

        if cnam in mccd_dfct:
            # plot any pre-existing apertures, keeping track of
            # the plot objects
            pobjs[cnam] = hcam.mpl.pCcdDefect(axes, mccd_dfct[cnam])

        else:
            # add in an empty CcdApers for any CCD not already present
            mccd_dfct[cnam] = hcam.CcdDefect()

            # and an empty container for any new plot objects
            pobjs[cnam] = hcam.Group(tuple)

    # create the aperture picker (see below for class def)
    picker = PickDefect(
        mccd, cnams, anams, toolbar, fig, mccdaper, linput,
        rtarg, rsky1, rsky2, profit, method, beta,
        fwhm, fwhm_min, fwhm_fix, shbox, smooth, fhbox,
        read, gain, thresh, dfct, pobjs
    )

    plt.tight_layout()
    PickStar.action_prompt(False)

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
            self, mccd, cnams, anams, toolbar, fig, mccd_dfct,
            dfctnam, pobjs):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect('key_press_event', self._keyPressEvent)
        self.mccd = mccd
        self.cnams = cnams
        self.anams = anams
        self.toolbar = toolbar
        self.mccd_dfct = mccd_dfct
        self.dfctnam = dfctnam
        self.pobjs = pobjs

        # then mutually exclusive flags to indicate the action we are in
        # for actions that require extra input. We are not in these at the
        # start so we set them False
        self._pixel_mode = False

    @staticmethod
    def action_prompt(cr):
        """Prompts user for an action. cr controls whether there is an intial
        carriage return or not.  It leaves the cursor at the end of the
        line. This appears in multiple places hence why it is a method.

        """
        if cr:
            print()

        print(
            'd(elete), h(elp), l(ine), p(ixel), q(uit): ', end='', flush=True
        )


    def _keyPressEvent(self, event):
        """
        This is where we do the hard work. Every key press event is diverted
        to this method. It either takes an action based on the input, such as
        removing a defect, or sometimes it causes a state change to get other
        input.
        """

        if self._pixel_mode:

            if event.key == 'q':
                self._pixel_mode = False
                print('no pixel defect added')
                PickDefect.action_prompt(True)
                return

            elif event.key == 'm':
                self._severity = defect.Severity.MODERATE

            elif event.key == 's':
                self._severity = defect.Severity.SEVERE

            self._pixel()

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

            if key == 'h':
                # help text
                print(key)
                print("""

Help on the actions available in 'setaper':

  d(elete)   : delete a defect
  h(elp)     : print this help text
  l(ine)     : add a line defect (straight lines only)
  p(ixel)    : add a pixel defect
  q(uit)     : quit 'setdefect' and save the defects to disk

Hitting 'd' will delete the defect nearest to the cursor, as long as it is
close enough.
""")

            elif key == 'p':

                # add a pixel defect
                print(key)

                # Try to calculate the largest number, label the new aperture
                # with one more
                high = 0
                for dfct in self.mccd_dfct[self._cnam]:
                    try:
                        high = max(high, int(aper))
                    except ValueError:
                        pass
                self._buffer = str(high+1)
                self._pixel_stage = 0
                self._pixel_mode = True
                self._pixel()

                print(PickDefect.ADD_PROMPT, end='',flush=True)

            elif key == 'd':
                # delete an defect
                print(key)
                self._delete()

            elif key == 'q':
                print(key)
                # quit and clear up
                plt.close()

                # old files are over-written at this point
                self.mccd_dfct.write(self.dfctnam)
                print('\nDefects saved to {:s}.\nBye'.format(self.dfctnam))

            elif key == 'enter':
                PickDefect.action_prompt(True)

            elif key == 'shift' or key == 'alt' or key == 'control' or \
                 key == 'pagedown' or key == 'pageup':
                # trap some special keys to avoid irritating messages
                pass

            else:
                print('\nNo action is defined for key = "{:s}"'.format(key))
                PickDefect.action_prompt(False)

    def _pixel(self):
        """Once all set to add a Pixel defect, this routine actually carries out the
        necessary operations

        """

        self._pixel_stage += 1

        if self._pixel_stage == 1:

            # store the CCD, the defect label and the position
            self._pixel_cnam = self._cnam
            self._pixel_dnam = self._buffer
            self._pixel_x = self._x
            self._pixel_y = self._y

            # prompt stage 2
            print(" Defect level: m(oderate), s(evere), q(uit)")

        elif self._pixel_stage == 2:

            self._pixel_mode = False

            dfct = defect.Pixel(self._severity, self._pixel_x, self._pixel_y)
            self.mccd_dfct[self._cnam][self._buffer] = dfct

            # add defect to the plot, store plot objects
            self.pobjs[self._cnam][self._buffer] = hcam.mpl.pDefect(
                self._axes, dfct, self._buffer
            )

            # make sure it appears
            plt.draw()

            print('added defect {:s} to CCD {:s} at x,y = {:.2f},{:.2f}'.format(
                self._buffer,self._cnam,self._pixel_x,self._pixel_y)
              )
            PickDefect.action_prompt(True)

    def _delete(self):
        """This deletes the nearest defect to the currently selected
        position, if it is near enough. 

        'Near enough' is defined as within max(rtarg,min(100,max(20,2*rsky2)))
        of the aperture centre.

        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is not None and \
           dmin < max(self.rtarg,min(100,max(20.,2*self.rsky2))):
            # near enough for deletion
            for obj in self.pobjs[self._cnam][apnam]:
                obj.remove()

            # delete Aperture from containers
            del self.pobjs[self._cnam][apnam]
            del self.mccdaper[self._cnam][apnam]

            # break any links to the deleted aperture
            for anam, aper in self.mccdaper[self._cnam].items():
                if aper.link == apnam:
                    aper.link = ''
                    print('  removed link to aperture "{:s}"'
                          ' from aperture "{:s}"'.format(
                              apnam, anam))

                    # remove the plot of the aperture that was linked in to
                    # wipe the link
                    for obj in self.pobjs[self._cnam][anam]:
                        obj.remove()

                    # then re-plot
                    self.pobjs[self._cnam][anam] = hcam.mpl.pAper(
                        self._axes, aper, anam)

            # update plot
            plt.draw()
            print('  deleted aperture "{:s}"'.format(apnam))

        else:
            print('  found no aperture near enough '
                  'the cursor position for deletion')

        PickDefect.action_prompt(True)

    def _mask(self):
        """
        Adds a sky mask to an aperture
        """

        # count the stage we are at
        self._mask_stage += 1

        if self._mask_stage == 1:
            # Stage 1 first see if there is an aperture near enough the
            # selected position
            aper, apnam, dmin = self._find_aper()

            if dmin is None or \
               dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):
                print('  *** found no aperture near to the'
                      ' cursor position to mask; nothing done'
                )
                PickDefect.action_prompt(True)

            else:

                # ok, we have an aperture. store the CCD, aperture label and
                # aperture for future ref.
                self._mask_cnam = self._cnam
                self._mask_aper = aper
                self._mask_apnam = apnam

                # prompt stage 2
                print(" 'm' at the centre of the region to mask ['q' to quit]")

        elif self._mask_stage == 2:

            if self._cnam != self._mask_cnam:
                print('  *** cannot add sky mask across CCDs; no mask added')
                self._mask_mode = False
            else:
                # store mask centre
                self._mask_xcen = self._x
                self._mask_ycen = self._y

                # prompt stage 3
                print(" 'm' at the edge of the region to mask ['q' to quit]")

        elif self._mask_stage == 3:

            # final stage of mask mode
            self._mask_mode = False

            if self._cnam != self._mask_cnam:
                print('  *** cannot add sky mask across CCDs; no mask added')
            else:
                # compute radius
                radius = np.sqrt((self._x-self._mask_xcen)**2 +
                                 (self._y-self._mask_ycen)**2)

                # add mask to the aperture
                self._mask_aper.add_mask(
                    self._mask_xcen-self._mask_aper.x,
                    self._mask_ycen-self._mask_aper.y, radius
                )

                # delete the aperture from the plot
                for obj in self.pobjs[self._cnam][self._mask_apnam]:
                    obj.remove()

                # re-plot new version, over-writing plot objects
                self.pobjs[self._cnam][self._mask_apnam] = hcam.mpl.pAper(
                    self._axes, self._mask_aper, self._mask_apnam,
                    self.mccdaper[self._mask_cnam]
                )
                plt.draw()

                print(
                    '  added mask to aperture {:s} in CCD {:s}'.format(
                            self._mask_apnam, self._mask_cnam)
                )
                PickDefect.action_prompt(True)

    def _find_aper(self):
        """Finds the nearest aperture to the currently selected position,

        It returns (aper, apnam, dmin) where aper is the Aperture, apnam its
        label, and dmin is the minimum distance. These are all returned as
        None if no suitable Aperture is found.

        """

        dmin = None
        apmin = None
        anmin = None
        for anam, aper in self.mccdaper[self._cnam].items():
            dist = np.sqrt((aper.x-self._x)**2+(aper.y-self._y)**2)
            if dmin is None or dist < dmin:
                dmin = dist
                apmin = aper
                anmin = anam

        return (apmin, anmin, dmin)

    def _add_input(self, key):
        """Accumulates input to label an aperture
        """

        if key == 'enter':
            # trap 'enter'
            print()

            if self._buffer in self.mccdaper[self._cnam]:
                print(
                    'label={:s} already in use; please try again'.format(self._buffer),
                    file=sys.stderr
                )
                print(PickDefect.ADD_PROMPT, end='',flush=True)
                self._buffer = ''

            elif self._buffer == '':
                print(
                    'label blank; please try again'.format(self._buffer),
                    file=sys.stderr
                )
                print(PickDefect.ADD_PROMPT, end='',flush=True)

            else:
                # add & plot aperture
                self._add_aperture()
                self._add_mode = False

        elif key == '!' and self._buffer == '':
            # terminate accumulation mode without bothering to wait for an 'enter'
            print('\n*** no aperture added')
            PickDefect.action_prompt(True)
            self._add_mode = False

        elif key == 'backspace' or key == 'delete':
            # remove a character 
            self._buffer = self._buffer[:-1]
            print('{:s}{:s} '.format(PickDefect.ADD_PROMPT, self._buffer))

        elif key in '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ':
            # accumulate input and add to the buffer
            self._buffer += key

            if self.linput == 's':
                # single character input. bail out immediately (inside _add_aperture)
                print(key)

                if key == '0':
                    print("Apertures cannot be labelled just '0'")
                    print(PickDefect.ADD_PROMPT, end='',flush=True)
                    self._buffer = ''

                elif self._buffer in self.mccdaper[self._cnam]:
                    print(
                        'label={:s} already in use; please try again'.format(self._buffer),
                        file=sys.stderr
                    )
                    print(PickDefect.ADD_PROMPT, end='',flush=True)
                    self._buffer = ''

                elif self._buffer == '':
                    print(
                        'label blank; please try again'.format(self._buffer),
                        file=sys.stderr
                    )
                    print(PickDefect.ADD_PROMPT, end='',flush=True)

                else:
                    # add & plot aperture
                    self._add()
                    self._add_mode = False

            else:
                # multi character input. just accumulate characters
                print(key, end='', flush=True)
