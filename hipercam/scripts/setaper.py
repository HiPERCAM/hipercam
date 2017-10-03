import sys
import os

import numpy as np
import matplotlib.pyplot as plt

from trm.pgplot import *
import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

#############################################
#
# setaper -- defines apertures given an image
#
#############################################

def setaper(args=None):
    """Interactive definition of photometric extraction apertures. This is
    a matplotlib-based routine allowing you to place apertures on targets
    using the cursor.

    Arguments::

      mccd   : (string)
         name of an MCCD file, as produced by e.g. 'grab'

      aper   : (string)
         the name of an aperture file. If it exists it will be read so that apertures
         can be added to it. If it does not exist, it will be created on exiting the
         routine. The aperture files are is a fairly readable / editiable text format

      ccd    : (string)
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      schar  : (bool) [hidden]
         flag determining whether single character input without needing to
         hit <CR> is used for aperture labelling or not. Single character
         input is probably what is wanted when observing (allowed characters
         0-9, a-z, A-Z, no spaces / punctuation) as it's fastest.

      width  : (float) [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : (float) [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

      rtarg  : (float) [unbinned pixels]
         radius of target aperture. The exact value of this does not matter
         too much since it is normally overridden in 'reduce', but typically
         one aims for 1.5 to 2.5 x FWHM, seeing, depending upon the target
         brightness.

      rsky1  : (float) [unbinned pixels]
         inner radius of sky aperture.

      rsky2  : (float) [unbinned pixels]
         radius of target aperture

      nx     : (int)
         number of panels across to display, prompted if more than one CCD is
         to be plotted.

      msub   : (bool)
         True/False to subtract median from each window before scaling

      iset   : (string) [single character]
         determines how the intensities are determined. There are three
         options: 'a' for automatic simply scales from the minimum to the
         maximum value found on a per CCD basis. 'd' for direct just takes two
         numbers from the user. 'p' for percentile dtermines levels based upon
         percentiles determined from the entire CCD on a per CCD bais.

      ilo    : (float) [if iset=='d']
         lower intensity level

      ihi    : (float) [if iset=='d']
         upper intensity level

      plo    : (float) [if iset=='p']
         lower percentile level

      phi    : (float) [if iset=='p']
         upper percentile level

      profit : (bool) [hidden]
         start aperture selection with profile fitting-based position refinement (or not). Can be toggled
         at any point.

      method : (string) [hidden]
         this defines the profile fitting method, if profile fitting is used to refine
         the aperture position. Either a gaussian or a moffat profile, 'g' or 'm'.
         The latter is usually best.

      beta   : (float) [if method == 'm'; hidden]
         default Moffat exponent

      fwhm   : (float) [hidden]
         default FWHM, unbinned pixels.

      fwhm_min : (float) [hidden]
         minimum FWHM to allow, unbinned pixels.

      shbox  : (float) [hidden]
         half width of box for searching for a star, unbinned pixels. The
         brightest target in a region +/- shbox around an intial position will
         be found. 'shbox' should be large enough to allow for likely changes
         in position from frame to frame, but try to keep it as small as you
         can to avoid jumping to different targets and to reduce the chances
         of interference by cosmic rays.

      smooth : (float) [hidden]
         FWHM for gaussian smoothing, binned pixels. The initial position for
         fitting is determined by finding the maximum flux in a smoothed
         version of the image in a box of width +/- shbox around the starter
         position. Typically should be comparable to the stellar width. Its
         main purpose is to combat cosmi rays which tend only to occupy a
         single pixel.

      splot  : (bool) [hidden]
         Controls whether an outline of the search box and a target number
         is plotted (in red) or not.

      fhbox  : (float) [hidden]
         half width of box for profile fit, unbinned pixels. The fit box is centred
         on the position located by the initial search. It should normally be
         > ~2x the expected FWHM.

      read   : (float) [hidden]
         readout noise, RMS ADU, for assigning uncertainties

      gain   : (float) [hidden]
         gain, ADU/count, for assigning uncertainties

      sigma  : (float) [hidden]
         sigma rejection threshold


    Notes. There are a few conveniences to make setaper easier::

      1. You can toggle between the pan/zoom and input modes using the 'Alt'
         button as well as the usual icon in the toolbar.

      2. Clicking any mouse button acts like hitting 'a' for adding an aperture.

      3. The label input can be switched between single- and multi-character
         input ('schar'). Single character input requires no <enter>.

    Various standard keyboard shortcuts (e.g. 's' to save) are disabled as
    they just confuse things and are of limited use in setaper in any case.

    Some aspects of the usage of matplotlib in setaper are tricky. It is possible
    that particular 'backends' will cause problems. I have tested this with Qt4Agg.
    """

    if args is None:
        args = sys.argv[1:]

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', 'setaper', args) as cl:

        # register parameters
        cl.register('mccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('aper', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('schar', Cline.LOCAL, Cline.HIDE)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)
        cl.register('rtarg', Cline.LOCAL, Cline.PROMPT)
        cl.register('rsky1', Cline.LOCAL, Cline.PROMPT)
        cl.register('rsky2', Cline.LOCAL, Cline.PROMPT)
        cl.register('nx', Cline.LOCAL, Cline.PROMPT)
        cl.register('msub', Cline.GLOBAL, Cline.PROMPT)
        cl.register('iset', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ilo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ihi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('plo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('phi', Cline.GLOBAL, Cline.PROMPT)
#       cl.register('fwidth', Cline.LOCAL, Cline.HIDE)
#       cl.register('fheight', Cline.LOCAL, Cline.HIDE)
        cl.register('profit', Cline.LOCAL, Cline.HIDE)
        cl.register('method', Cline.LOCAL, Cline.HIDE)
        cl.register('beta', Cline.LOCAL, Cline.HIDE)
        cl.register('fwhm', Cline.LOCAL, Cline.HIDE)
        cl.register('fwhm_min', Cline.LOCAL, Cline.HIDE)
        cl.register('shbox', Cline.LOCAL, Cline.HIDE)
        cl.register('smooth', Cline.LOCAL, Cline.HIDE)
#        cl.register('splot', Cline.LOCAL, Cline.HIDE)
        cl.register('fhbox', Cline.LOCAL, Cline.HIDE)
        cl.register('read', Cline.LOCAL, Cline.HIDE)
        cl.register('gain', Cline.LOCAL, Cline.HIDE)
        cl.register('sigma', Cline.LOCAL, Cline.HIDE)

        # get inputs
        mccd = cl.get_value('mccd', 'frame to plot',
                            cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(mccd)

        aper = cl.get_value('aper', 'frame to plot',
                            cline.Fname('hcam', hcam.APER, exist=False))

        if os.path.exists(aper):
            # read in old apertures
            with open(aper) as fp:
                mccdaper = hcam.MccdAper.fromJson(fp)
            print('Loaded existing file = {:s}'.format(aper))
        else:
            # create empty container
            mccdaper = hcam.MccdAper()
            print('No file called {:s} exists; will create from scratch'.format(aper))

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

        # next three are usually hidden
        schar = cl.get_value('schar', 'single character input for aperture names?', True)
        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

        # aperture radii
        rtarg = cl.get_value('rtarg', 'target aperture radius [unbinned pixels]', 10., 0.)
        rsky1 = cl.get_value('rsky1', 'inner sky aperture radius [unbinned pixels]', 15., 0.)
        rsky2 = cl.get_value('rsky2', 'outer sky aperture radius [unbinned pixels]', 25., 0.)

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
            'iset', 'set intensity a(utomatically), d(irectly) or with p(ercentiles)?',
            'a', lvals=['a','A','d','D','p','P'])
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == 'd':
            ilo = cl.get_value('ilo', 'lower intensity limit', 0.)
            ihi = cl.get_value('ihi', 'upper intensity limit', 1000.)
        elif iset == 'p':
            plo = cl.get_value('plo', 'lower intensity limit percentile', 5., 0., 100.)
            phi = cl.get_value('phi', 'upper intensity limit percentile', 95., 0., 100.)

        nxmax, nymax = 0, 0
        for cnam in ccds:
            nxmax = max(nxmax, mccd[cnam].nxtot)
            nymax = max(nymax, mccd[cnam].nytot)

        # might be worth trying to improve this at some point
        xlo, xhi, ylo, yhi = 0, nxmax+1, 0, nymax+1

#        fwidth = cl.get_value('fwidth', 'fit plot width (inches)', 0.)
#        fheight = cl.get_value('fheight', 'fit plot height (inches)', 0.)
        profit = cl.get_value('profit', 'use profile fits to refine the aperture positions?', True)
        method = cl.get_value('method', 'fit method g(aussian) or m(offat)', 'm', lvals=['g','m'])
        if method == 'm':
            beta = cl.get_value('beta', 'initial exponent for Moffat fits', 5., 0.5)
        else:
            beta = 0
        fwhm_min = cl.get_value('fwhm_min', 'minimum FWHM to allow [unbinned pixels]', 1.5, 0.01)
        fwhm = cl.get_value('fwhm', 'initial FWHM [unbinned pixels] for profile fits', 6., fwhm_min)
        shbox = cl.get_value(
            'shbox', 'half width of box for initial location of target [unbinned pixels]', 11., 2.)
        smooth = cl.get_value(
            'smooth', 'FWHM for smoothing for initial object detection [binned pixels]', 6.)
#        splot = cl.get_value(
#            'splot', 'plot outline of search box?', True)
        fhbox = cl.get_value(
            'fhbox', 'half width of box for profile fit [unbinned pixels]', 21., 3.)
        read = cl.get_value('read', 'readout noise, RMS ADU', 3.)
        gain = cl.get_value('gain', 'gain, ADU/e-', 1.)
        sigma = cl.get_value('sigma', 'readout noise, RMS ADU', 3.)

    # Inputs obtained.

    # re-configure keyboard shortcuts to avoid otherwise confusing behaviour
    plt.rcParams['keymap.back'] = ''
    plt.rcParams['keymap.forward'] = ''
    plt.rcParams['keymap.fullscreen'] = ''
    plt.rcParams['keymap.grid'] = ''
    plt.rcParams['keymap.home'] = ''
    plt.rcParams['keymap.pan'] = ['alt']
    plt.rcParams['keymap.pan'] = ['alt']
    plt.rcParams['keymap.save'] = ''
    plt.rcParams['keymap.xscale'] = ''
    plt.rcParams['keymap.yscale'] = ''
    plt.rcParams['keymap.zoom'] = ''

    # start plot
    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width,height))
    else:
        fig = plt.figure()

    # get the navigation toolbar which is used to check the pan/zoom mode
    toolbar = fig.canvas.manager.toolbar

    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    # we need to store some stuff
    ax = None
    cnams = {}

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

        if cnam in mccdaper:
            # plot any pre-existing apertures, keeping track of
            # the plot objects
            pobjs[cnam] = hcam.mpl.pCcdAper(axes, mccdaper[cnam])

        else:
            # add in an empty CcdApers for any CCD not already present
            mccdaper[cnam] = hcam.CcdAper()

            # and an empty container for any new plot objects
            pobjs[cnam] = hcam.Group(tuple)

    # create the aperture picker (see below for class def)
    picker = PickStar(mccd, cnams, toolbar, fig, mccdaper, schar,
                      rtarg, rsky1, rsky2, profit, method, beta,
                      fwhm_min, fwhm, shbox, smooth, fhbox,
                      read, gain, sigma, aper, pobjs)

    plt.tight_layout()
    picker.action_prompt(False)

    # finally show stuff ....
    plt.show()


# the next class is where all the action occurs. A rather complicated matter
# of handling events. note that standard terminal input with 'input' becomes
# impossible, explaining some of the weirdness (accumulation mode)

class PickStar:
    """Class to pick targets for apertures.
    """

    ADD_PROMPT = "enter a label for the aperture, '!' to abort: "

    def __init__(self, mccd, cnams, toolbar, fig, mccdaper, schar,
                 rtarg, rsky1, rsky2, profit, method, beta,
                 fwhm_min, fwhm, shbox, smooth, fhbox,
                 read, gain, sigma, apernam, pobjs):
        """axes is dictionary keyed by CCD label of the axes. all the args are
        stored as attributes to allow "global"-type access to a call back
        added to fig.canvas to handle KeyPress events from the user.

          mccd     : (hipercam.MCCD)
            the MCCD to define the apertures

          cnams    : (dict)
            dictionary keyed on Axes instances giving corresponding CCD labels.
            They are in 1-to-1 correspondence.

          toolbar  : (NavigationToolbar2)
             represents the plot toolbar

          fig      : (Figure)
             the Figure object of the whole plot

          mccdaper : (MccdAper)
             contains whatever Apertures have been loaded from file

          schar    : (bool)
             schar == True implies apertures labels are input as single characters
             without needing <CR>.

          rtarg    : (float)
             target aperture radius, unbinned pixels

          rsky1    : (float)
             inner sky aperture radius, unbinned pixels

          rsky2    : (float)
             outer sky aperture radius, unbinned pixels

          profit   : (bool)
             whether to apply profile fitting to refine the aperture positions. This carries
             out a source search before the fit.

          apernam  : (string)
             name of file to store final aperture set

          pobjs    : (Group(Group))
             container for all plot objects for the apertures to allow
             deletion. Specifically it stores tuples of plot objects so the
             inner Groups must be created starting Group(tuple,...)

        """

        # save the inputs, tack on event handlers.
        self.mccd = mccd
        self.cnams = cnams
        self.toolbar = toolbar
        self.fig = fig
        self.fig.canvas.mpl_connect('key_press_event',
                                    self._keyPressEvent)
        self.fig.canvas.mpl_connect('button_press_event',
                                    self._buttonPressEvent)
        self.mccdaper = mccdaper
        self.schar = schar
        self.rtarg = rtarg
        self.rsky1 = rsky1
        self.rsky2 = rsky2
        self.rsky2 = rsky2
        self.profit = profit
        self.method = method
        self.beta = beta
        self.fwhm_min = fwhm_min
        self.fwhm = fwhm
        self.shbox = shbox
        self.smooth = smooth
        self.fhbox = fhbox
        self.read = read
        self.gain = gain
        self.sigma = sigma
        self.apernam = apernam
        self.pobjs = pobjs

        # the next are to do with entering / leaving accumulation mode and
        # storing the state at the point when we enter accumulation mode
        self._accum = False
        self._action = None
        self._buffer = None
        self._cnam = None
        self._key = None
        self._x = None
        self._y = None
        self._axes = None

    def action_prompt(self, cr):
        """
        Prompts user for an action. cr controls whether there is an intial carriage return or not.
        It leaves the cursor at the end of the line. This appears in multiple places hence why it
        is a method.
        """
        if cr:
            print()

        print('a(dd), c(entre) d(elete), h(elp), p(rofit), r(eference), q(uit): ', end='',flush=True)


    def _keyPressEvent(self, event):
        """
        This is where we do the hard work. Every key press event is diverted
        to this method. It either takes an action based on the input, such as
        removing an Aperture, or sometimes it causes a state change such that
        input is diverted to and accumulated in a buffer until 'enter' is hit.
        The latter stage comes first.
        """

        if self._accum:
            # accumulation mode action
            self._accumulate(event.key)

        else:
            # standard mode action
            self._standard(event.key, event.xdata, event.ydata, event.inaxes)


    def _buttonPressEvent(self, event):
        """
        Treat a mouse button press event as equivalent to 'enter' as a convenience, if we are
        not in pan/zoom mode that is
        """
        if self.toolbar.mode != 'pan/zoom':
            if self._accum:
                # accumulation mode action, pretend we have hit 'enter'
                self._accumulate('enter')

            elif event.inaxes is not None:
                # standard mode action, either pretend we have hit 'a' for add
                # or 'd' for delete if the cursor is inside an existing aperture

                for anam, aper in self.mccdaper[self._cnam].items():
                    if np.sqrt((aper.x-event.xdata)**2+(aper.y-event.ydata)**2) < self.rtarg:
                        self._standard('d', event.xdata, event.ydata, event.inaxes)
                        break
                else:
                    self._standard('a', event.xdata, event.ydata, event.inaxes)

    def _accumulate(self, key):
        """Carries out the work needed when we are in accumulation mode. Just pass
        through the value of the key pressed and this will handle the rest. It
        defines how we treat characters.

        """

        if key == 'enter':
            # trap 'enter'
            print()
            if self._action == 'a':
                if self._buffer in self.mccdaper[self._cnam]:
                    print(
                        'label={:s} already in use; please try again'.format(self._buffer),
                        file=sys.stderr
                        )
                    print(PickStar.ADD_PROMPT, end='',flush=True)
                    self._buffer = ''

                elif self._buffer == '':
                    print(
                        'label blank; please try again'.format(self._buffer),
                        file=sys.stderr
                        )
                    print(PickStar.ADD_PROMPT, end='',flush=True)

                else:
                    # add & plot aperture
                    self._add_aperture()
                    self._accum = False

        elif key == '!' and self._buffer == '':
            # terminate accumulation mode without bothering to wait for an 'enter'
            print('\n*** no aperture added')
            self.action_prompt(True)
            self._accum = False

        elif key == 'backspace' or key == 'delete':
            # remove a character 
            self._buffer = self._buffer[:-1]
            print('{:s}{:s} '.format(PickStar.ADD_PROMPT, self._buffer))

        elif key in '!0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ':
            # accumulate input and add to the
            self._buffer += key
            if self.schar:
                # single character input. bail out immediately (inside _add_aperture)
                print(key)
                if self._action == 'a':
                    if self._buffer in self.mccdaper[self._cnam]:
                        print(
                            'label={:s} already in use; please try again'.format(self._buffer),
                            file=sys.stderr
                            )
                        print(PickStar.ADD_PROMPT, end='',flush=True)
                        self._buffer = ''

                    elif self._buffer == '':
                        print(
                            'label blank; please try again'.format(self._buffer),
                            file=sys.stderr
                            )
                        print(PickStar.ADD_PROMPT, end='',flush=True)

                    else:
                        # add & plot aperture
                        self._add_aperture()
                        self._accum = False

            else:
                # multi character input. just accumulate characters
                print(key, end='', flush=True)

    def _standard(self, key, x, y, axes):
        """
        Carries out the work needed when we are in the standard mode. Just pass through
        the value of the key pressed, the associated x, y position and the axes instance
        (all from the event) and this will handle the rest.
        """

        if key == 'alt':
            # toggle pan/zoom mode in which the plots can be zoomed and panned
            # versus the mode which accepts character input to add apertures etc
            if self.toolbar.mode == 'pan/zoom':
                print('\nswitched to pan/zoom mode')
            else:
                print('\nswitched to aperture update mode')
                self.action_prompt(False)

        elif self.toolbar.mode != 'pan/zoom' and axes is not None:

            # store information in attributes accessible to all methods, sort
            # of restricted global variables in effect.
            self._cnam = self.cnams[axes]
            self._key = key
            self._axes = axes
            self._axes = axes
            self._x = x
            self._y = y

            if key == 'h':
                print(key)
                print("""

Help on the actions available in 'setaper'. Actions:

 a(dd)      : add an aperture
 c(entre)   : centres an aperture by fitting nearby star
 d(elete)   : delete an aperture
 h(elp)     : print this help text
 p(rofit)   : toggles between fitting+position correction and no fits
 r(eference): marks an aperture as a reference aperture
 q(uit)     : quit setaper and save apertures

A few shortcuts exist. The 'alt' key can be used to toggle between the
aperture input mode where one starts, and 'pan & zoom' mode in which the plots
can be zoomed in and out and moved (same action as clicking the crossed arrows
in the top tool bar). Clicking the mouse buttons is equivalent to hitting 'a'
to add an aperture unless you are inside the circle of the target aperture in
which case it is taken as a 'd' to delete it. Thus clicking in an empty area
then entering a character will label a new Aperture and another click without
having moved the cursor will remove it. Clicking the mouse with one hand and
entering labels with the other can be the fastest way to add a series of
apertures. During multi-character entry, clicking the mouse buttons is
equivalent to hitting the 'enter' key to terminate input which can free up the
need to move your hand from mouse to keyboard.

Hitting 'd' will delete the aperture nearest to the cursor, as long as it is
'close' defined as being within max(rtarg,min(100,max(20,2*rsky2))) pixels of
the aperture centre.
""")

            elif key == 'a':
                print(key)
                print(PickStar.ADD_PROMPT, end='',flush=True)
                self._accum = True
                self._buffer = ''
                self._action = 'a'

            elif key == 'c':
                print(key)
                self._centre_aperture()

            elif key == 'd':
                print(key)
                self._delete_aperture()

            elif key == 'p':
                print(key)
                self.profit = not self.profit
                if self.profit:
                    print(' turned on profile fitting & re-positioning when adding apertures')
                else:
                    print(' switched off profile fitting & re-positioning when adding apertures')
                self.action_prompt(True)

            elif key == 'q':
                print(key)
                # quit and clear up
                plt.close()

                # old files are over-written at this point
                with open(self.apernam,'w') as fp:
                    self.mccdaper.toJson(fp)
                print(' apertures saved to {:s}.\nBye'.format(self.apernam))

            elif key == 'r':
                print(key)
                self._toggle_reference()

            elif key == 'enter':
                self.action_prompt(True)

            elif key == 'shift':
                # trap some special keys to avoid irritating messages
                pass

            else:
                print('\nNo action is defined for key = "{:s}"'.format(key))
                self.action_prompt(False)

    def _add_aperture(self):
        """
        Once all set to add an aperture, and relevant attribute values are all
        set, this routine actually carries out the necessary operations which
        are (i) refine the position of the aperture, (ii) create the aperture,
        (iii) store in the multiaperture object, (iv) plot it.
        """

        if self.profit:

            print ('  fitting ...')

            # extract the CCD
            ccd = self.mccd[self._cnam]

            # check that the selected position is inside a window
            wnam = ccd.inside(self._x, self._y, 2)
            if wnam is None:
                print('  *** selected position ({:.1f},{:.1f}) not in a window; should not occur'.format(self._x,self._y), file=sys.stderr)
                self.action_prompt(True)
                return

            # get Windat around the selected position
            wind = ccd[wnam].window(self._x-self.shbox, self._x+self.shbox, self._y-self.shbox, self._y+self.shbox)

            # carry out initial search
            x,y,peak = wind.find(self.smooth, False)

            # now for a more refined fit. First extract fit Windat
            fwind = ccd[wnam].window(x-self.fhbox, x+self.fhbox, y-self.fhbox, y+self.fhbox)
            sky = np.percentile(fwind.data, 25)

            # refine the Aperture position by fitting the profile
            if self.method == 'g':
                # gaussian fit
                (sky, peak, x, y, fwhm), sigs, (fit, X, Y, weights) = \
                    fwind.fitGaussian(sky, peak-sky, x, y, self.fwhm, self.fwhm_min, self.read, self.gain, self.sigma)

                if sigs is None:
                    print('  *** Aperture {:s}: fit failed'.format(self._buffer), file=sys.stderr)
                    self.action_prompt(True)
                    return
                else:
                    esky, epeak, ex, ey, efwhm = sigs 
                    print('  x,y = {:.1f}({:.1f}),{:.1f}({:.1f}), FWHM = {:.2f}({:.2f}), peak = {:.1f}({:.1f}), sky = {:.1f}({:.1f})'.format(
                            self._buffer,x,ex,y,ey,fwhm,efwhm,peak,epeak,sky,esky)
                          )
                    print(  'shifted by dx,dy = {:.2f},{.2f}'.format(x-self._x,y-self._y)) 
                    self._x = x
                    self._y = y

            else:
                raise NotImplementedError('{:s} fitting method not implemented'.format(method))

        # create and add aperture
        aper = hcam.Aperture(self._x, self._y, self.rtarg, self.rsky1, self.rsky2, False)
        self.mccdaper[self._cnam][self._buffer] = aper

        # add aperture to the plot, store plot objects
        self.pobjs[self._cnam][self._buffer] = hcam.mpl.pAper(self._axes, aper, self._buffer)

        # make sure it appears
        plt.draw()

        print('added aperture {:s} to CCD {:s} at x,y = {:.2f},{:.2f}'.format(
                self._buffer,self._cnam,self._x,self._y)
              )

        self.action_prompt(True)


    def _delete_aperture(self):
        """
        This deletes the nearest aperture to the currently selected
        position, if it is near enough. It will also break any links
        from other Apertures to the Aperture that is deleted.

        'Near enough' is defined as within max(rtarg,min(100,max(20,2*rsky2))) of
        the aperture centre.
        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is not None and dmin < max(self.rtarg,min(100,max(20.,2*self.rsky2))):
            # near enough for deletion
            for obj in self.pobjs[self._cnam][apnam]:
                obj.remove()

            # update plot
            plt.draw()
            print('  deleted aperture "{:s}"'.format(apnam))

            # delete Aperture from containers
            del self.pobjs[self._cnam][apnam]
            del self.mccdaper[self._cnam][apnam]

            # break any links to the deleted aperture
            for anam, aper in self.mccdaper[self._cnam].items():
                if aper.link == apnam:
                    aper.link = ''
                    print('  removed link to aperture "{:s}" from aperture "{:s}"'.format(
                            apnam, anam))

        else:
            print('  found no aperture near enough the cursor position for deletion')

        self.action_prompt(True)

    def _centre_aperture(self):
        """
        Centres an aperture using the current fit parameters.
        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is None or dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):
            print('  *** found no aperture near enough the cursor position to re-centre')
        else:
            # OK, there is one near enough
            print ('  fitting ...')

            # extract the CCD
            ccd = self.mccd[self._cnam]

            # check that the selected position is inside a window
            wnam = ccd.inside(self._x, self._y, 2)
            if wnam is None:
                print('  *** selected position ({:.1f},{:.1f}) not in a window; should not occur'.format(self._x,self._y), file=sys.stderr)
            else:
                # get Windat around the selected position
                wind = ccd[wnam].window(self._x-self.shbox, self._x+self.shbox, self._y-self.shbox, self._y+self.shbox)

                # carry out initial search
                x,y,peak = wind.find(self.smooth, False)

                # now for a more refined fit. First extract fit Windat
                fwind = ccd[wnam].window(x-self.fhbox, x+self.fhbox, y-self.fhbox, y+self.fhbox)
                sky = np.percentile(fwind.data, 25)

                # refine the Aperture position by fitting the profile
                if self.method == 'g':
                    # gaussian fit
                    (sky, peak, x, y, fwhm), sigs, (fit, X, Y, weights) = \
                        fwind.fitGaussian(sky, peak-sky, x, y, self.fwhm, self.fwhm_min,
                                          self.read, self.gain, self.sigma)

                    if sigs is None:
                        print('  *** Aperture {:s}: fit failed; current aperture left as is'.format(anmin), file=sys.stderr)
                    else:
                        esky, epeak, ex, ey, efwhm = sigs 
                        print('  x,y = {:.1f}({:.1f}),{:.1f}({:.1f}), FWHM = {:.2f}({:.2f}), peak = {:.1f}({:.1f}), sky = {:.1f}({:.1f})'.format(
                                x,ex,y,ey,fwhm,efwhm,peak,epeak,sky,esky)
                              )
                        print('  shifted by dx,dy = {:.2f},{:.2f}'.format(x-self._x,y-self._y)) 

                        # update aperture position
                        aper.x = x
                        aper.y = y

                        # remove old aperture from from plot
                        for obj in self.pobjs[self._cnam][apnam]:
                            obj.remove()

                        # plot in new position, over-writing the plot objects
                        self.pobjs[self._cnam][apnam] = hcam.mpl.pAper(self._axes, aper, apnam)
                        plt.draw()

                else:
                    raise NotImplementedError('{:s} fitting method not implemented'.format(method))

        self.action_prompt(True)

    def _toggle_reference(self):
       """
        Toggles the reference status of an aperture
        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is None or dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):
            print('  *** found no aperture near enough the cursor position')
        else:
            aper.ref = not aper.ref
            if aper.ref:
                print('  aperture {:s} is now a reference aperture'.format(apnam))
            else:
                print('  aperture {:s} is no longer a reference aperture'.format(apnam))

            # remove aperture from plot
            for obj in self.pobjs[self._cnam][apnam]:
                obj.remove()

            # re-plot new version, over-writing plot objects
            self.pobjs[self._cnam][apnam] = hcam.mpl.pAper(self._axes, aper, apnam)

        self.action_prompt(True)


    def _find_aper(self):
        """
        Finds the nearest aperture to the currently selected position,

        It returns (aper, apnam, dmin) where aper is the Aperture, apnam its label,
        and dmin is the minimum distance. These are all returned as None if no
        suitable Aperture is found.
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
