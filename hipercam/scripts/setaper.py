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
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ = ['setaper',]

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
         the name of an aperture file. If it exists it will be read so that
         apertures can be added to it. If it does not exist, it will be
         created on exiting the routine. The aperture files are is a fairly
         readable / editiable text format

      ccd    : (string)
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      linput  : (string) [hidden]
         sets the way in which the apertures are labelled. 'n' = numerical
         input, with the program just incrementing by 1 for each successive
         aperture; 's' = single character (without requiring the user to hit
         hit <CR>); 'm' = multi-character, ending with <CR>.
         Allowed characters are 0-9, a-z, A-Z, no spaces or punctuation, but a
         single '0' on its own is not permitted.

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
         start aperture selection with profile fitting-based position
         refinement (or not). Can be toggled at any point.

      method : (string) [hidden]
         this defines the profile fitting method, if profile fitting is used
         to refine the aperture position. Either a gaussian or a moffat
         profile, 'g' or 'm'.  The latter should usually be best.

      beta   : (float) [if method == 'm'; hidden]
         default Moffat exponent

      fwmin  : (float) [hidden]
         minimum FWHM to allow, unbinned pixels.

      fwhm   : (float) [hidden]
         default FWHM, unbinned pixels.

      fwfix  : (bool) [hidden]
         don't fit the FWHM. Can be more robust; the position is still fitted.

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
         main purpose is to combat cosmic rays which tend only to occupy a
         single pixel.

      splot  : (bool) [hidden]
         Controls whether an outline of the search box and a target number
         is plotted (in red) or not.

      fhbox  : (float) [hidden]
         half width of box for profile fit, unbinned pixels. The fit box is
         centred on the position located by the initial search. It should
         normally be > ~2x the expected FWHM.

      read   : (float) [hidden]
         readout noise, RMS ADU, for assigning uncertainties

      gain   : (float) [hidden]
         gain, ADU/count, for assigning uncertainties

      sigma  : (float) [hidden]
         sigma rejection threshold


    Notes. There are a few conveniences to make setaper easier::

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

    """

    if args is None:
        args = sys.argv[1:]

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', 'setaper', args) as cl:

        # register parameters
        cl.register('mccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('aper', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('linput', Cline.LOCAL, Cline.HIDE)
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
        cl.register('fwmin', Cline.LOCAL, Cline.HIDE)
        cl.register('fwfix', Cline.LOCAL, Cline.HIDE)
        cl.register('shbox', Cline.LOCAL, Cline.HIDE)
        cl.register('smooth', Cline.LOCAL, Cline.HIDE)
        cl.register('fhbox', Cline.LOCAL, Cline.HIDE)
        cl.register('read', Cline.LOCAL, Cline.HIDE)
        cl.register('gain', Cline.LOCAL, Cline.HIDE)
        cl.register('sigma', Cline.LOCAL, Cline.HIDE)

        # get inputs
        mccd = cl.get_value('mccd', 'frame to plot',
                            cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.read(mccd)

        aper = cl.get_value('aper', 'frame to plot',
                            cline.Fname('hcam', hcam.APER, exist=False))

        if os.path.exists(aper):
            # read in old apertures
            mccdaper = hcam.MccdAper.read(aper)
            print('Loaded existing file = {:s}'.format(aper))
        else:
            # create empty container
            mccdaper = hcam.MccdAper()
            print(
                'No file called {:s} exists; '
                'will create from scratch'.format(aper)
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

        # next three are usually hidden
        linput = cl.get_value(
            'linput', 'n(umerical), s(ingle) or m(ulti)-character label input',
            'n', lvals=('n','s','m')
        )
        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

        # aperture radii
        rtarg = cl.get_value(
            'rtarg', 'target aperture radius [unbinned pixels]', 10., 0.)
        rsky1 = cl.get_value(
            'rsky1', 'inner sky aperture radius [unbinned pixels]', 15., 0.)
        rsky2 = cl.get_value(
            'rsky2', 'outer sky aperture radius [unbinned pixels]', 25., 0.)

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

#        fwidth = cl.get_value('fwidth', 'fit plot width (inches)', 0.)
#        fheight = cl.get_value('fheight', 'fit plot height (inches)', 0.)
        profit = cl.get_value(
            'profit', 'use profile fits to refine'
            ' the aperture positions?', True
        )
        method = cl.get_value(
            'method', 'fit method g(aussian) or m(offat)', 'm', lvals=['g','m']
        )
        if method == 'm':
            beta = cl.get_value(
                'beta', 'initial exponent for Moffat fits', 5., 0.5)
        else:
            beta = 0
        fwhm_min = cl.get_value(
            'fwmin', 'minimum FWHM to allow [unbinned pixels]', 1.5, 0.01)
        fwhm = cl.get_value(
            'fwhm', 'initial FWHM [unbinned pixels] for profile fits',
            6., fwhm_min
        )
        fwhm_fix = cl.get_value('fwfix', 'fix the FWHM at start value?', False)

        shbox = cl.get_value(
            'shbox', 'half width of box for initial'
            ' location of target [unbinned pixels]', 11., 2.
        )
        smooth = cl.get_value(
            'smooth', 'FWHM for smoothing for initial object'
            ' detection [binned pixels]', 6.
        )
#        splot = cl.get_value(
#            'splot', 'plot outline of search box?', True)
        fhbox = cl.get_value(
            'fhbox', 'half width of box for profile fit'
            ' [unbinned pixels]', 21., 3.)
        read = cl.get_value('read', 'readout noise, RMS ADU', 3.)
        gain = cl.get_value('gain', 'gain, ADU/e-', 1.)
        sigma = cl.get_value('sigma', 'readout noise, RMS ADU', 3.)

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
    picker = PickStar(
        mccd, cnams, anams, toolbar, fig, mccdaper, linput,
        rtarg, rsky1, rsky2, profit, method, beta,
        fwhm, fwhm_min, fwhm_fix, shbox, smooth, fhbox,
        read, gain, sigma, aper, pobjs
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

class PickStar:
    """Class to pick targets for apertures.
    """

    ADD_PROMPT = "enter a label for the aperture, '!' to abort: "

    def __init__(
            self, mccd, cnams, anams, toolbar, fig, mccdaper, linput,
            rtarg, rsky1, rsky2, profit, method, beta, fwhm, fwhm_min,
            fwhm_fix, shbox, smooth, fhbox, read, gain, sigma, apernam,
            pobjs):

        # save the inputs, tack on event handlers.
        self.fig = fig
        self.fig.canvas.mpl_connect('key_press_event', self._keyPressEvent)
        self.mccd = mccd
        self.cnams = cnams
        self.anams = anams
        self.toolbar = toolbar
        self.mccdaper = mccdaper
        self.linput = linput
        self.rtarg = rtarg
        self.rsky1 = rsky1
        self.rsky2 = rsky2
        self.rsky2 = rsky2
        self.profit = profit
        self.method = method
        self.beta = beta
        self.fwhm = fwhm
        self.fwhm_min = fwhm_min
        self.fwhm_fix = fwhm_fix
        self.shbox = shbox
        self.smooth = smooth
        self.fhbox = fhbox
        self.read = read
        self.gain = gain
        self.sigma = sigma
        self.apernam = apernam
        self.pobjs = pobjs

        # then mutually exclusive flags to indicate the action we are in
        # for actions that require extra input. We are not in these at the
        # start so we set them False
        self._add_mode = False
        self._link_mode = False
        self._mask_mode = False
        self._extra_mode = False

    @staticmethod
    def action_prompt(cr):
        """Prompts user for an action. cr controls whether there is an intial
        carriage return or not.  It leaves the cursor at the end of the
        line. This appears in multiple places hence why it is a method.

        """
        if cr:
            print()

        print(
            'a(dd), b(reak), c(entre), d(elete), e(xtra) h(elp), '
            ' l(ink), m(ask), p(rofit), q(uit), r(eference), '
            'C(opy): ', end='', flush=True
        )


    def _keyPressEvent(self, event):
        """
        This is where we do the hard work. Every key press event is diverted
        to this method. It either takes an action based on the input, such as
        removing an Aperture, or sometimes it causes a state change such that
        input is diverted to and accumulated in a buffer until 'enter' is hit.
        The latter stage comes first.
        """

        if self._add_mode:
            # accumulate input
            self._add_input(event.key)

        elif self._link_mode:
            # if in link mode, we should be selecting another aperture, the
            # one to link to
            if event.key == 'q':
                # trap 'q' for quit during linking
                self._link_mode = False
                print('no extra aperture added')
                PickStar.action_prompt(True)

            elif event.key == 'l':
                # link mode. this should be the second aperture
                # store essential data
                self._cnam = self.cnams[event.inaxes]
                self._axes = event.inaxes
                self._x = event.xdata
                self._y = event.ydata
                self._link()

        elif self._extra_mode:
            # add an extra aperture
            if event.key == 'q':
                self._extra_mode = False
                PickStar.action_prompt(True)

            elif event.key == 'e':

                # extra mode. store essential data
                self._cnam = self.cnams[event.inaxes]
                self._axes = event.inaxes
                self._x = event.xdata
                self._y = event.ydata
                self._extra()

        elif self._mask_mode:
            if event.key == 'q':
                self._mask_mode = False
                PickStar.action_prompt(True)

            elif event.key == 'm':

                # mask mode. store essential data
                self._cnam = self.cnams[event.inaxes]
                self._axes = event.inaxes
                self._x = event.xdata
                self._y = event.ydata
                self._mask()

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

  a(dd)      : add an aperture
  b(reak)    : break the link on an aperture
  c(entre)   : centre an aperture by fitting nearby star
  d(elete)   : delete an aperture
  e(xtra)    : add extra pixels to the target aperture
  h(elp)     : print this help text
  l(ink)     : link one aperture to another in the same CCD for re-positioning
  m(ask)     : add a mask to an aperture to ignore regions of sky
  p(rofit)   : toggle between fitting+position correction and no fits
  r(eference): toggle whether an aperture is a reference aperture
  C(opy)     : copy apertures of the CCD the cursor is in to all others (overwrites)
  q(uit)     : quit 'setaper' and save the apertures to disk

Hitting 'd' will delete the aperture nearest to the cursor, as long as it is
close enough. Note that the 'extra' option aperture will be scaled to have the
same size as the main target aperture. The 'mask' apertures have a fixed size.
""")

            elif key == 'a':

                # add an aperture
                print(key)

                if self.linput == 'n':

                    # numerical sequence input. Try to calculate
                    # the largest number, label the new aperture
                    # with one more
                    high = 0
                    for aper in self.mccdaper[self._cnam]:
                        try:
                            high = max(high, int(aper))
                        except ValueError:
                            pass
                    self._buffer = str(high+1)
                    self._add()
                    self._add_mode = False

                else:
                    # switch to add mode, initialise buffer for label input
                    self._add_mode = True
                    self._buffer = ''
                    print(PickStar.ADD_PROMPT, end='',flush=True)

            elif key == 'b':
                # break a link
                print(key)
                self._break()

            elif key == 'c':
                # centre an aperture
                print(key)
                self._centre()

            elif key == 'd':
                # delete an aperture
                print(key)
                self._delete()

            elif key == 'e':
                # add extra target pixels to an aperture
                print(key)

                # switch to extra mode
                self._extra_mode = True
                self._extra_stage = 0
                self._extra()

            elif key == 'l':
                # link an aperture
                print(key)
                if len(self.mccdaper[self._cnam]) < 2:
                    print('need at least 2 apertures in a CCD to be able to make links')
                    PickStar.action_prompt(True)
                else:
                    # switch to link mode, set number of apertures picked
                    self._link_stage = 0
                    self._link_mode = True
                    self._link()

            elif key == 'm':
                # add a sky mask to an aperture
                print(key)

                # switch to mask mode
                self._mask_mode = True
                self._mask_stage = 0
                self._mask()

            elif key == 'p':
                print(key)
                self.profit = not self.profit
                if self.profit:
                    print(' turned on profile fitting & re-positioning'
                          ' when adding apertures')
                else:
                    print(' switched off profile fitting & re-positioning'
                          ' when adding apertures')
                PickStar.action_prompt(True)

            elif key == 'q':
                print(key)
                # quit and clear up
                plt.close()

                # old files are over-written at this point
                self.mccdaper.write(self.apernam)
                print('\nApertures saved to {:s}.\nBye'.format(self.apernam))

            elif key == 'r':
                print(key)
                self._reference()

            elif key == 'C':
                print(key)

                # clone the apertures from this CCD
                if len(self.mccdaper[self._cnam]):

                    # only if there are some to clone ...
                    for cnam in self.mccdaper:

                        if cnam != self._cnam and cnam in self.anams:
                            # ... it is not the CCD to be cloned and it
                            # is being displayed.

                            # remove all existing aperture on this CCD from
                            # the plot
                            for apnam in self.pobjs[cnam]:
                                for obj in self.pobjs[cnam][apnam]:
                                    obj.remove()

                            # copy over the apertures of the CCD to be cloned
                            self.mccdaper[cnam] = self.mccdaper[self._cnam].copy()

                            # plot them, storing the plot objects
                            self.pobjs[cnam] = hcam.mpl.pCcdAper(
                                self.anams[cnam], self.mccdaper[cnam]
                            )

                            plt.draw()

                    print('Copied apertures of CCD {:s}'
                          ' to all other displayed CCDs'.format(self._cnam)
                      )

                PickStar.action_prompt(True)

            elif key == 'enter':
                PickStar.action_prompt(True)

            elif key == 'shift' or key == 'alt' or key == 'control' or \
                 key == 'pagedown' or key == 'pageup':
                # trap some special keys to avoid irritating messages
                pass

            else:
                print('\nNo action is defined for key = "{:s}"'.format(key))
                PickStar.action_prompt(False)

    def _add(self):
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
                print(
                    '  *** selected position ({:.1f},{:.1f}) not in a window;'
                    ' should not occur'.format(self._x,self._y), file=sys.stderr
                )
                PickStar.action_prompt(True)
                return

            # get Windat around the selected position
            wind = ccd[wnam].window(
                self._x-self.shbox, self._x+self.shbox,
                self._y-self.shbox, self._y+self.shbox
            )

            # carry out initial search
            x,y,peak = wind.find(self.smooth, False)

            # now for a more refined fit. First extract fit Windat
            fwind = ccd[wnam].window(x-self.fhbox, x+self.fhbox,
                                     y-self.fhbox, y+self.fhbox)
            sky = np.percentile(fwind.data, 25)

            # refine the Aperture position by fitting the profile
            try:
                (sky, height, x, y, fwhm, beta), epars, \
                    (X, Y, message) = hcam.combFit(
                        fwind, self.method, sky, peak-sky,
                        x, y, self.fwhm, self.fwhm_min, self.fwhm_fix,
                        self.beta, self.read, self.gain
                    )

                print('Aperture {:s}: {:s}'.format(self._buffer,message))
                self._x = x
                self._y = y

            except hcam.HipercamError as err:
                print(err, file=sys.stderr)
                # fit failed.
                PickStar.action_prompt(True)
                return

        # create and add aperture
        aper = hcam.Aperture(
            self._x, self._y, self.rtarg, self.rsky1, self.rsky2, False
        )
        self.mccdaper[self._cnam][self._buffer] = aper
        print('added aperture=',aper)

        # add aperture to the plot, store plot objects
        self.pobjs[self._cnam][self._buffer] = hcam.mpl.pAper(
            self._axes, aper, self._buffer
        )

        # make sure it appears
        plt.draw()

        print('added aperture {:s} to CCD {:s} at x,y = {:.2f},{:.2f}'.format(
            self._buffer,self._cnam,self._x,self._y)
        )
        PickStar.action_prompt(True)

    def _break(self):
        """
        Breaks the link on an aperture (if any)
        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is None or dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):
            print('  *** found no aperture near to the cursor position')

        elif aper.link == '':
            print('  *** there is no link on aperture {:s}'.format(apnam))

        else:

            # cancel the link on the aperture
            self.mccdaper[self._cnam][apnam].break_link()

            # delete the aperture
            for obj in self.pobjs[self._cnam][apnam]:
                obj.remove()

            # re-plot new version, over-writing plot objects
            self.pobjs[self._cnam][apnam] = hcam.mpl.pAper(self._axes, aper, apnam)
            plt.draw()

            print('  cancelled link on aperture {:s} in CCD {:s}'.format(
                    apnam, self._cnam))

        PickStar.action_prompt(True)

    def _centre(self):
        """
        Centres an aperture using the current fit parameters.
        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is None or dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):
            print('  *** found no aperture near to the cursor position to re-centre')
        else:
            # OK, there is one near enough
            print ('  fitting ...')

            # extract the CCD
            ccd = self.mccd[self._cnam]

            # check that the selected position is inside a window
            wnam = ccd.inside(self._x, self._y, 2)
            if wnam is None:
                print('  *** selected position ({:.1f},{:.1f}) not in a window;'
                      ' should not occur'.format(self._x,self._y), file=sys.stderr)
            else:
                # get Windat around the selected position
                wind = ccd[wnam].window(
                    self._x-self.shbox, self._x+self.shbox,
                    self._y-self.shbox, self._y+self.shbox)

                # carry out initial search
                x,y,peak = wind.find(self.smooth, False)

                # now for a more refined fit. First extract fit Windat
                fwind = ccd[wnam].window(
                    x-self.fhbox, x+self.fhbox,
                    y-self.fhbox, y+self.fhbox)
                sky = np.percentile(fwind.data, 25)

                # refine the Aperture position by fitting the profile
                try:
                    (sky, height, x, y, fwhm, beta), epars, \
                        (X, Y, message) = hcam.combFit(
                            fwind, self.method, sky, peak-sky,
                            x, y, self.fwhm, self.fwhm_min, self.fwhm_fix,
                            self.beta, self.read, self.gain
                        )

                    print('Aperture {:s}: {:s}'.format(apnam,message))
                    dx = x - aper.x
                    dy = y - aper.y
                    aper.x = x
                    aper.y = y

                    # remove old aperture from from plot
                    for obj in self.pobjs[self._cnam][apnam]:
                        obj.remove()

                    # plot in new position, over-writing the plot objects
                    self.pobjs[self._cnam][apnam] = hcam.mpl.pAper(
                        self._axes, aper, apnam)

                    # look for any apertures linked to this one and adjust
                    # them so the links are still valid
                    for lanam, laper in self.mccdaper[self._cnam].items():
                        if laper.link == apnam:
                            laper.x += dx
                            laper.y += dy

                            # remove the plot of the aperture that was linked
                            # in to wipe the link
                            for obj in self.pobjs[self._cnam][lanam]:
                                obj.remove()

                            # then re-plot
                            self.pobjs[self._cnam][lanam] = hcam.mpl.pAper(
                                self._axes, laper, lanam, self.mccdaper[self._cnam]
                            )

                    # finally update the plot
                    plt.draw()

                except hcam.HipercamError as err:
                    print(err, file=sys.stderr)
                    print('  *** aperture left as is')

        PickStar.action_prompt(True)

    def _delete(self):
        """This deletes the nearest aperture to the currently selected
        position, if it is near enough. It will also break any links
        from other Apertures to the Aperture that is deleted.

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

        PickStar.action_prompt(True)

    def _extra(self):
        """
        Adds extra target aperture
        """

        # count the stage we are at
        self._extra_stage += 1

        if self._extra_stage == 1:
            # Stage 1 first see if there is an aperture near enough the
            # selected position
            aper, apnam, dmin = self._find_aper()

            if dmin is None or \
               dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):
                print('  *** found no aperture near to the'
                      ' cursor position to add extra pixels; nothing done'
                )
                PickStar.action_prompt(True)

            else:

                # ok, we have an aperture. store the CCD, aperture label and
                # aperture for future ref.
                self._extra_cnam = self._cnam
                self._extra_aper = aper
                self._extra_apnam = apnam

                # prompt stage 2
                print(" 'e' at the centre of the place to add an"
                      " extra aperture ['q' to quit]")

        elif self._extra_stage == 2:

            self._extra_mode = False

            if self._cnam != self._extra_cnam:
                print('  *** cannot add extra apertures across'
                      ' CCDs; no extra aperture added')
            else:
                # add extra to the aperture
                self._extra_aper.add_extra(
                    self._x-self._extra_aper.x,
                    self._y-self._extra_aper.y
                )

                # delete the aperture from the plot
                for obj in self.pobjs[self._cnam][self._extra_apnam]:
                    obj.remove()

                # re-plot new version, over-writing plot objects
                self.pobjs[self._cnam][self._extra_apnam] = hcam.mpl.pAper(
                    self._axes, self._extra_aper, self._extra_apnam,
                    self.mccdaper[self._extra_cnam]
                )
                plt.draw()

                print(
                    '  added extra target aperture to aperture {:s} in CCD {:s}'.format(
                            self._extra_apnam, self._extra_cnam)
                )

            PickStar.action_prompt(True)

    def _link(self):
        """Links one aperture to another. The other one is used when
        re-positioning. The actions of this depend upon the link status ._link
        which is used to separate the first from the second aperture.

        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is None or \
           dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):

            print('  *** found no aperture near to the cursor position to link')
            if self._link_stage == 1:
                print(" 'l' to select the aperture to link"
                      " aperture {:s} from ['q' to quit]".format(apnam))
            else:
                print('  *** no link made.')
                PickStar.action_prompt(True)

        else:

            # ok, we have an aperture
            self._link_stage += 1

            if self._link_stage == 1:

                if aper.ref:
                    print('  *** cannot link a reference aperture')
                    self._link_mode = False

                else:
                    # first time through, store the CCD, aperture label and
                    # aperture of the first aperture which is the one that
                    # will contain the link, set the link status to True and
                    # prompt the user. these values are used second time round
                    self._link_cnam = self._cnam
                    self._link_aper = aper
                    self._link_apnam = apnam
                    print(" 'l' to select the aperture to link"
                          " aperture {:s} from ['q' to quit]".format(apnam))

            else:
                # second (or more) time through. Check we are in the same CCD
                # but on a different aperture first change the link status
                self._link_mode = False

                # then see if we can make a valid link.
                if self._cnam != self._link_cnam:
                    print('  *** cannot link across CCDs; no link made')
                elif apnam == self._link_apnam:
                    print('  *** cannot link an aperture to itself; '
                          'no link made')
                elif self._link_aper.is_linked():
                    print('  *** cannot link an aperture to an aperture'
                          ' that is itself linked; no link made')
                else:
                    # add link to the first aperture
                    self._link_aper.set_link(apnam)

                    # delete the first aperture
                    for obj in self.pobjs[self._cnam][self._link_apnam]:
                        obj.remove()

                    # re-plot new version, over-writing plot objects
                    self.pobjs[self._cnam][self._link_apnam] = hcam.mpl.pAper(
                        self._axes, self._link_aper, self._link_apnam,
                        self.mccdaper[self._link_cnam])
                    plt.draw()

                    print('  linked aperture {:s} to aperture {:s}'
                          ' in CCD {:s}'.format(
                              self._link_apnam, apnam, self._link_cnam))
                    PickStar.action_prompt(True)

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
                PickStar.action_prompt(True)

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
                PickStar.action_prompt(True)

    def _reference(self):
        """
        Toggles the reference status of an aperture
        """

        # first see if there is an aperture near enough the selected position
        aper, apnam, dmin = self._find_aper()

        if dmin is None or \
           dmin > max(self.rtarg,min(100,max(20.,2*self.rsky2))):
            print('  *** found no aperture near to the cursor position')
        else:

            if aper.is_linked():
                print('  *** a linked aperture cannot become a '
                      'reference aperture')
            else:
                aper.ref = not aper.ref
                if aper.ref:
                    print('  aperture {:s} is now a '
                          'reference aperture'.format(apnam))
                else:
                    print('  aperture {:s} is no longer a '
                          'reference aperture'.format(apnam))

                # remove aperture from plot
                for obj in self.pobjs[self._cnam][apnam]:
                    obj.remove()

                # re-plot new version, over-writing plot objects
                self.pobjs[self._cnam][apnam] = hcam.mpl.pAper(
                    self._axes, aper, apnam)

        PickStar.action_prompt(True)

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
                self._add_mode = False

        elif key == '!' and self._buffer == '':
            # terminate accumulation mode without bothering to wait for an 'enter'
            print('\n*** no aperture added')
            PickStar.action_prompt(True)
            self._add_mode = False

        elif key == 'backspace' or key == 'delete':
            # remove a character 
            self._buffer = self._buffer[:-1]
            print('{:s}{:s} '.format(PickStar.ADD_PROMPT, self._buffer))

        elif key in '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ':
            # accumulate input and add to the buffer
            self._buffer += key

            if self.linput == 's':
                # single character input. bail out immediately (inside _add_aperture)
                print(key)

                if key == '0':
                    print("Apertures cannot be labelled just '0'")
                    print(PickStar.ADD_PROMPT, end='',flush=True)
                    self._buffer = ''

                elif self._buffer in self.mccdaper[self._cnam]:
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
                    self._add()
                    self._add_mode = False

            else:
                # multi character input. just accumulate characters
                print(key, end='', flush=True)


