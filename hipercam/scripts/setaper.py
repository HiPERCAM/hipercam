import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from matplotlib.backend_bases import NavigationToolbar2

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
    """Interactive definition of photometric extraction apertures

    Arguments::

      mccd   : (string)
         name of MCCD file

      ccd    : (string)
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      rtarg  : (float) [unbinned pixels]
         radius of target aperture. The exact value of this does not matter too much since
         it is normally overridden in 'reduce', but typically one aims for 1.5 to 2.5 x FWHM,
         seeing, depending upon the target brightness.

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

      xlo    : (float)
         left X-limit, for PGPLOT and matplotlib harcopy plots only as the
         matplotlib interactive plot is zoomable

      xhi    : (float)
         right X-limit

      ylo    : (float)
         bottom Y-limit

      yhi    : (float)
         top Y-limit

      width  : (float) [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : (float) [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

    Notes. There are a few conveniences to make setaper easier::

      1. You can toggle between the pan/zoom and input modes using the 'Alt' 
         button as well as the usual icon in the toolbar. 

      2. Clicking any mouse button acts like hitting 'a' for adding an aperture.

    Various standard keyboard shortcuts are disabled as they just confuse things
    and are of limited use in setaper in any case.
    """

    if args is None:
        args = sys.argv[1:]

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', 'setaper', args) as cl:

        # register parameters
        cl.register('mccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('aper', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
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
        cl.register('xlo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('xhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ylo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('yhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)

        # get inputs
        mccd = cl.get_value('mccd', 'frame to plot',
                            cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(mccd)

        aper = cl.get_value('aper', 'frame to plot',
                            cline.Fname('hcam', hcam.APER, exist=False))

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

        xlo, xhi, ylo, yhi = 0, nxmax+1, 0, nymax+1

        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)


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

    # all inputs obtained, plot
    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width,height))
    else:
        fig = plt.figure()

    # get the navigation toolbar which is used to check the pan/zoom mode
    toolbar = fig.canvas.manager.toolbar

    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    ax = None
    pstars = {}
    cnams = {}
    for n, cnam in enumerate(ccds):
        if ax is None:
            axes = ax = fig.add_subplot(ny, nx, n+1)
            axes.set_aspect('equal', adjustable='box')
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

    # define the multi apertures and the picker
    ccdaper = hcam.CcdAper()
    picker = PickStar(cnams, toolbar, fig, ccdaper, rtarg, rsky1, rsky2, aper)

    plt.tight_layout()
    print('a(dd), d(elete), h(elp), Q(uit)')

    # finally show stuff ....
    plt.show()


# the next class is really where all the action occurs. A rather complicated
# matter of handling events. note that standard terminal input with 'input'
# becomes impossible, explaining some of the weirdness

class PickStar:
    """Class to pick targets for apertures.
    """

    def __init__(self, cnams, toolbar, fig, ccdaper, rtarg, rsky1, rsky2, apernam):
        """
        axes is dictionary keyed by CCD label of the axes. all the args are
        stored as attributes to allow "global"-type access to a call back
        added to fig.canvas to handle KeyPress events from the user.
        """
        self.cnams = cnams
        self.toolbar = toolbar
        self.fig = fig
        self.fig.canvas.mpl_connect('key_press_event', self._keyPressEvent)
        self.fig.canvas.mpl_connect('button_press_event', self._buttonPressEvent)
        self.ccdaper = ccdaper
        self.rtarg = rtarg
        self.rsky1 = rsky1
        self.rsky2 = rsky2
        self.rsky2 = rsky2
        self.apernam = apernam

        # the next are to do with entering / leaving accumulation mode and storing 
        # the state at the point when we enter accumulation mode
        self._accum = False
        self._action = None
        self._buffer = None
        self._cnam = None
        self._key = None
        self._x = None
        self._y = None
        self._axes = None

    def _keyPressEvent(self, event):
        """
        This is where we do the hard work. Every key press event is diverted
        to this method. It either takes an action based on the input, such as
        removing an Aperture, or sometimes it causes a state change such that
        input is diverted to and accumulated in a buffer until 'enter' is hit.
        The latter stage comes first.
        """
        ADD_PROMPT = "enter a label for the aperture, '!' to abort: "

        if self._accum:
            # accumulation mode action
            self._accumulate(event.key)

        else:
            # standard mode action
            self._accumulate(event.key, event.xdata, event.ydata, event.inaxes)


    def _buttonPressEvent(self, event):
        """Treat a button press event as equivalent to 'enter' as a convenience"""
        if self._accum:
            # accumulation mode action, pretend we have hit 'enter'
            self._accumulate('enter')

        else:
            # standard mode action, pretend we have hit 'a' for add
            self._standard('a', event.xdata, event.ydata, event.inaxes)

    def _accumulate(self, key):
        """
        Carries out the work needed when we are in accumulation mode. Just pass through
        the value of the key pressed and this will handle the rest. It defines how we
        treat characters.
        """

        if key == 'enter':
            # trap 'enter'
            print()
            if self._action == 'a':
                if self._buffer in self.ccdaper:
                    print(
                        'label={:s} already in use; please try again'.format(self._buffer),
                        file=sys.stderr
                        )
                    print(ADD_PROMPT, end='',flush=True)
                    self._buffer = ''

                elif self._buffer == '':
                    print(
                        'label blank; please try again'.format(self._buffer),
                        file=sys.stderr
                        )
                    print(ADD_PROMPT, end='',flush=True)

                else:
                    # add & plot aperture
                    aper = hcam.Aperture(
                        self._x, self._y, self.rtarg, self.rsky1, self.rsky2, False)
                    self.ccdaper[self._buffer] = aper

                    # add aperture to the plot
                    hcam.mpl.pAper(self._axes, aper)

                    # make sure it appears
                    self.fig.canvas.draw()

                    print('added aperture {:s} to CCD {:s}'.format(self._buffer,self._cnam))
                    print('\na(dd), d(elete), h(elp), Q(uit)')

                    # terminate accumulation mode
                    self._accum = False

        elif key == '!' and self._buffer == '':
            # terminate accumulation mode without bothering to wait for an 'enter'
            print('\n*** no aperture added')
            print('\na(dd), d(elete), h(elp), Q(uit)')
            self._accum = False

        elif key == 'backspace' or key == 'delete':
            # remove a character 
            self._buffer = self._buffer[:-1]
            print('{:s}{:s} '.format(ADD_PROMPT, self._buffer))

        elif key in '!0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ':
            # accumulate input and add to the
            self._buffer += key
            print(key, end='', flush=True)

    def _standard(self, key, x, y, axes):
        """
        Carries out the work needed when we are in the standard mode. Just pass through
        the value of the key pressed and this will handle the rest. It defines how we
        treat characters.
        """

        if key == 'alt':
            if self.toolbar.mode == 'pan/zoom':
                print('switched to pan/zoom mode')
            else:
                print('switched to aperture update mode')
                print('a(dd), d(elete), h(elp), Q(uit)')

        elif self.toolbar.mode != 'pan/zoom' and axes is not None:

            # save information which is useful for accumulation mode
            # with the axes
            self._cnam = self.cnams[axes]
            self._key = key
            self._axes = axes
            self._axes = axes


            if key == 'h':
                print('a(dd)      : add an aperture')
                print('p(an/zoom) : toggle between pan/zoom vs selection mode')
                print('d(elete)   : delete an aperture')
                print('h(elp)     : print this list')
                print('Q(uit)     : quite setaper')

            elif key == 'a':
                print(ADD_PROMPT, end='',flush=True)
                self._accum = True
                self._buffer = ''
                self._action = 'a'

            elif key == 'p':
                pass

            elif key == 'd':
                print('will delete aperture from CCD',cnam)

            elif key == 'Q':
                # quit and clear up
                plt.close()
                print('apertures saved. bye')

            elif key == 'enter':
                print('a(dd), r(emove) p(an/zoom), h(elp), Q(uit)')

            else:
                print('There is no action defined for key =',key,'(xdata,ydata,axes =',
                      x,y,axes,')')
