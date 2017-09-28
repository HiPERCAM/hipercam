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

      input  : (string)
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

    """

    if args is None:
        args = sys.argv[1:]

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', 'setaper', args) as cl:

        # register parameters
        cl.register('input', Cline.LOCAL, Cline.PROMPT)
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
        frame = cl.get_value('input', 'frame to plot',
                             cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(frame)

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

    # define the aperture picker
    ccdaper = hcam.CcdAper()
    picker = PickStar(cnams, toolbar, fig, ccdaper)

    plt.tight_layout()
    print('a(dd), r(emove) p(an/zoom), h(elp), q(uit)')

    # once this goes we enter into a loop where we can select stars
    # using the 
    plt.show()


class PickStar:
    """Class to pick targets for apertures.
    """

    def __init__(self, cnams, toolbar, fig, rtarg, rsky1, rsky2, ccdaper):
        """
        axes is dictionary keyed by CCD label of the axes
        """
        self.cnams = cnams
        self.toolbar = toolbar
        self.fig = fig
        self.fig.canvas.mpl_connect('key_press_event', self._keyPressEvent)
        self.rtarg = rtarg
        self.rsky1 = rsky1
        self.rsky2 = rsky2
        self.ccdaper = ccdaper

    def _keyPressEvent(self, event):
        """
        This is where we do most of the hard work
        """
        if self.toolbar.mode != 'pan/zoom' and event.inaxes is not None:
            cnam = self.cnams[event.inaxes]
            if event.key == 'h':
                print('a(dd)      : add an aperture')
                print('p(an/zoom) : toggle between pan/zoom vs selection mode')
                print('r(emove)   : remove an aperture')
                print('h(elp)     : print this list')

            elif event.key == 'a':
                # add a new aperture
                while 1:
                    label = input("enter a label for the aperture, '!' to abort: ")
                    if label in ccdaper:
                        print(
                            'label={:s} already in use; please try again'.format(label),
                            file=sys.stderr
                            )
                    elif label == '':
                        print(
                            'label blank; please try again'.format(label),
                            file=sys.stderr
                            )
                    elif label != '!':
                        # add & plot aperture
                        aper = hcam.Aperture(event.x, event.y, rtarg, rsky1, rsky2)
                        ccdaper[label] = aper
                        hcam.mpl.pAper(event.inaxes, aper)
                        
                print('will add an aperture to CCD',cnam)
            elif event.key == 'p':
                pass
            elif event.key == 'q':
                # quit and clear up
                plt.close()
                print('apertures saved. bye')
            elif event.key == 'r':
                print('will remove an aperture from CCD',cnam)
            elif event.key == 'enter':
                print('a(dd), r(emove) p(an/zoom), h(elp), q(uit)')
            else:
                print('There is no action for key =',event.key,'(x,y,axes =',
                      event.x,event.y,event.inaxes,')')

