import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = ['hist',]

##########################################
#
# hist -- plots histograms of MCCD objects
#
##########################################

def hist(args=None):
    """Plots histograms of :class:`MCCD` objects. Histograms can be powerful
    diagnostics of CCD problems.

    Arguments::

      input  : string
         name of the MCCD file

      ccd    : string
         CCD or CCDs to plot in histogram form, one CCD per panel. Can be '0'
         for all of them or a specific set '2 3 5'

      window : string
         the window label to show. '0' will merge all windows on each CCD.

      x1     : float
         lower end of histoigram range

      x2     : float
         upper end of histoigram range. Set == x1 to get automatic min to max.

      nbins  : int
         number of bins to use for the histogram

      gplot  : bool
         plot a gaussian of the same mean and RMS as the data. Note
         this is not a fit.

      nx     : int
         number of panels across to display, prompted if more than one CCD is
         to be plotted.

      msub   : bool
         True/False to subtract median from each window before scaling

      nint   : bool [if msub]
         True/False to take the nearest integer of the median or not. Can help
         with histograms of digitised data.

      device : string [hidden]
         Plot device name. Either 'term' for an interactive plot or a
         name like 'plot.pdf' for a hard copy.

      width  : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('input', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('window', Cline.LOCAL, Cline.PROMPT)
        cl.register('x1', Cline.LOCAL, Cline.PROMPT)
        cl.register('x2', Cline.LOCAL, Cline.PROMPT)
        cl.register('nbins', Cline.LOCAL, Cline.PROMPT)
        cl.register('gplot', Cline.LOCAL, Cline.PROMPT)
        cl.register('nx', Cline.LOCAL, Cline.PROMPT)
        cl.register('msub', Cline.GLOBAL, Cline.PROMPT)
        cl.register('nint', Cline.GLOBAL, Cline.PROMPT)
        cl.register('device', Cline.LOCAL, Cline.HIDE)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)


        # get inputs
        frame = cl.get_value('input', 'frame to plot',
                             cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.read(frame)

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

        wnam = cl.get_value('window', 'window to plot [0 for all]', '0')

        x1 = cl.get_value('x1', 'left-hand value of histogram', 0.)
        x2 = cl.get_value('x2', 'right-hand value of histogram', 0.)
        if x1 != x2:
            range = (x1,x2)
        else:
            range = None
        nbins = cl.get_value('nbins', 'number of bins for the histogram', 100, 2)
        gplot = cl.get_value('gplot', 'overplot gaussians?', True)

        if len(ccds) > 1:
            nxdef = min(len(ccds), nxdef)
            cl.set_default('nx', nxdef)
            nx = cl.get_value('nx', 'number of panels in X', 3, 1)
        else:
            nx = 1

        msub = cl.get_value('msub', 'subtract median from each window?', True)
        if msub:
            nint = cl.get_value('nint', 'take nearest integer median?', True)

        device = cl.get_value('device', "plot device name ['term' for interactive]", 'term')
        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

    # OK, inputs obtained, plot
    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width,height))
    else:
        fig = plt.figure()

    mpl.rcParams['xtick.labelsize'] = hcam.mpl.Params['axis.number.fs']
    mpl.rcParams['ytick.labelsize'] = hcam.mpl.Params['axis.number.fs']

    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    ax = None
    for n, cnam in enumerate(ccds):
        if ax is None:
            axes = ax = fig.add_subplot(ny, nx, n+1)
        else:
            axes = fig.add_subplot(ny, nx, n+1, sharex=ax, sharey=ax)

        ccd = mccd[cnam]

        if msub:
            # subtract median from each window
            for wind in ccd.values():
                wmed = wind.median()
                if nint:
                    wmed = round(wmed)
                wind -= wmed

        if wnam == '0':
            arr = ccd.flatten()
        else:
            arr = ccd[wnam].flatten()

        n, bins, ps = plt.hist(arr,bins=nbins,range=range)

        num = len(arr)
        mean = arr.mean()
        std = arr.std()
        if wnam == '0':
            plt.title('CCD = {:s}, $\sigma = {:.2f}$'.format(cnam,std))
        else:
            plt.title('CCD = {:s}, window = {:s}, $\sigma = {:.2f}$'.format(
                cnam,wnam,std)
            )

        if gplot:
            x1, x2 = bins.min(), bins.max()
            x = np.linspace(x1, x2, 500)
            y = abs(x2-x1)/nbins*num*np.exp(
                -((x-mean)/std)**2/2.)/np.sqrt(2.*np.pi)/std
            plt.plot(x,y,'r--',lw=0.5)

    plt.ylabel('Number of pixels')
    plt.xlabel('Counts')
    plt.tight_layout()
    if device == 'term':
        plt.show()
    else:
        plt.savefig(device)

