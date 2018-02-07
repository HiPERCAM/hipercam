import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

###########################################################
#
# pfolder -- plots phase-folded versions of reduce log data
#
###########################################################

def pfolder(args=None):
    """``pfolder log [device width height] ccd aper``

    Folds data from reduce logs on user specified ephemeris (in MJD)
    and plots the results. Rather specific routine to help with timing
    tests where the LED switches on on the second and off on the half second.

    Parameters:

      log     : string
          name of reduce log file (text file with loads of columns)

      device  : string [hidden, defaults to 'term']
         'term' for interactive plot, file name such as 'plot.pdf'
         for a hardcopy.

      width   : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height  : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH
         width AND height must be non-zero to have any effect

      ccd     : string
         the CCD to consider, e.g. '1'

      aper    : string
         the aperture to consider

      t0      : float
         zero point of the ephemeris in MJD (i.e. no light travel correction)

      period  : float
         period to fold on [seconds]

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('log', Cline.LOCAL, Cline.PROMPT)
        cl.register('device', Cline.LOCAL, Cline.HIDE)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('aper', Cline.LOCAL, Cline.PROMPT)
        cl.register('t0', Cline.LOCAL, Cline.PROMPT)
        cl.register('period', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        log = cl.get_value(
            'log', 'reduce log file to plot',
            cline.Fname('hcam', hcam.LOG)
        )

        device = cl.get_value('device', 'plot device name', 'term')
        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

        ccd = cl.get_value('ccd', 'first CCD to plot', '1')
        aper = cl.get_value('aper', 'first aperture', '1')
        t0 = cl.get_value('t0', 'zero point of ephemeris [MJD]', 55000.)
        period = cl.get_value('period', 'period of ephemeris [seconds]',
                              1., 1e-6)

    # load the reduce log
    hlog = hcam.hlog.Hlog.fromLog(log)

    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width,height))
    else:
        fig = plt.figure()

    # load counts data, fold, plot two cycles
    data = hlog.tseries(ccd, aper, 'counts')
    data.t = np.mod(86400.*(data.t-t0)/period,1)
    xlabel = 'Phase [cycles]'
    ylabel = 'Counts'
    data.mplot(plt)
    data.t += 1
    data.mplot(plt)

    plt.title('{:s}, CCD {:s}, Aperture {:s}'.format(log,ccd,aper))

    if device == 'term':
        plt.show()
    else:
        plt.savefig(device)
