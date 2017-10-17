import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

##########################################
#
# plog -- simple plots of reduce log files
#
##########################################

def plog(args=None):
    """Provides quick-look plots of HiPERCAM reduce logs.

    Arguments::

      log    : (string)
          name of reduce log file (text file with loads of columns)

      device : (string) [hidden, defaults to 'term']
         'term' for interactive plot, file name such as 'plot.pdf'
         for a hardcopy.

      width  : (float) [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : (float) [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect

      ccd    : (string)
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. 

      targ   : (string)
         target aperture

      comp   : (string)
         comparison aperture, 'none' to ignore.

      param  : (string)
         parameter to plot

      scheme : (string) [if 'comp' not none]
         how to plot 'targ' and 'comp'. 'd' = difference;
         'i' = ignore; 'r' = ratio; 's' = scatter plot.
    """

    if args is None:
        args = sys.argv[1:]

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', 'plog', args) as cl:

        # register parameters
        cl.register('log', Cline.LOCAL, Cline.PROMPT)
        cl.register('device', Cline.LOCAL, Cline.HIDE)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('targ', Cline.LOCAL, Cline.PROMPT)
        cl.register('comp', Cline.LOCAL, Cline.PROMPT)
        cl.register('param', Cline.LOCAL, Cline.PROMPT)
        cl.register('scheme', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        log = cl.get_value(
            'log', 'reduce log file to plot',
            cline.Fname('hcam', hcam.LOG)
        )

        device = cl.get_value('device', 'plot device name', 'term')
        print('device = ',device)
        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

        ccd = cl.get_value('ccd', 'CCD(s) to plot [0 for all]', '0')

        targ = cl.get_value('targ', 'target aperture', '1')

        comp = cl.get_value(
            'comp', "comparison aperture ['!' to ignore]", '2'
        )

        param = cl.get_value(
            'param', 'parameter [x, y, f(whm), b(eta), c(ounts), s(ky)]',
            'c', lvals=('x','y','f','b','c','s')
        )
        MAP = {
            'x' : 'x',
            'y' : 'y',
            'f' : 'fwhm',
            'b' : 'beta',
            'c' : 'counts',
            's' : 'sky'
        }
        pname = MAP[param]

        if comp != '!':
            scheme = cl.get_value(
                'scheme', 'd(ifference), i(gnore), r(atio), s(catter)',
                'i', lvals=('d','i','r','s')
            )

    # load the reduce log
    hlog = hcam.hlog.Hlog.fromLog(log)

    if ccd == '0':
        ccds = list(hlog.keys())
    else:
        ccds = ccd.split()

    if width > 0 and height > 0:
        fig = plt.figure(figsize=(width,height))
    else:
        fig = plt.figure()

#    mpl.rcParams['xtick.labelsize'] = hcam.mpl.Params['axis.number.fs']
#    mpl.rcParams['ytick.labelsize'] = hcam.mpl.Params['axis.number.fs']

    for cnam in ccds:
        tdat = hlog.tseries(cnam, targ, pname)
        if comp != '!':
            cdat = hlog.tseries(cnam, comp, pname)
            if scheme == 'i':
                tdat.mplot(plt)
                cdat.mplot(plt)

            elif scheme == 'r':
                ratio = tdat / cdat
                ratio.mplot(plt)

            elif scheme == 'd':
                diff = tdat - cdat
                diff.mplot(plt)

            elif scheme == 's':
                ok = (tdat.ye > 0) & (cdat.ye > 0)
                plt.errorbar(
                    tdat.y[ok], cdat.y[ok], cdat.ye[ok], tdat.ye[ok],
                    '.', capsize=0
                )
        else:
            tdat.mplot(plt)

    if device == 'term':
        plt.show()
    else:
        plt.savefig(device)
