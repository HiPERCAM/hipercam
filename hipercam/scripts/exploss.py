import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

##########################################
#
# exploss -- computes exposure loss factor
#
##########################################


def exploss(args=None):
    """``exploss log [device] aper readout gain``

    Computes the equivalent exposure time needed to match the
    signal-to-noise ratio of the current input log file as a fraction
    of the actual exposure time, assuming zero readout noise. The
    purpose of the routine is to help in judging how readout noise
    limited a given run is. The result is a loss factor between 0 and
    1. If close to 1, readout noise is not significant. As well as the
    loss factor for the supplied run (plotted in blue), the case for
    double the exposure time is plotted (in red). If this is
    sinificantly larger than the blue values, it may be advisable to
    increase the exposure time or add to NBLUE / NSKIP.

    Parameters:

      log : string
          name of |reduce| ASCII log file (text file with loads of columns)

      device : string [hidden, defaults to 'term']
         'term' for interactive plot, file name such as 'plot.pdf'
         for a hardcopy.

      aper : str
         aperture to consider

      readout : float
         readout noise, RMS ADU.

      gain : float
         Gain,  electrons/ADU.

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("log", Cline.LOCAL, Cline.PROMPT)
        cl.register("device", Cline.LOCAL, Cline.HIDE)
        cl.register("aper", Cline.LOCAL, Cline.PROMPT)
        cl.register("readout", Cline.LOCAL, Cline.PROMPT)
        cl.register("gain", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        log = cl.get_value(
            "log", "reduce log file to plot", cline.Fname("hcam", hcam.LOG)
        )

        device = cl.get_value("device", "plot device name", "term")
        aper = cl.get_value("aper", "aperture", "1")
        readout = cl.get_value("readout", "readout noise [ADU]", 2.5)
        gain  = cl.get_value("gain", "gain [electrons/ADU]", 1.0)

    # load the reduce log
    hlog = hcam.hlog.Hlog.rascii(log)

    ccds = []
    for ccd in hlog:
        if aper in hlog.apnames[ccd]:
            ccds.append(ccd)

    fig,axs = plt.subplots(len(ccds),1,sharex=True)

    bitmask = hcam.BAD_TIME|hcam.JUNK

    for ccd,ax in zip(ccds,axs):
        tims = hlog[ccd][f'MJD']
        obj = hlog[ccd][f'counts_{aper}']
        nsky = hlog[ccd][f'nsky_{aper}']
        sky = hlog[ccd][f'sky_{aper}']

        objsky = obj + nsky*sky
        lfac = objsky/(nsky*readout**2 + objsky)
        lfac2 = 2*objsky/(nsky*readout**2 + 2*objsky)
        ax.plot(tims,lfac,'.b')
        ax.plot(tims,lfac2,'.r')
        ax.set_title(f'CCD {ccd}')
        ax.set_ylabel('Loss factor')
        ax.axhline(1,ls='--',color='k')
        ax.set_ylim(0,1.1)
        
    ax.set_xlabel('Time [MJD]')
 
    if device == "term":
        plt.show()
    else:
        plt.savefig(device)
