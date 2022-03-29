#!/usr/bin/env python

import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

from astropy.timeseries import LombScargle

def pbands(args=None):
    """``pbands log [device width height] norm zero [plo phi increment]
    pgram (flo fhi over) title``

    Plots a HiPERCAM |reduce| log as a single light curve in one panel
    per CCD. Can optionally be normalised and the y-axis origin set to
    zero, and can also include periodogram plots if wanted. The
    purpose of the routine is to facilitate quick-look plots during
    observing runs.

    Parameters:

      log : string
          name of |reduce| ASCII log file (text file with loads of columns)

      device : string [hidden, defaults to 'term']
         'term' for interactive plot, file name such as 'plot.pdf'
         or 'object.png' for a hardcopy.

      width : float [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : float [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH
         width AND height must be non-zero to have any effect

      norm : bool
         normalise each light curve to a unit median (or not)

      zero : bool
         set y-axis origin to zero

      plo : float [hidden]
         percentile to use to set the lower limit of the lightcurve plots
         (e.g. 5). Takes account of error bars. 

      phi : float [hidden]
         percentile to use to set the upper limit of the lightcurve plot
         (e.g. 95). Must be > plo. Takes account of error bars.

      increment : float [hidden]
         multiplier to increment the raw numbers derived via plo and phi.
         If they are y1 and y2, then each is corrected up and down by
         increment*(y2-y1).

      pgram : bool
         plot periodogram panels (or not). They appear on the right-hand
         side if yes.

      flo : float (if pgram)
         If pgram, the mimimum frequency for the periodogram plots,
         cycles per day. Typically there is little point pushing flo
         very far below 1 cycle over the timespan of the data.

      fhi : float (if pgram)
         If pgram, the maximum frequency for the periodogram plots
         cycles per day. The maximum useful is set by the Nyquist
         limit (1 cycle per 2 exposures). e.g. For 10 second cadence,
         this would be be 86400/10/2 = 4320 cycles/day.

      over : float (if pgram)
         If pgram, the oversampling factor to determine the frequency grid.
         Corresponds to the "sample_per_peak" parameter of the astropy
         LombScargle routine.

      title : str
         plot title. Defaults to name of log file

      colours : list of strings
         colours for each CCD, space or comma-separated.
         e.g. "r,g,b" or "r g b" (without quotes).

    """

    command, args = utils.script_args(sys.argv)
    if command.endswith('.py'):
        command = command[:-3]

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("log", Cline.LOCAL, Cline.PROMPT)
        cl.register("device", Cline.LOCAL, Cline.HIDE)
        cl.register("width", Cline.LOCAL, Cline.HIDE)
        cl.register("height", Cline.LOCAL, Cline.HIDE)
        cl.register("aper1", Cline.LOCAL, Cline.PROMPT)
        cl.register("aper2", Cline.LOCAL, Cline.PROMPT)
        cl.register("norm", Cline.LOCAL, Cline.PROMPT)
        cl.register("zero", Cline.LOCAL, Cline.PROMPT)
        cl.register("plo", Cline.LOCAL, Cline.HIDE)
        cl.register("phi", Cline.LOCAL, Cline.HIDE)
        cl.register("increment", Cline.LOCAL, Cline.HIDE)
        cl.register("pgram", Cline.LOCAL, Cline.PROMPT)
        cl.register("flo", Cline.LOCAL, Cline.PROMPT)
        cl.register("fhi", Cline.LOCAL, Cline.PROMPT)
        cl.register("over", Cline.LOCAL, Cline.PROMPT)
        cl.register("title", Cline.LOCAL, Cline.PROMPT)
        cl.register("colours", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        log = cl.get_value(
            "log", "reduce log file to plot", cline.Fname("hcam", hcam.LOG)
        )
        # load the reduce log
        hlog = hcam.hlog.Hlog.rascii(log)

        cl.set_default("device","term")
        device = cl.get_value("device", "plot device name", "term")
        width = cl.get_value("width", "plot width (inches)", 0.0)
        height = cl.get_value("height", "plot height (inches)", 0.0)
        aper1 = cl.get_value("aper1", "target aperture", "1")
        aper2 = cl.get_value("aper2", "comparison aperture ('!' to ignore)", "2")

        ccds = []
        for ccd in hlog.apnames:
            if aper1 in hlog.apnames[ccd] and \
               (aper2 == '!' or aper2 in hlog.apnames[ccd]):
                ccds.append(ccd)

        if len(ccds) == 0:
            raise hcam.HipercamError(
                f'Failed to find any CCDs with apertures {aper1} and {aper2}'
            )

        norm = cl.get_value("norm", "normalise to median", True)
        zero = cl.get_value("zero", "set y-axis origin to zero", True)
        plo = cl.get_value("plo", "percentile for plot lower limit", 0.0, 0., 100.)
        phi = cl.get_value("phi", "percentile for plot upper limit", max(100.0,plo), plo, 100.)
        increment = cl.get_value("increment", "plot range increment factor", 0.01, 0.)
        pgram = cl.get_value("pgram", "plot periodogram panels too", False)
        if pgram:
            flo = cl.get_value("flo", "lower frequency limit (cycles/day)", 5.0, 0.)
            fhi = cl.get_value(
                "fhi", "upper frequency limit (cycles/day)", max(flo,2000.), flo
            )
            over = cl.get_value("over", "oversampling factor", 5., 1.)

        if norm:
            cl.set_default("title", f'{log} [median normalised]')
        else:
            cl.set_default("title", f'{log}')
        title = cl.get_value("title", "plot title", "Plot Title")

        # Need to modify default according to numbers of CCDs
        colours = cl.get_default("colours",('r','g','b'))
        if len(colours) > len(ccds):
            colours = colours[:len(ccds)]
            cl.set_default("colours",colours)
        elif len(colours) < len(ccds):
            colours = tuple(colours) + (len(ccds)-len(colours))*('k',)
            cl.set_default("colours",colours)
        colours = cl.get_value(
            "colours", "colours [strings, one per CCD]",
            ('r','g','b')
        )

    if width > 0 and height > 0:
        fig,axs = plt.subplots(
            len(ccds),2 if pgram else 1,figsize=(width, height),
            sharex='col', squeeze=False
        )
    else:
        fig,axs = plt.subplots(
            len(ccds),2 if pgram else 1,
            sharex='col',squeeze=False
        )

    # light curves
    for ccd,ax,col in zip(ccds,axs[:,0],colours):
        targ = hlog.tseries(ccd, aper1, 'counts')
        if aper2 != "!":
            comp = hlog.tseries(ccd, aper2, 'counts')
            targ /= comp
        if norm:
            targ.normalise()
        targ.mplot(
            ax, col, bitmask=hcam.BAD_TIME|hcam.JUNK,
        )

        # plot limits
        _d,y1,_d = targ.percentile(plo,hcam.BAD_TIME|hcam.JUNK)
        _d,_d,y2 = targ.percentile(phi,hcam.BAD_TIME|hcam.JUNK)
        yrange = y2-y1
        y1 -= increment*yrange
        y2 += increment*yrange

        if zero:
            ax.set_ylim(0,y2)
        else:
            ax.set_ylim(y1,y2)
        if aper2 == "!":
            ax.set_ylabel('Targ')
        else:
            ax.set_ylabel('Targ/Comp')

    if pgram:
        # periodograms
        for ccd,ax,col in zip(ccds,axs[:,1],colours):
            targ = hlog.tseries(ccd, aper1, 'counts')
            if aper2 != "!":
                comp = hlog.tseries(ccd, aper2, 'counts')
                targ /= comp

            # extract data
            ts,_d,ys,yes = targ.get_data(hcam.BAD_TIME|hcam.JUNK)
            ls = LombScargle(ts,ys,yes)
            f,p = ls.autopower(
                minimum_frequency=flo,maximum_frequency=fhi,
                samples_per_peak=over
            )
            ax.plot(f,p,col)
            ax.set_ylabel('Power')

        # plot limits
        _d,_d,y1 = targ.percentile(plo,hcam.BAD_TIME|hcam.JUNK)
        _d,_d,y2 = targ.percentile(phi,hcam.BAD_TIME|hcam.JUNK)
        yrange = y2-y1
        y1 -= increment*yrange
        y2 += increment*yrange

        if zero:
            ax.set_ylim(0,y2)
        else:
            ax.set_ylim(y1,y2)

    if pgram:
        axs[0,0].set_title(f'{title} [light curves]')
        axs[0,1].set_title(f'{title} [periodograms]')
        axs[-1,0].set_xlabel('Time [MJD]')
        axs[-1,1].set_xlabel('Frequency [cycles/day]')
    else:
        axs[0,0].set_title(f'{title}')
        axs[-1,0].set_xlabel('Time [MJD]')

    plt.tight_layout()
    if device == "term":
        plt.show()
    else:
        plt.savefig(device)
