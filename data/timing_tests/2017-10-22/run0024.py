#!/usr/bin/env python

"""
Script to calculate and plot timing offset and exposure length
of run0024 from 2017-10-22, a no clear mode, 2-windowed run with
4x4 binning.

Simply invoke with no arguments and it will generate plots for 
each CCD of interest with names like '2017-10-22-0024-4.png'
"""

import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import hipercam as hcam

if __name__ == '__main__':

    # load the reduce log
    log = 'run0024.log'
    hlog = hcam.hlog.Hlog.fromLog(log)

    # range of times to include in analysis and
    # in the final plot, fractional difference
    # to identify on / off
    TRANGE, PRANGE, FDIFF = 0.03, 0.015, 0.01

    # CCDs to plot
    ccds = ('4','5')

    aper = '1'
    for ccd in ccds:

        fig = plt.figure(figsize=(8,5))

        # load counts data
        data = hlog.tseries(ccd, aper, 'counts')

        # Zero point = start of day
        T0 = int(data.t[0])

        # Fold on a 1-second period
        data.t = np.mod(86400.*(data.t-T0),1)

        t = np.concatenate((data.t-1,data.t))
        c = np.concatenate((data.y,data.y))
        ce = np.concatenate((data.ye,data.ye))

        ok = (t > -TRANGE) & (t < TRANGE)
        t = t[ok]
        c = c[ok]
        ce = ce[ok]

        # First estimates of LED 'off' and 'on' levels
        off, on = np.percentile(c,(5,95))

        # split into 3 classes
        MAXDIFF = FDIFF*(on-off)

        low = c - off < MAXDIFF
        low = t <= t[low].max()
        high = on - c < MAXDIFF
        high = t >= t[high].min()
        mid = ~low & ~high

        # better estimate of on and off levels
        off = c[low].mean()
        on  = c[high].mean()

        # Linear fit to transition points
        grad, const = np.polyfit(t[mid],c[mid],1,w=1/ce[mid])
        toff = (off-const)/grad
        ton  = (on-const)/grad
        tmid = ((off+on)/2-const)/grad
        plt.plot([toff,ton],[off,on],'--',color='0.7')
        plt.plot([tmid,tmid],[off,on],'--',color='0.7')

        # plot
        plt.errorbar(t[low],c[low],ce[low],fmt='.r',ms=1)
        plt.errorbar(t[mid],c[mid],ce[mid],fmt='.g',ms=1)
        plt.errorbar(t[high],c[high],ce[high],fmt='.b',ms=1)

        plt.xlabel('1-second phase')
        plt.ylabel('Counts in aperture')
        plt.title(
            ('2017-10-17/{:s}, CCD {:s},'
             ' off = {:+d} $\mu$s, exp = {:.5f} s').format(
                 log[:-4], ccd, int(round(1.e6*tmid)),
                 ton-toff)
        )
        plt.xlim(-PRANGE, PRANGE)
#        plt.show()
        pname = '2017-10-22-{:s}-{:s}.png'.format(log[3:-4],ccd)
        plt.savefig(pname)
