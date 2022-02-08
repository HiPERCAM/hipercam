#!/usr/bin/env python

"""
Script for measuring extinction starting from a hipercam pipeline log file.
Not general; needs adapting. See comments. See part following "if __name__"
etc for main bits to adapt. The log file should have been created using a large
aperture scale factor (~2.5) to grab most of the flux.

Written by TRM, 2022-02-08
"""

import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import least_squares
import hipercam as hcam

def model(p, airms):
    # Linear model extinction in mags vs airmass
    return p[0] + p[1]*airms

def res(p, airms, ms, mes):
    return (ms-model(p,airms))/mes

def rejfit(airms, ms, mes, thresh):
    """
    Carries out fit with one-by-one rejection.
    Returns with little vector of fit coeffients,
    array of ok points [logical, True=not-rejected],
    and the total number of points rejected.

    "thresh" is maximum deviation in terms of sigma,
    scaled by sqrt(chisq/ndof) of fit, accounting for
    previously rejected points.
    """
    ok = ms == ms
    x = [-12,0.2]
    indices = np.arange(len(ok),dtype=int)
    nrej = 0
    while True:
        results = least_squares(res,x,args=(airms[ok],ms[ok],mes[ok]),method='lm')
        x = results.x
        fit = model(x,airms)
        dev = (ms-fit)/mes
        chisq = (dev[ok]**2).sum()
        ndof = len(dev[ok])-2
        imax = np.abs(dev[ok]).argmax()
        dmax = dev[ok][imax]
        if dmax > thresh*np.sqrt(chisq/ndof):
            ok[indices[ok][imax]] = False
            nrej += 1
        else:
            break
    return (x,ok,nrej)

if __name__ == '__main__':

    # Target position, ICRS
    position = '16:21:17.36 +44:12:54.1'

    # site / telescope [astropy]
    telescope = 'lapalma'

    # The log file to analyse
    hlog = 'run0014-extinct.log'

    # The apertures containing comparison stars to show,
    # along with colours to plot them, and offsets for residual
    # plots in mags
    comps = (('1','b',0.05),('2','g',-0.05))

    # CCDs to analyse
    ccds = ('1','2','3','4','5')

    # Amount to bin by to clean data and speed up computations
    nbin = 20

    # Threshold to use to attempt to clean data
    thresh = 3.0

    # Load data
    hlg = hcam.hlog.Hlog.rascii(hlog)


    # Should be mostly fixed from here

    # Create plot
    fig,axs = plt.subplots(len(ccds),2,sharex=True,figsize=(8,14))

    for n,ccd in enumerate(ccds):
        for aper,col,yoff in comps:
            # extract light curve, convert fluxes to mags, and times
            # to airmass
            lc = hlg.tseries(ccd,aper)
            lc.bin(nbin)
            lc.to_mag()
            lc.to_airmass(position,telescope)

            # fit linear extinction with rejection
            airm = np.array(lc.t)
            x,ok,nrej = rejfit(airm,lc.y,lc.ye,thresh)

            # Flag rejected ppoints
            lc.set_bitmask(hcam.CLOUDS,~ok)

            # Plot without rejected points
            lc.mplot(axs[n,0],col,bitmask=hcam.CLOUDS)
            fit = model(x,airm)
            axs[n,0].plot(airm,fit,'r--',zorder=10)

            # Report results
            print(f'CCD {ccd}, aperture {aper}: constant,gradient = {x}, nrej = {nrej}')

            # Plot residuals
            lc -= fit
            (lc+yoff).mplot(axs[n,1],col,bitmask=hcam.CLOUDS)
            axs[n,1].axhline(yoff,ls='--',color='r')

        # Invert y-axis limits (mags)
        y1,y2=axs[n,0].set_ylim()
        axs[n,0].set_ylim(y2,y1)
        axs[n,0].set_ylabel(f'CCD {ccd} [mags]')
        y1,y2=axs[n,1].set_ylim()
        axs[n,1].set_ylim(y2,y1)

    # Finish off
    axs[len(axs)-1,0].set_xlabel('Airmass')
    axs[len(axs)-1,1].set_xlabel('Airmass')
    plt.tight_layout()
    plt.show()

