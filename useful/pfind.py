"""
Script to search and measure periodic signals in ZTF data downloaded
from IRSA in ipac format. Exact format could easily be changed; see
script.

Does the following:

1) Loads data
2) Computes LS periodogram over defined frequency range
3) Plots it (blue solid), indicates frequency of peak power with
   green dashed line
4) Optionally plot reference frequency with red dashed line
5) Also allows injection of artificial signals
6) Fits a sinsuoid
7) If successful, will subtract sinusoid, re-compute LS power and
   plot as a blue dashed line.

Run as "python pfind.py <arguments>". "python pfind.py -h" will print
out some help.

Uses standard Python + numpy, scipy, matplotlib and astropy.

Written by T.R.Marsh, Warwick.
"""

import re
import argparse
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pylab as plt
from astropy.timeseries import LombScargle
from astropy.table import Table
from astropy import units as u

def model(p, ts):
    """
    Implements model of a sinsuoid plus constant, i.e.

    y = c + a*cos(2.*Pi*f*(t-T0))

    Arguments:

      p : array of 4 parameters, c,a,T0,f
      ts : numpy array of times

    Returns: equivalent y
    """
    c,a,t0,f = p
    return c + a*np.cos(2.*np.pi*f*(ts-t0))

def res(p, ts, ms, mes):
    """
    Computes residual vector (data-model)/uncertainty
    for sinusoid fitting. See "model". Adds magnitudes
    and their errors to list of arguments.
    """
    c,a,t0,f = p
    return (ms-model(p,ts))/mes

def jac(p, ts, ms, mes):
    """
    Computes partial derivatives with respect to parameters
    equivalent to "res" to help sinusoid fitting
    """
    c,a,t0,f = p
    phase = 2.*np.pi*f*(ts-t0)
    dc = np.ones_like(ts)
    da = np.cos(phase)
    fixed = -2.*np.pi*a*np.sin(phase)
    dt0 = -fixed*f
    df = fixed*(ts-t0)
    return np.column_stack(
        [-dc/mes,-da/mes,-dt0/mes,-df/mes]
    )

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=""" Searches for a periodic signal in a ZTF data file given a
        range of frequencies, fits a sinusoid and reports parameters.  """
    )

    parser.add_argument(
        'ztf', help='Data file from ZTF'
    )
    parser.add_argument(
        'flo', type=float,
        help='lower frequency limit, cycles/day'
    )
    parser.add_argument(
        'fhi', type=float,
        help='upper frequency limit, cycles/day'
    )
    parser.add_argument(
        '--filter', default='zr|zg',
        help='filter codes to analyse. Separate filters with | if more than one'
    )
    parser.add_argument(
        '--over', type=float, default=1.,
        help='oversampling factor'
    )
    parser.add_argument(
        '--fref', type=float, default=0.,
        help='reference frequency, cycles/day'
    )
    parser.add_argument(
        '--ffake', type=float, default=0.,
        help='frequency of artificial injection signal, cycles/day'
    )
    parser.add_argument(
        '--afake', type=float, default=0.05,
        help='amplitude of artificial injection signal, mags'
    )

    args = parser.parse_args()

    tab = Table.read(args.ztf, format='ascii.ipac')

    # extract data of interest
    ts = tab['hjd']
    ms = tab['mag']
    mes = tab['magerr']
    filts = tab['filtercode']
    ok = filts == filts

    # Report the filter codes loaded
    print(f'Loaded {len(ts)} points from {args.ztf}')
    ufilts = np.unique(filts)
    for filt in ufilts:
        ok = filts == filt
        print(f'{len(ts[ok])} point have filter code "{filt}"')

    # Identify the ones of interest
    refilt = re.compile(f'^{args.filter}$')
    for n, filt in enumerate(filts):
        if not refilt.match(filt):
            ok[n] = False

    print(f'{len(ts[ok])} points matched the search code(s) "{args.filter}"')
    print(f'Oversampling factor = {args.over}')
    if args.fref > 0:
        print(f'Will show frequency reference [red dashes] at {args.fref} cycles/day')

    # refine the data
    ts,ms,mes = np.array(ts[ok]),np.array(ms[ok]),np.array(mes[ok])

    if args.ffake > 0:
        # signal injection
        print(
            'Will inject artificial signal with frequency = ' +
            f'{args.ffake} cycles/day and amplitude {args.afake} mags'
        )
        ms += args.afake*np.cos(2.*np.pi*args.ffake*(ts-ts.mean()))

    if len(ts) < 6:
        print(f'Too few points ({len(ts)}) for valid operation')
        exit(1)
    elif len(ts) < 20:
        print(f'Very few points ({len(ts)}); unlikely to work well')

    # compute Lomb-Scargle periodogram
    tbase = ts.max()-ts.min()
    nf = int(4*args.over*tbase*(args.fhi-args.flo))
    print(f'Searching {nf} frequencies')

    freq = np.linspace(args.flo,args.fhi,nf)
    power = LombScargle(ts,ms,mes).power(freq)

    # Report highest peak
    fmax = freq[power.argmax()]
    print(f'Maximum power found at frequency = {fmax} cycles/day')

    # Sinusoid fit. When fitting, it's best to shift the time
    # zeropoint for numerical reasons
    tref = ts.mean()
    ts -= tref
    x0 = (ms.mean(),ms.std(),0.,fmax)
    results = least_squares(res, x0, jac, method='lm', args=(ts,ms,mes))
    print(results.x,results.cost)
    J = np.matrix(results.jac)
    c,a,t0,f = results.x
    try:
        # Report fit values, subtract fit and re-compute
        # power spectrum. Main peak should be removed if
        # ok. Trap errors from bar covariances.

        # scale uncertainties to give chi**2/dof = 1.
        sfac = np.sqrt(results.cost/(len(ts)-4))
        ce,ae,t0e,fe = sfac*np.sqrt(np.diag((J.T * J).I))
        print('Best fit sinusoid parameters:')
        print(f' c  = {c} +/- {ce} mags')
        print(f' a  = {a} +/- {ce} mags')
        print(f' t0 = {tref+t0} +/- {t0e}')
        print(f' f  = {f} +/- {fe}')
        mms = ms - model(results.x, ts)
        powersub = LombScargle(ts,mms,mes).power(freq)
        plt.plot(freq,powersub,'b--')

    except (np.linalg.LinAlgError, ValueError) as err:
        print('Failed to compute uncertainties; probably a poor fit')
        print('Will not attempt to subtract fitted sinusoid')
        print(f'Error returned: {err}')

    # Plot power in raw data
    plt.plot(freq,power,'b')

    if args.fref > 0.:
        # Indicate reference frequency
        plt.axvline(args.fref,ls='--',color='r')

    # Indicate frequency of maximum power
    plt.axvline(fmax,ls='--',color='g')

    # Finish up
    plt.xlabel('Frequency [cycles/day]')
    plt.ylabel('Power')
    plt.title(f'Lomb-Scargle periodogram of {args.ztf}')
    plt.show()


