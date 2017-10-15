#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import time

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
from numpy.lib import recfunctions
from sklearn.cluster import KMeans
import scipy.interpolate as interpolate
import scipy.spatial as sp
import sep
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clip
from astropy.convolution import Gaussian2DKernel

import hipercam as hcam
from trm.pgplot import *
import hipercam.cline as cline
from hipercam.cline import Cline
from hipercam.extraction import findStars


def findBestRigidTransform(x, y, xref, yref):
    """
    Given a set of points and *matching* reference points, find the optimal Rigid transform.

    A Rigid transform is a rotation and translation, i.e:

    B = R.A + t

    See https://scicomp.stackexchange.com/questions/6878/fitting-one-set-of-points-to-another-by-a-rigid-motion
    or http://nghiaho.com/?page_id=671

    for method.
    """

    A = np.column_stack((xref, yref))
    B = np.column_stack((x, y))

    # find centre of points
    centroidA = np.mean(A, axis=0)
    centroidB = np.mean(B, axis=0)

    # shift centres to origin
    Ashifted = A - centroidA
    Bshifted = B - centroidB

    # find dot product
    H = np.dot(Bshifted.T, Ashifted)

    # use SVD to find rotation
    u, s, v = np.linalg.svd(H)
    R = np.dot(v, u.T)
    if np.linalg.det(R) < 0:
        R = np.dot(v, np.array([1, -1])*u.T)

    T = np.median(B - np.dot(A, R), axis=0)
    theta1 = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
    theta2 = np.degrees(-np.arctan2(R[0, 1], R[1, 1]))
    theta = 0.5*(theta1 + theta2)
    return theta, T


def matchSpots(x, y, xref, yref):
    """
    Match pairs of spots, looking for the nearest neighbour.

    Parameters
    ----------
    x, y: `numpy.ndarray`
        spots to match
    xref, yref: `numpy.ndarray`
        reference positions

    Returns
    -------
    theta, offset: `numpy.ndarray`
        rotation and offset that matches points to reference points
    """
    tree = sp.cKDTree(np.column_stack((xref, yref)))
    seps, indices = tree.query(np.column_stack((x, y)))
    if np.all(seps < 0.01):
        # probably the same data
        mask = np.zeros_like(seps).astype('bool')
    else:
        sep_sigma = (seps-seps.mean())/seps.std()
        mask = sep_sigma < 1.5
    return findBestRigidTransform(x[mask], y[mask], xref[indices][mask], yref[indices][mask])


def displayFWHM(xpos, ypos, xref, yref, fwhm, peak, fix_fwhm_scale=False, bins=20):

    mask = ~sigma_clip(peak).mask
    if len(mask) > 0 and len(mask) == len(xpos):
        xpos = xpos[mask]
        ypos = ypos[mask]
        fwhm = fwhm[mask]
        peak = peak[mask]

    # remove correlation between flux and spot fwhm
    fitOK = True
    try:
        f = np.poly1d(np.polyfit(peak, fwhm, 1))
        meanF = np.mean(f(peak))
        fwhm_prime = fwhm - f(peak) + meanF
    except:
        fitOK = False

    imageOK = True
    try:
        xi = np.linspace(xpos.min(), xpos.max(), bins)
        yi = np.linspace(ypos.min(), ypos.max(), bins)
        xx, yy = np.meshgrid(xi, yi)
        func = interpolate.Rbf(xpos, ypos, fwhm_prime, function='quintic')
        image = func(xx, yy)
    except:
        imageOK = False

    # definitions for the axes
    left, width = 0.05, 0.38
    bottom, height = 0.1, 0.5
    right = left + width + 0.02
    left_h = right + width + 0.02
    bottom_h = bottom + height + 0.02

    rect_im = [right, bottom, width, height]
    rect_xcut = [right, bottom_h, width, 0.2]
    rect_ycut = [left_h, bottom, 0.1, height]
    rect_pos = [left, bottom, width, height]
    rect_corr = [left, bottom_h+0.1, width, 0.1]

    axIm = plt.axes(rect_im)
    axX = plt.axes(rect_xcut)
    axY = plt.axes(rect_ycut)
    axPos = plt.axes(rect_pos, sharex=axIm, sharey=axIm)
    axCorr = plt.axes(rect_corr)

    nullfmt = NullFormatter()
    axY.yaxis.set_major_formatter(nullfmt)
    axX.xaxis.set_major_formatter(nullfmt)
    axX.yaxis.set_label_position("right")
    axX.yaxis.tick_right()

    if imageOK:
        axIm.pcolormesh(xx, yy, image, cmap='magma')
    axX.plot(xpos, fwhm, 'r.')
    axY.plot(fwhm, ypos, 'r.')
    if fitOK:
        axX.plot(xpos, fwhm_prime, 'b+')
        axY.plot(fwhm_prime, ypos, 'b+')
        axCorr.plot(peak, fwhm, 'r.')
        idx = np.argsort(peak)
        axCorr.plot(peak[idx], f(peak[idx]), 'b--')

    axPos.plot(xpos, ypos, 'r.')
    axPos.plot(xref, yref, 'b+', alpha=0.5)

    axPos.set_xlabel('x')
    axPos.set_ylabel('y')
    axIm.set_xlabel('x')
    axX.set_ylabel('FWHM')
    axY.set_xlabel('FWHM')
    if fix_fwhm_scale:
        axY.set_xlim(1.9, 2.8)
        axX.set_ylim(1.9, 2.8)
    theta, offset = matchSpots(xpos, ypos, xref, yref)
    xoff, yoff = offset
    textstr = r'$\theta = %.1f$; $\mathrm{xoff}=%.1f$, $\mathrm{yoff}=%.1f$' % (theta, xoff, yoff)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axPos.text(left, 1.02, textstr, horizontalalignment='left',
               verticalalignment='bottom', bbox=props, fontsize=10,
               transform=axPos.transAxes)


def measureSpots(mccd, ccdname, thresh, use_hfd=True):
    """
    Detect spots in all windows of a ccd.
    """
    ccd = mccd[ccdname]
    objects = np.concatenate([findStars(wind, thresh, 3) for (wname, wind) in ccd.items()])
    xpos = objects['x']
    ypos = objects['y']
    fwhm = objects['fwhm']
    hfd = objects['hfd']

    if use_hfd:
        mask = np.logical_and(np.isfinite(hfd).data, objects['peak'] < 20000)
        return xpos[mask], ypos[mask], hfd[mask], objects['cflux'][mask]
    else:
        mask = fwhm <= 20
        mask = np.logical_and(mask, objects['peak'] < 20000)
        return xpos[mask], ypos[mask], fwhm[mask], objects['cflux'][mask]


def separateBigSmallSpots(xpos, ypos, fwhm, flux):
    """
    Use K-Means clustering to seperate large and small starfield tube spots based on
    FWHM and flux.

    Arguments
    ---------
    xpos, ypos : np.ndarray
        spot locations
    fwhm : np.ndarray
        fwhm of spots
    flux : np.ndarray
        flux contained within spots

    Returns
    -------
    lx, ly, lfw, lfl, bx, by, bfw, bfl : np.ndarray
        positions, fwhm and flux of small and large spots, respectively
    """
    X = np.array((fwhm, flux)).T
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(X)
    size_ids = kmeans.predict(X)
    spot_mask0 = size_ids == 0
    spot_mask1 = size_ids == 1
    if fwhm[spot_mask0].mean() > fwhm[spot_mask1].mean():
        big_spot_mask = spot_mask0
        little_spot_mask = spot_mask1
    else:
        big_spot_mask = spot_mask1
        little_spot_mask = spot_mask0
    return (xpos[little_spot_mask], ypos[little_spot_mask], fwhm[little_spot_mask], flux[little_spot_mask],
            xpos[big_spot_mask], ypos[big_spot_mask], fwhm[big_spot_mask], flux[big_spot_mask])


def main(args=None):
    """
    Measures alignment of HiperCAM CCDs w.r.t a reference image of the blue CCD.

    Arguments::

        ref     : (string)
           name of an MCCD file to use as reference , as produced
           by e.g. 'grab'

        rccd    : (string)
           CCD to use as reference

        thresh  : (float)
            threshold for object detection, in multiples of background RMS

        source  : (string) [hidden]
           's' = server, 'l' = local, 'f' = file list (ucm for ULTRACAM /
           ULTRASPEC, FITS files for HiPERCAM)

        inst    : (string) [hidden]
           If 's' = server, 'l' = local, name of instrument. Choices: 'u' for
           ULTRACAM / ULTRASPEC, 'h' for HiPERCAM. This is needed because of the
           different formats.

        device  : (string) [hidden]
          Plot device. PGPLOT is used so this should be a PGPLOT-style name,
          e.g. '/xs', '1/xs' etc. At the moment only ones ending /xs are
          supported.

        width   : (float) [hidden]
           plot width (inches). Set = 0 to let the program choose.

        height  : (float) [hidden]
           plot height (inches). Set = 0 to let the program choose. BOTH width
           AND height must be non-zero to have any effect

        run     : (string) [if source == 's' or 'l']
           run number to access, e.g. 'run034'

        flist   : (string) [if source == 'f']
           name of file list

        first   : (int) [if source='s' or 'l']
           exposure number to start from. 1 = first frame; set = 0 to
           always try to get the most recent frame (if it has changed)

        twait   : (float) [if source == 's'; hidden]
           time to wait between attempts to find a new exposure, seconds.

        tmax    : (float) [if source == 's'; hidden]
           maximum time to wait between attempts to find a new exposure, seconds.

        pause   : (float) [hidden]
           seconds to pause between frames (defaults to 0)

        ccd     : (string)
           CCD to measure aligment of.

        bias    : (string)
           Name of bias frame to subtract, 'none' to ignore.

        msub    : (bool)
           subtract the median from each window before scaling for the
           image display or not. This happens after any bias subtraction.

        iset    : (string) [single character]
           determines how the intensities are determined. There are three
           options: 'a' for automatic simply scales from the minimum to the
           maximum value found on a per CCD basis. 'd' for direct just takes
           two numbers from the user. 'p' for percentile dtermines levels
           based upon percentiles determined from the entire CCD on a per CCD
           basis.

        ilo     : (float) [if iset='d']
           lower intensity level

        ihi     : (float) [if iset='d']
           upper intensity level

        plo     : (float) [if iset='p']
           lower percentile level

        phi     : (float) [if iset='p']
           upper percentile level

        small_spots : (bool) [hidden]
            use the small spots for alignment instead of the big spots
    """
    if args is None:
        args = sys.argv[1:]

    with Cline('HIPERCAM_ENV', '.hipercam', 'rtplot', args) as cl:
        # register parameters
        cl.register('ref', Cline.LOCAL, Cline.PROMPT)
        cl.register('thresh', Cline.LOCAL, Cline.PROMPT)
        cl.register('small_spots', Cline.LOCAL, Cline.HIDE)
        cl.register('source', Cline.LOCAL, Cline.HIDE)
        cl.register('inst', Cline.GLOBAL, Cline.HIDE)
        cl.register('device', Cline.LOCAL, Cline.HIDE)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('twait', Cline.LOCAL, Cline.HIDE)
        cl.register('tmax', Cline.LOCAL, Cline.HIDE)
        cl.register('flist', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('rccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('pause', Cline.LOCAL, Cline.HIDE)
        cl.register('bias', Cline.GLOBAL, Cline.PROMPT)
        cl.register('msub', Cline.GLOBAL, Cline.PROMPT)
        cl.register('iset', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ilo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ihi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('plo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('phi', Cline.LOCAL, Cline.PROMPT)
        cl.register('xlo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('xhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ylo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('yhi', Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        reference = cl.get_value('ref', 'reference frame to align to',
                                 cline.Fname('hcam', hcam.HCAM))

        ref_mccd = hcam.MCCD.rfits(reference)

        thresh = cl.get_value('thresh', 'detection threshold (sigma)', 5.0)
        small_spots = cl.get_value('small_spots', 'use small spots for alignment', False)

        source = cl.get_value('source', 'data source [s(erver), l(ocal), f(ile list)]',
                              'l', lvals=('s', 'l', 'f'))
        if source == 's' or source == 'l':
            inst = cl.get_value('inst', 'instrument [h(ipercam), u(ltracam/spec)]',
                                'h', lvals=('h', 'u'))

        # plot device stuff
        device = cl.get_value('device', 'plot device', '1/xs')
        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

        if source == 's' or source == 'l':
            run = cl.get_value('run', 'run name', 'run005')
            first = cl.get_value('first', 'first frame to plot', 1, 1)

            if source == 's':
                twait = cl.get_value('twait', 'time to wait for a new frame [secs]', 1., 0.)
                tmax = cl.get_value('tmax', 'maximum time to wait for a new frame [secs]', 10., 0.)

        else:
            # set inst = 'h' as only lists of HiPERCAM files are supported
            inst = 'h'
            run = cl.get_value('flist', 'file list', cline.Fname('files.lis',hcam.LIST))
            first = 1

        flist = source == 'f'
        server = source == 's'
        if inst == 'u':
            instrument = 'ULTRA'
        elif inst == 'h':
            instrument = 'HIPER'

        # get the labels and maximum dimensions
        ccdinf = hcam.get_ccd_pars(instrument, run, flist)
        if len(ccdinf) > 1:
            ccdnam = cl.get_value('ccd', 'CCD to plot alignment of', '1')
            rccdnam = cl.get_value('rccd', 'CCD to use as reference', '5')

        cl.set_default('pause', 0.)
        pause = cl.get_value('pause', 'time delay to add between frame plots [secs]', 0., 0.)

        # bias frame (if any)
        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        if bias is not None:
            # read the bias frame
            bframe = hcam.MCCD.rfits(bias)

        msub = cl.get_value('msub', 'subtract median from each window?', True)

        iset = cl.get_value(
            'iset', 'set intensity a(utomatically), d(irectly) or with p(ercentiles)?',
            'a', lvals=['a', 'A', 'd', 'D', 'p', 'P'])
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == 'd':
            ilo = cl.get_value('ilo', 'lower intensity limit', 0.)
            ihi = cl.get_value('ihi', 'upper intensity limit', 1000.)
        elif iset == 'p':
            plo = cl.get_value('plo', 'lower intensity limit percentile', 5., 0., 100.)
            phi = cl.get_value('phi', 'upper intensity limit percentile', 95., 0., 100.)

        # region to plot
        nxmax, nymax = 0, 0
        nxtot, nytot = ccdinf[ccdnam]
        nxmax = max(nxmax, nxtot)
        nymax = max(nymax, nytot)

        xlo = cl.get_value('xlo', 'left-hand X value', 0., 0., float(nxmax+1))
        xhi = cl.get_value('xhi', 'right-hand X value', float(nxmax), 0., float(nxmax+1))
        ylo = cl.get_value('ylo', 'lower Y value', 0., 0., float(nymax+1))
        yhi = cl.get_value('yhi', 'upper Y value', float(nymax), 0., float(nymax+1))

    # arguments defined, let's do something!
    # open image plot device
    imdev = hcam.pgp.Device(device)
    if width > 0 and height > 0:
        pgpap(width, height/width)

    pgsubp(1, 2)
    for i in range(2):
        pgsci(hcam.pgp.Params['axis.ci'])
        pgsch(hcam.pgp.Params['axis.number.ch'])
        pgenv(xlo, xhi, ylo, yhi, 1, 0)
        pglab('X', 'Y', ('reference', 'data')[i])

    total_time = 0  # time waiting for new frame

    with hcam.data_source(instrument, run, flist, server, first) as spool:
        # 'spool' is an iterable source of MCCDs
        for n, frame in enumerate(spool):
            # None objects are returned from failed server reads. This could
            # be because the file is still exposing, so we hang about.
            if server and frame is None:
                if tmax < total_time + twait:
                    print('Have waited for {:.1f} sec. cf tmax = {:.1f}; will wait no more'.format(total_time, tmax))
                    print('alignment tool stopped.')
                    break

                print('Have waited for {:.1f} sec. cf tmax = {:.1f}; will wait another twait = {:.1f} sec.'.format(
                        total_time, tmax, twait
                        ))

                # pause
                time.sleep(twait)
                total_time += twait
                # have another go
                continue
            elif server:
                # reset the total time waited when we have a success
                total_time = 0

            # indicate progress
            if flist:
                print('File {:d}: '.format(n+1), end='')
            else:
                print('Frame {:d}: '.format(frame.head['NFRAME']), end='')

            ccd = frame[ccdnam]
            ref_ccd = ref_mccd[rccdnam]

            if bias is not None:
                ccd -= bframe[ccdnam]

            if msub:
                for wind in ccd.values():
                    wind -= wind.median()
                for wind in ref_ccd.values():
                    wind -= wind.median()

            # plot images
            message = ''
            for i in range(2):
                this_ccd = (ref_ccd, ccd)[i]
                name = ('ref: ', 'data: ')[i]
                pgpanl(1, i+1)
                vmin, vmax = hcam.pgp.pCcd(this_ccd, iset, plo, phi, ilo, ihi)
                message += ' {:s}: {:.2f} to {:.2f}'.format(name, vmin, vmax)
            print(message)

            # time for measurement of spots!
            try:
                xpos, ypos, fwhm, flux = measureSpots(frame, ccdnam, thresh)
                lx, ly, lf, lp, bx, by, bf, bp = separateBigSmallSpots(xpos, ypos, fwhm, flux)
                if not small_spots:
                    # use big spots
                    x, y, f, flux = bx, by, bf, bp
                else:
                    x, y, f, flux = lx, ly, lf, lp
            except Exception as err:
                print('Failed to find spots in frame: ', end=' ')
                print(str(err))
                continue

            # also measure spots in reference image
            try:
                xpos, ypos, fwhm, rflux = measureSpots(ref_mccd, ccdnam, thresh)
                lx, ly, lf, lp, bx, by, bf, bp = separateBigSmallSpots(xpos, ypos, fwhm, rflux)
                if not small_spots:
                    # use big spots
                    xref, yref, fref = bx, by, bf
                else:
                    xref, yref, fref = lx, ly, lf
            except Exception as err:
                print('Failed to find spots in reference: ', end=' ')
                print(str(err))
                continue

            plt.clf()
            displayFWHM(x, y, xref, yref, f, flux, fix_fwhm_scale=False)
            plt.draw()
            plt.pause(0.05)

    plt.ioff()
    plt.show()
