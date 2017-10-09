import sys
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
import configparser

from trm.pgplot import *
import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ = ['reduce',]

################################################
#
# reduce -- reduces multi-CCD imaging photometry
#
################################################

def reduce(args=None):
    """Reduces a sequence of multi-CCD images, plotting lightcurves as they
    come in.

    reduce can source data from both the ULTRACAM and HiPERCAM servers, from
    local 'raw' ULTRACAM and HiPERCAM files (i.e. .xml + .dat for ULTRACAM, 3D
    FITS files for HiPERCAM) and from lists of HiPERCAM '.hcm' files.

    reduce is primarily configured from a file, typically "reduce.red". This
    is closely related to the format of the ULTRACAM files except it uses
    ConfigParser which allows sections

    Arguments::

        source  : (string) [hidden]
           's' = server, 'l' = local, 'f' = file list (ucm for ULTRACAM /
           ULTRASPEC, FITS files for HiPERCAM)

        inst    : (string) [hidden]
           If 's' = server, 'l' = local, name of instrument. Choices: 'u' for
           ULTRACAM / ULTRASPEC, 'h' for HiPERCAM. This is needed because of the
           different formats.

        ldevice  : (string) [hidden]
          Plot device for the light curves. PGPLOT is used so this should be a PGPLOT-style name,
          e.g. '/xs', '1/xs' etc. At the moment only ones ending /xs are supported.

        lwidth  : (float) [hidden]
           plot width (inches). Set = 0 to let the program choose.

        lheight : (float) [hidden]
           plot height (inches). Set = 0 to let the program choose. BOTH width
           AND height must be non-zero to have any effect

        rfile   : (string)
          the "reduce" file, i.e. ASCII text file suitable for reading by ConfigParser. Best
          seen by example as it has many parts.

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

    """

    if args is None:
        args = sys.argv[1:]

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'rtplot', args) as cl:

        # register parameters
        cl.register('source', Cline.LOCAL, Cline.HIDE)
        cl.register('inst', Cline.GLOBAL, Cline.HIDE)
        cl.register('ldevice', Cline.LOCAL, Cline.HIDE)
        cl.register('lwidth', Cline.LOCAL, Cline.HIDE)
        cl.register('lheight', Cline.LOCAL, Cline.HIDE)
        cl.register('rfile', Cline.GLOBAL, Cline.PROMPT)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('twait', Cline.LOCAL, Cline.HIDE)
        cl.register('tmax', Cline.LOCAL, Cline.HIDE)
        cl.register('flist', Cline.LOCAL, Cline.PROMPT)

        # get inputs

        # image plot
        source = cl.get_value('source', 'data source [s(erver), l(ocal), f(ile list)]',
                              'l', lvals=('s','l','f'))

        if source == 's' or source == 'l':
            inst = cl.get_value('inst', 'instrument [h(ipercam), u(ltracam/spec)]',
                                'h', lvals=('h','u'))

        # plot device stuff
        ldevice = cl.get_value('ldevice', 'light curve plot device', '1/xs')
        lwidth = cl.get_value('lwidth', 'light curve plot width (inches)', 0.)
        lheight = cl.get_value('lheight', 'light curve plot height (inches)', 0.)

        # the reduce file
        rfilen = cl.get_value('rfile', 'reduce file', cline.Fname('reduce.red',hcam.RED))
        rfile = Rfile.fromFile(rfilen)
        print(rfile)

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

        # define the panel grid. first get the labels and maximum dimensions
        ccdinf = hcam.get_ccd_pars(instrument, run, flist)

    ################################################################
    #
    # all the inputs have now been obtained. Get on with doing stuff

    # open the light curve plot
    lcdev = hcam.pgp.Device(ldevice)
    if lwidth > 0 and lheight > 0:
        pgpap(lwidth,lheight/lwidth)

    pgsci(hcam.pgp.Params['axis.ci'])
    pgsch(hcam.pgp.Params['axis.number.ch'])
    pgenv(0, 10, 0, 10, 1, 0)
    pglab('Time [mins]','Mag','')

    # a couple of initialisations
    total_time = 0 # time waiting for new frame
    fpos = [] # list of target positions to fit

    # get frames
    with hcam.data_source(instrument, run, flist, server, first) as spool:

        # 'spool' is an iterable source of MCCDs
        for nframe, mccd in enumerate(spool):

            # None objects are returned from failed server reads. This could
            # be because the file is still exposing, so we hang about.
            if server and mccd is None:

                if tmax < total_time + twait:
                    print('Have waited for {:.1f} sec. cf tmax = {:.1f}; will wait no more'.format(total_time, tmax))
                    print('reduce stopped.')
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
                print('Frame {:d}: '.format(mccd.head['NFRAME']), end='')

            if nframe == 0:
                # This is the very first data frame read in. We need to get the 
                rfile.coerce()

            # do something else ...
            mccdwins = {}
            for cnam in mccd:
                # get the CCD and the apertures
                ccd = mccd[cnam]
                ccdaper = rfile.aper[cnam]
                if nframe == 0:
                    # first time through, work out an array of which window each aperture lies in
                    # we will assume this is fixed for the whole run, i.e. that apertures do not
                    # drift from one window to another. Set to None if no window found
                    mccdwins[cnam] = {}
                    for apnam, aper in ccdaper.items():
                        for wnam, wind in ccd.items():
                            if wnam.distance(aper.x,aper.y) > 0:
                                mccdwins[cnam][apnam] = wnam
                                break
                        else:
                            mccdwins[cnam][apnam] = None

                    store ={}

                ccdwins = mccdwins[cnam]

                if isinstance(rfile.read, hcam.MCCD):
                    read = rfile.read[cnam]
                else:
                    read = rfile.read

                if isinstance(rfile.gain, hcam.MCCD):
                    gain = rfile.gain[cnam]
                else:
                    gain = rfile.gain

                # at this point 'ccd' contains all the Windats of a CCD, ccdaper all of its apertures
                # ccdwins the label of the Windat relevant for each aperture, rfile contains some
                # control parameters, read contains the readout noise, either a float or a CCD,
                # gain contains the gain, either a float or a CCD.

                moveApers(cnam, ccd, ccdaper, ccdwins, rfile, read, gain, store)
                print('frame=',nframe,'\n',ccdaper,'\n')

# From here is support code not visible outside

class Rfile(configparser.ConfigParser):
    """Class to read and interpret reduce files."""

    # Data
    VERSION = '2017-10-05'

    @classmethod
    def fromFile(self, filename):

        # read it in, stripping off trailing comments
        rfile = configparser.ConfigParser(inline_comment_prefixes='#')
        with open(filename) as fp:
            rfile.read_file(fp)

        # process it

        # the version
        if rfile['general']['version'] != Rfile.VERSION:
            # check the version
            raise ValueError(
                'Version mismatch: file = {:s}, reduce = {:s}'.format(
                    rfile['general']['version'], Rfile.VERSION)
            )

        # the apertures
        self.aper = hcam.MccdAper.fromJson(
            hcam.add_extension(rfile['apertures']['aperfile'],hcam.APER)
            )
        if rfile['apertures']['fit_method'] == 'moffat':
            self.method = 'm'
        elif  rfile['apertures']['fit_method'] == 'gaussian':
            self.method = 'g'
        else:
            raise NotImplementedError('apertures.fit_method = {:s} not recognised'.format(
                    rfile['apertures']['fit_method']))

        if rfile['apertures']['reposition'] != 'yes' and rfile['apertures']['reposition'] != 'no':
            raise ValueError("apertures.reposition: 'yes' or 'no' are the only supported values")

        if rfile['apertures']['fit_fwhm_fixed'] != 'yes' and rfile['apertures']['fit_fwhm_fixed'] != 'no':
            raise ValueError("apertures.reposition: 'yes' or 'no' are the only supported values")

        # the calibrations
        if rfile['calibration']['bias'] != '':
            self.bias == hcam.MCCD.rfits(
                hcam.add_extension(rfile['calibration']['bias'],hcam.HCAM)
                )
        else:
            self.bias = None

        if rfile['calibration']['dark'] != '':
            self.dark == hcam.MCCD.rfits(
                hcam.add_extension(rfile['calibration']['dark'],hcam.HCAM)
                )
        else:
            self.dark = None

        if rfile['calibration']['flat'] != '':
            self.flat == hcam.MCCD.rfits(
                hcam.add_extension(rfile['calibration']['flat'],hcam.HCAM)
                )
        else:
            self.flat = None

        try:
            self.readout = float(rfile['calibration']['readout'])
        except TypeError:
            self.readout = hcam.MCCD.rfits(
                hcam.add_extension(rfile['calibration']['readout'],hcam.HCAM)
                )

        try:
            self.readout = float(rfile['calibration']['gain'])
        except TypeError:
            self.gain = hcam.MCCD.rfits(
                hcam.add_extension(rfile['calibration']['gain'],hcam.HCAM)
                )

        return rfile

    def coerce(self, mccd):
        """
        This uses a template file 'mccd' to try to get the calibration files into the same format.
        """
        if self.bias is not None:
            self.bias = self.bias.crop(mccd)

        if self.dark is not None:
            self.dark = self.dark.crop(dark)

        if self.flat is not None:
            self.flat = self.flat.crop(flat)

        if isinstance(self.read, hcam.MCCD):
            self.read = self.read.crop(read)

        if isinstance(self.gain, hcam.MCCD):
            self.read = self.gain.crop(read)

def moveApers(cnam, ccd, ccdaper, ccdwin, rfile, read, gain, store):
    """
    Encapsulates aperture re-positioning. 'store'
    is a dictionary of results that will be used to start
    the fits from one frame to the next. It must start
    as {}.
    """

    if rfile['apertures']['reposition'] == 'no':
        # do nothing
        return

    # first of all try to get a mean shift from the reference apertures.
    xsum, ysum = 0., 0.
    wxsum, wysum = 0., 0.
    ref = False

    for apnam, aper:
        if aper.ref:
            ref = True
            # extract Windat for this reference aperture
            wind = ccd[ccdwin[apnam]]

            if isinstance(read, hcam.CCD):
                rd = read[ccdwin[apnam]].data
            else:
                rd = read

            if isinstance(gain, hcam.CCD):
                gn = gain[ccdwin[apnam]].data
            else:
                gn = gain

            # get sub-windat around start position
            shbox = int(rfile['apertures']['search_half_width_ref'])
            swind = wind.window(aper.x-shbox, aper.x+shbox, aper.y-shbox, aper.y+shbox)

            # carry out initial search
            x,y,peak = swind.find(float(rfile['apertures']['search_smooth_fwhm']), False)

            # now for a more refined fit. First extract fit Windat
            fhbox = int(rfile['apertures']['fit_half_width'])
            fwind = wind.window(x-fhbox, x+fhbox, y-fhbox, y+fhbox)
            sky = np.percentile(fwind.data, 25)

            # get some parameters from previous run where possible
            if apnam in store and store[apnam]['efwhm'] > 0.:
                fit_fwhm = store[apnam]['fwhm']
            else:
                fit_fwhm = float(rfile['apertures']['fit_fwhm'])

            if apnam in store and store[apnam]['ebeta'] > 0.:
                fit_beta = store[apnam]['beta']
            else:
                fit_beta = float(rfile['apertures']['fit_beta'])

            # refine the Aperture position by fitting the profile
            try:
                (sky, height, x, y, fwhm, beta), \
                    (esky, eheight, ex, ey, efwhm, ebeta), \
                    extras = hcam.combFit(
                    fwind, rfile.method, sky, peak-sky,
                    x, y, fit_fwhm,
                    float(rfile['apertures']['fit_fwhm_min']),
                    rfile['apertures']['fit_fwhm_fixed'] == 'yes',
                    fit_beta, rd, gn
                    )

                if height > float(rfile['apertures']['height_min']):
                    dx = x - aper.x
                    wx = 1./ex**2
                    wxsum += wx
                    xsum += wx*dx

                    dy = y - aper.y
                    wy = 1./ey**2
                    wysum += wy
                    ysum += wy*dy

                    # store some stuff for next time
                    store[apnam] = {'fwhm' : fwhm, 'efwhm' : efwhm,
                                    'beta' : beta, 'ebeta' : ebeta}
                else:
                    print('CCD {:s}, reference aperture {:s}, peak = {:.1f} < {:s}'.format(
                            cnam, apnam, height, rfile['apertures']['height_min']),
                          file=sys.stderr)
                    store[apnam] = {'efwhm' : -1, 'ebeta' : -1}


            except hcam.HipercamError as err:
                print('CCD {:s}, reference aperture {:s}, fit failed'.format(cnam, apnam), file=sys.stderr)
                store[apnam] = {'efwhm' : -1, 'ebeta' : -1}

    if ref:
        if wxsum > 0 and wysum > 0:
            xshift = xsum / wxsum
            yshift = ysum / wysum
            print('Mean x,y shift from reference aperture(s) = {:.2f}, {:.2f}'.format(xshift, yshift))

        else:
            raise hcam.HipercamError('reference aperture fit(s) failed; giving up.')

    else:
        # no reference apertures. All individual
        xshift, yshift = 0., 0.


    # now go over all apertures
    for apnam, aper:
        if aper.ref:
            if store['apnam']['efwhm'] <= 0.:
                # Move failed reference fit to the mean shift
                aper.x += xshift
                aper.y += yshift

        else:

            # extract Windat for this reference aperture
            wind = ccd[ccdwin[apnam]]

            if isinstance(read, hcam.CCD):
                rd = read[ccdwin[apnam]].data
            else:
                rd = read

            if isinstance(gain, hcam.CCD):
                gn = gain[ccdwin[apnam]].data
            else:
                gn = gain

            # get sub-windat around start position
            shbox = int(rfile['apertures']['search_half_width_non'])
            swind = wind.window(aper.x+xshift-shbox, aper.x+xshift+shbox, aper.y+yshift-shbox, aper.y+yshift+shbox)

            # carry out initial search
            x,y,peak = swind.find(float(rfile['apertures']['search_smooth_fwhm']), False)

            # now for a more refined fit. First extract fit Windat
            fhbox = int(rfile['apertures']['fit_half_width'])
            fwind = wind.window(x-fhbox, x+fhbox, y-fhbox, y+fhbox)
            sky = np.percentile(fwind.data, 25)

            # get some parameters from previous run where possible
            if apnam in store and store[apnam]['efwhm'] > 0.:
                fit_fwhm = store[apnam]['fwhm']
            else:
                fit_fwhm = float(rfile['apertures']['fit_fwhm'])

            if apnam in store and store[apnam]['ebeta'] > 0.:
                fit_beta = store[apnam]['beta']
            else:
                fit_beta = float(rfile['apertures']['fit_beta'])

            # refine the Aperture position by fitting the profile
            try:
                (sky, height, x, y, fwhm, beta), \
                    (esky, eheight, ex, ey, efwhm, ebeta), \
                    extras = hcam.combFit(
                    fwind, rfile.method, sky, peak-sky,
                    x, y, fit_fwhm,
                    float(rfile['apertures']['fit_fwhm_min']),
                    rfile['apertures']['fit_fwhm_fixed'] == 'yes',
                    fit_beta, rd, gn
                    )

                if height > float(rfile['apertures']['height_min']):
                    aper.x = x
                    aper.y = y

                    # store some stuff for next time
                    store[apnam] = {'fwhm' : fwhm, 'efwhm' : efwhm,
                                    'beta' : beta, 'ebeta' : ebeta}
                else:
                    print('CCD {:s}, aperture {:s}, peak = {:.1f} < {:s}'.format(
                            cnam, apnam, height, rfile['apertures']['height_min']),
                          file=sys.stderr)
                    store[apnam] = {'efwhm' : -1, 'ebeta' : -1}


            except hcam.HipercamError as err:
                print('CCD {:s}, reference aperture {:s}, fit failed'.format(cnam, apnam), file=sys.stderr)
                store[apnam] = {'efwhm' : -1, 'ebeta' : -1}


