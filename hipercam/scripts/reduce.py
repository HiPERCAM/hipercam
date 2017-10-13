import sys
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
import configparser
from collections import OrderedDict

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
           ULTRACAM / ULTRASPEC, 'h' for HiPERCAM. This is needed because of
           the different formats.

        ldevice  : (string) [hidden]
           Plot device for the light curves. PGPLOT is used so this should be
           a PGPLOT-style name, e.g. '/xs', '1/xs' etc. At the moment only
           ones ending /xs are supported.

        lwidth  : (float) [hidden]
           plot width (inches). Set = 0 to let the program choose.

        lheight : (float) [hidden]
           plot height (inches). Set = 0 to let the program choose. BOTH width
           AND height must be non-zero to have any effect

        rfile   : (string)
           the "reduce" file, i.e. ASCII text file suitable for reading by
           ConfigParser. Best seen by example as it has many parts.

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
           maximum time to wait between attempts to find a new exposure,
           seconds.

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
        source = cl.get_value(
            'source', 'data source [s(erver), l(ocal), f(ile list)]',
            'l', lvals=('s','l','f'))

        if source == 's' or source == 'l':
            inst = cl.get_value(
                'inst', 'instrument [h(ipercam), u(ltracam/spec)]',
                'h', lvals=('h','u'))

        # plot device stuff
        ldevice = cl.get_value('ldevice', 'light curve plot device', '1/xs')
        lwidth = cl.get_value('lwidth', 'light curve plot width (inches)', 0.)
        lheight = cl.get_value('lheight', 'light curve plot height (inches)', 0.)

        # the reduce file
        rfilen = cl.get_value(
            'rfile', 'reduce file', cline.Fname('reduce.red',hcam.RED))
        rfile = Rfile.fromFile(rfilen)

        if source == 's' or source == 'l':
            run = cl.get_value('run', 'run name', 'run005')
            first = cl.get_value('first', 'first frame to plot', 1, 1)

            if source == 's':
                twait = cl.get_value(
                    'twait', 'time to wait for a new frame [secs]', 1., 0.)
                tmax = cl.get_value(
                    'tmax', 'maximum time to wait for a new frame [secs]',
                    10., 0.)

        else:
            # set inst = 'h' as only lists of HiPERCAM files are supported
            inst = 'h'
            run = cl.get_value(
                'flist', 'file list', cline.Fname('files.lis',hcam.LIST))
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

    # define and draw panels. 'theight' is the total height of all panel which
    # will be used to work out the fraction occupied by each one
    theight = 1.0

    # define standard viewport. We set the character height to ensure there is
    # enough room around the edges.
    pgsch(max(hcam.pgp.Params['axis.label.ch'],hcam.pgp.Params['axis.number.ch']))
    pgvstd()
    xv1, xv2, yv1, yv2 = pgqvp()

    # Light curve is always the top panel. 
    yv1_lc = yv1 + (theight-1.)/theight*(yv2-yv1)

    # set up the light curve panel
    lcpanel = Panel(
        lcdev, xv1, xv2, yv1_lc, yv2, 'Time [mins]',
        'Flux' if rfile['light']['linear'] else 'Magnitudes', '',
        'bcnst', 'bcnst', 0, rfile['light']['extend_xrange'],
        rfile['light']['y1'], rfile['light']['y2']
        )

    # plot it
    lcpanel.plot()

    # a couple of initialisations
    total_time = 0 # time waiting for new frame
    fpos = [] # list of target positions to fit

    # dictionary of dictionaries for looking up the
    # window associated with a given aperture, i.e.
    # mccdwins[cnam][apnam] give the name of the Windat.
    mccdwins = {}

    # create buffers to store the points plotted in each panel
    lcbuffer = []
    for tconf in rfile['light']['plot']:
        lcbuffer.append(LC(tconf, rfile['light']['linear']))

    # get frames
    with hcam.data_source(instrument, run, flist, server, first) as spool:

        # 'spool' is an iterable source of MCCDs
        for nframe, mccd in enumerate(spool):

            # None objects are returned from failed server reads. This could
            # be because the file is still exposing, so we hang about.
            if server and mccd is None:

                if tmax < total_time + twait:
                    print(
                        'Have waited for {:.1f} sec. cf tmax = {:.1f}; will wait no more'.format(
                            total_time, tmax)
                    )
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
                print('File {:d}: '.format(n+1))
            else:
                print('Frame {:d}: '.format(mccd.head['NFRAME']))

            if nframe == 0 and rfile['calibration']['crop']:
                # This is the very first data frame read in. We need to trim
                # the calibrations, on the assumption that all data frames
                # will have the same format
                rfile.crop(mccd)

            # container for the results from each CCD
            results = {}
            for cnam in mccd:
                # get the apertures
                ccdaper = rfile.aper[cnam]
                if len(ccdaper) == 0:
                    continue

                # get the CCD and the apertures
                ccd = mccd[cnam]

                if nframe == 0:
                    # first time through, work out an array of which window
                    # each aperture lies in we will assume this is fixed for
                    # the whole run, i.e. that apertures do not drift from one
                    # window to another. Set to None if no window found
                    mccdwins[cnam] = {}
                    for apnam, aper in ccdaper.items():
                        for wnam, wind in ccd.items():
                            if wind.distance(aper.x,aper.y) > 0:
                                mccdwins[cnam][apnam] = wnam
                                break
                        else:
                            mccdwins[cnam][apnam] = None

                    # initialisations
                    store = {}
                    mfwhm, mbeta = -1, -1

                ccdwins = mccdwins[cnam]

                if isinstance(rfile.readout, hcam.MCCD):
                    read = rfile.readout[cnam]
                else:
                    read = rfile.readout

                if isinstance(rfile.gain, hcam.MCCD):
                    gain = rfile.gain[cnam]
                else:
                    gain = rfile.gain

                # at this point 'ccd' contains all the Windats of a CCD,
                # ccdaper all of its apertures ccdwins the label of the Windat
                # relevant for each aperture, rfile contains some control
                # parameters, read contains the readout noise, either a float
                # or a CCD, gain contains the gain, either a float or a CCD.

                mfwhm, mbeta = moveApers(
                    cnam, ccd, ccdaper, ccdwins, rfile,
                    read, gain, mfwhm, mbeta, store
                )

                # extract flux from all apertures of each CCD
                results[cnam] = extractFlux(
                    cnam, ccd, ccdaper, ccdwins, rfile,
                    read, gain, store, mfwhm
                )

            # plot the light curve
            t = nframe / 10 # temporary test
            replot = plotLight(lcpanel, t, results, rfile, lcbuffer)

            # plot some other stuff ...

            # re-plot
            if replot:
                # start buffering, erase, re-plot, end buffering
                pgbbuf()
                pgeras()

                # re-draw the light curve panel
                lcpanel.plot()

                for lc in lcbuffer:
                    # convert the buffered data into float32 ndarrays
                    t = np.array(lc.t, dtype=np.float32)
                    f = np.array(lc.f, dtype=np.float32)
                    fe = np.array(lc.fe, dtype=np.float32)

                    # Plot the error bars
                    pgsci(lc.ecol)
                    pgerry(t, f-fe, f+fe, 0)

                    # Plot the data
                    pgsci(lc.dcol)
                    pgpt(t, f, 17)

                pgebuf()

# END OF MAIN SECTION

#
# From here is support code not exported outside
#

class Rfile(OrderedDict):
    """
    Class to read and interpret reduce files. Similar
    to configparser but a bit freer.
    """

    # Data
    VERSION = '2017-10-05'

    @classmethod
    def fromFile(cls, filename):
        """
        Builds an Rfile from a reduce file
        """

        rfile = cls()
        insection = False
        with open(filename) as fp:

            for line in fp:

                if not line.startswith('#') and line != '' and not line.isspace():

                    # strip trailing comments
                    comm = line.find('#')
                    if comm != -1:
                        line = line[:comm].strip()

                    if line.startswith('['):
                        # new section
                        section = line[1:line.find(']')].strip()
                        sec = rfile[section] = cls()
                        insection = True

                    elif insection:
                        # another entry to current section
                        key = line[:line.find('=')].strip()
                        val = line[line.find('=')+1:].strip()

                        if key in sec:
                            if isinstance(sec[key],list):
                                sec[key].append(val)
                            else:
                                sec[key] = [sec[key], val]
                        else:
                            sec[key] = val

                    else:
                        raise hcam.HipercamError(
                            'found entry line before any section = \n{:s}'.format(line)
                        )

        # process it. this is a matter of checking entries and
        # in some cases converting them to correct or more
        # convenient forms

        #
        # general section
        #

        if rfile['general']['version'] != Rfile.VERSION:
            # check the version
            raise ValueError(
                'Version mismatch: file = {:s}, reduce = {:s}'.format(
                    rfile['general']['version'], Rfile.VERSION)
            )

        #
        # apertures section
        #
        apsec = rfile['apertures']

        rfile.aper = hcam.MccdAper.fromJson(
            hcam.add_extension(apsec['aperfile'],hcam.APER)
            )
        if apsec['fit_method'] == 'moffat':
            rfile.method = 'm'
        elif  apsec['fit_method'] == 'gaussian':
            rfile.method = 'g'
        else:
            raise NotImplementedError('apertures.fit_method = {:s} not recognised'.format(
                    apsec['fit_method']))

        # type conversions
        toBool(rfile,'apertures','fit_fwhm_fixed')

        apsec['search_half_width_ref'] = int(apsec['search_half_width_ref'])
        apsec['search_half_width_non'] = int(apsec['search_half_width_non'])
        apsec['search_smooth_fwhm'] = float(apsec['search_smooth_fwhm'])
        apsec['fit_fwhm'] = float(apsec['fit_fwhm'])
        apsec['fit_fwhm_min'] = float(apsec['fit_fwhm_min'])
        apsec['fit_beta'] = float(apsec['fit_beta'])
        apsec['fit_half_width'] = int(apsec['fit_half_width'])
        apsec['fit_sigma'] = float(apsec['fit_sigma'])
        apsec['fit_height_min'] = float(apsec['fit_sigma'])

        #
        # calibration section
        #
        calsec = rfile['calibration']

        toBool(rfile,'calibration','crop')

        if calsec['bias'] != '':
            rfile.bias == hcam.MCCD.rfits(
                hcam.add_extension(calsec['bias'],hcam.HCAM)
                )
        else:
            rfile.bias = None

        if calsec['dark'] != '':
            rfile.dark == hcam.MCCD.rfits(
                hcam.add_extension(calsec['dark'],hcam.HCAM)
                )
        else:
            rfile.dark = None

        if calsec['flat'] != '':
            rfile.flat == hcam.MCCD.rfits(
                hcam.add_extension(calsec['flat'],hcam.HCAM)
                )
        else:
            rfile.flat = None

        try:
            rfile.readout = float(calsec['readout'])
        except TypeError:
            rfile.readout = hcam.MCCD.rfits(
                hcam.add_extension(calsec['readout'],hcam.HCAM)
                )

        try:
            rfile.gain = float(calsec['gain'])
        except TypeError:
            rfile.gain = hcam.MCCD.rfits(
                hcam.add_extension(calsec['gain'],hcam.HCAM)
                )

        # Extraction section

        # Separate the extraction entries into lists, check and convert some
        # entries
        extsec = rfile['extraction']

        for cnam in extsec:

            extsec[cnam] = lst = extsec[cnam].split()

            if lst[0] != 'variable' and lst[0] != 'fixed':
                raise ValueError(
                    "first entry of extraction lines must either be 'variable' or 'fixed'"
                    )

            if lst[1] != 'normal' and lst[1] != 'optimal':
                raise ValueError(
                    "second entry of extraction lines must either be 'normal' or 'optimal'"
                    )

            # type conversions
            for i in range(2,len(lst)):
                lst[i] = float(lst[i])

        #
        # sky section
        #
        skysec = rfile['sky']

        if skysec['error'] == 'variance':
            if skysec['method'] == 'median':
                raise ValueError(
                    'sky.error == variance requires sky.method == clipped'
                )
        elif skysec['error'] != 'photon':
            raise ValueError(
                "sky.error must be either 'variance' or 'photon'"
            )

        if skysec['method'] != 'clipped' and skysec['method'] != 'median':
            raise ValueError(
                "sky.method must be either 'clipped' or 'median'"
            )

        skysec['thresh'] = float(skysec['thresh'])

        #
        # Light curve plot section
        #

        ligsec = rfile['light']

        targ = ligsec['plot']
        if isinstance(targ, str):
            targ = [targ]

        # convert entries to the right type here and try colours to
        # PGPLOT colour indices.
        for n in range(len(targ)):
            cnam, tnm, cnm, off, fac, dcol, ecol = targ[n].split()
            targ[n] = {
                'ccd' : cnam,
                'targ' : tnm,
                'comp' : cnm,
                'fac' : fac,
                'off' : float(off),
                'fac' : float(fac),
                'dcol' : ctrans(dcol),
                'ecol' : ctrans(ecol)
                }

        ligsec['xrange'] = float(ligsec['xrange'])
        ligsec['extend_xrange'] = float(ligsec['extend_xrange'])
        if ligsec['extend_xrange'] <= 0:
            raise ValueError('light.extend_xrange must be > 0')

        toBool(rfile, 'light', 'linear')
        toBool(rfile, 'light', 'yrange_fixed')
        ligsec['y1'] = float(ligsec['y1'])
        ligsec['y2'] = float(ligsec['y2'])
        ligsec['extend_yrange'] = float(ligsec['extend_yrange'])
        if ligsec['extend_yrange'] <= 0:
            raise ValueError('light.extend_yrange must be > 0')

        # We are finally done reading and hecking the reduce script.
        # rfile[section][param] should from now on return something
        # useful and somewhat reliable.

        return rfile

    def crop(self, mccd):
        """
        This uses a template file 'mccd' to try to get the calibration files into the same format.
        """
        if self.bias is not None:
            self.bias = self.bias.crop(mccd)

        if self.dark is not None:
            self.dark = self.dark.crop(mccd)

        if self.flat is not None:
            self.flat = self.flat.crop(mccd)

        if isinstance(self.readout, hcam.MCCD):
            self.readout = self.readout.crop(mccd)

        if isinstance(self.gain, hcam.MCCD):
            self.gain = self.gain.crop(mccd)


def moveApers(cnam, ccd, ccdaper, ccdwin, rfile, read, gain, mfwhm, mbeta, store):
    """Encapsulates aperture re-positioning and resizing. 'store' is a
    dictionary of results that will be used to start the fits from one frame
    to the next. It must start as {}.

    It operates by first shifting any reference apertures, then non-linked
    apertures, and finally linked apertures.

    Finally it applies a resizing strategy deduced from the 'extraction'
    section of rfile.

    It returns a weighted mean FWHM suitable for re-sizing the apertures. This
    is only meaningful if the profile is fitted. It is returned as -1 if not.

    Arguments::

       cnam      : (string)
           CCD label

       ccd       : (CCD)
           the CCD

       ccdaper   : (CcdAper)
           the Apertures

       ccdwin    : 
           the Window label corresponding to each Aperture

       rfile     : (Rfile)
           reduce file configuration parameters

       read      : (float | CCD)
           readout noise

       gain      : (float | CCD)
           readout noise

       mfwhm     : (float)
           mean FWHM used as a starter for fits. Start at -1 and the reduce
           file starter value will be used.

       mbeta     : (float)
           mean beat used as a starter for future fits. Start at -1 and the reduce
           file starter value will be used.


    Returns: (mfwhm, mbeta) for use next time around; these are set to -1 if
    things go wrong. When OK, they are also stored internally in
    rfile['apertures']['fit_fwhm'] and rfile['apertures']['fit_beta'] for
    retrieval.

    """

    # short-hand that will used a lot
    apsec = rfile['apertures']

    # first of all try to get a mean shift from the reference apertures.  we
    # move any of these apertures that are fitted OK
    xsum, ysum = 0., 0.
    wxsum, wysum = 0., 0.
    ref = False

    # next used to work out weighted mean FWHM and beta values
    fsum, wfsum = 0, 0
    bsum, wbsum = 0, 0
    for apnam, aper in ccdaper.items():
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
            shbox = apsec['search_half_width_ref']
            swind = wind.window(aper.x-shbox, aper.x+shbox, aper.y-shbox, aper.y+shbox)

            # carry out initial search
            x,y,peak = swind.find(apsec['search_smooth_fwhm'], False)

            # now for a more refined fit. First extract fit Windat
            fhbox = apsec['fit_half_width']
            fwind = wind.window(x-fhbox, x+fhbox, y-fhbox, y+fhbox)

            sky = np.percentile(fwind.data, 25)

            # get some parameters from previous run where possible
            if mfwhm > 0.:
                fit_fwhm = mfwhm
            else:
                fit_fwhm = apsec['fit_fwhm']

            if mbeta > 0.:
                fit_beta = mbeta
            else:
                fit_beta = apsec['fit_beta']

            # refine the Aperture position by fitting the profile
            try:
                (sky, height, x, y, fwhm, beta), \
                    (esky, eheight, ex, ey, efwhm, ebeta), extras = \
                    hcam.combFit(
                    fwind, rfile.method, sky, peak-sky,
                    x, y, fit_fwhm, apsec['fit_fwhm_min'],
                    apsec['fit_fwhm_fixed'], fit_beta, rd, gn
                    )

                if height > apsec['height_min']:
                    dx = x - aper.x
                    wx = 1./ex**2
                    wxsum += wx
                    xsum += wx*dx

                    dy = y - aper.y
                    wy = 1./ey**2
                    wysum += wy
                    ysum += wy*dy

                    # store some stuff for next time
                    store[apnam] = {
                        'ex' : ex, 'ey' : ey,
                        'fwhm' : fwhm, 'efwhm' : efwhm,
                        'beta' : beta, 'ebeta' : ebeta,
                        'dx' : x-aper.x, 'dy' : y-aper.y
                    }

                    # update position
                    aper.x = x
                    aper.y = y

                    if efwhm > 0.:
                        # average FWHM computation
                        wf = 1./efwhm**2
                        fsum += wf*fwhm
                        wfsum += wf

                    if ebeta > 0.:
                        # average beta computation
                        wb = 1./ebeta**2
                        bsum += wb*beta
                        wbsum += wb

                else:
                    print('CCD {:s}, reference aperture {:s}, peak = {:.1f} < {:s}'.format(
                            cnam, apnam, height, apsec['height_min']),
                          file=sys.stderr)

                    store[apnam] = {
                        'ex' : -1, 'ey' : -1,
                        'fwhm' : 0, 'efwhm' : -1,
                        'beta' : 0, 'ebeta' : -1,
                        'dx' : 0, 'dy' : 0
                    }

            except hcam.HipercamError as err:
                print('CCD {:s}, reference aperture {:s}, fit failed'.format(cnam, apnam), file=sys.stderr)

                store[apnam] = {
                    'ex' : -1, 'ey' : -1,
                    'fwhm' : 0, 'efwhm' : -1,
                    'beta' : 0, 'ebeta' : -1,
                    'dx' : 0, 'dy' : 0
                }

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

    # now go over all apertures except linked ones. reference
    # apertures are skipped except any that failed are shifted
    # by the mean shift.
    for apnam, aper in ccdaper.items():
        if aper.ref:
            if store['apnam']['efwhm'] <= 0.:
                # Move failed reference fit to the mean shift
                aper.x += xshift
                aper.y += yshift
                store['apnam']['dx'] = xshift
                store['apnam']['dy'] = yshift

        elif not aper.is_linked():

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
            shbox = apsec['search_half_width_non']
            swind = wind.window(
                aper.x+xshift-shbox, aper.x+xshift+shbox,
                aper.y+yshift-shbox, aper.y+yshift+shbox
            )

            # carry out initial search
            x,y,peak = swind.find(apsec['search_smooth_fwhm'], False)

            # now for a more refined fit. First extract fit Windat
            fhbox = apsec['fit_half_width']
            fwind = wind.window(x-fhbox, x+fhbox, y-fhbox, y+fhbox)

            sky = np.percentile(fwind.data, 25)

            # get some parameters from previous run where possible
            if apnam in store and store[apnam]['efwhm'] > 0.:
                fit_fwhm = store[apnam]['fwhm']
            else:
                fit_fwhm = apsec['fit_fwhm']

            if apnam in store and store[apnam]['ebeta'] > 0.:
                fit_beta = store[apnam]['beta']
            else:
                fit_beta = apsec['fit_beta']

            # refine the Aperture position by fitting the profile
            try:
                (sky, height, x, y, fwhm, beta), \
                    (esky, eheight, ex, ey, efwhm, ebeta), extras = \
                    hcam.combFit(
                    fwind, rfile.method, sky, peak-sky,
                    x, y, fit_fwhm, apsec['fit_fwhm_min'],
                    apsec['fit_fwhm_fixed'], fit_beta, rd, gn
                    )

                if height > apsec['fit_height_min']:
                    # store some stuff for next time and for passing onto
                    # next routine
                    store[apnam] = {
                        'ex' : ex, 'ey' : ey,
                        'fwhm' : fwhm, 'efwhm' : efwhm,
                        'beta' : beta, 'ebeta' : ebeta,
                        'dx' : x-aper.x, 'dy' : y-aper.y
                    }
                    aper.x = x
                    aper.y = y

                    if efwhm > 0.:
                        # average FWHM computation
                        wf = 1./efwhm**2
                        fsum += wf*fwhm
                        wfsum += wf

                    if ebeta > 0.:
                        # average beta computation
                        wb = 1./ebeta**2
                        bsum += wb*beta
                        wbsum += wb

                else:
                    print('CCD {:s}, aperture {:s}, peak = {:.1f} < {:s}'.format(
                            cnam, apnam, height, apsec['fit_height_min']),
                          file=sys.stderr)
                    aper.x += xshift
                    aper.y += yshift

                    store[apnam] = {
                        'ex' : -1, 'ey' : -1,
                        'fwhm' : 0, 'efwhm' : -1,
                        'beta' : 0, 'ebeta' : -1,
                        'dx' : xshift, 'dy' : yshift
                    }

            except hcam.HipercamError as err:
                print(
                    'CCD {:s}, aperture {:s}, fit failed'.format(
                        cnam, apnam), file=sys.stderr
                )
                aper.x += xshift
                aper.y += yshift

                store[apnam] = {
                    'ex' : -1, 'ey' : -1,
                    'fwhm' : 0, 'efwhm' : -1,
                    'beta' : 0, 'ebeta' : -1,
                    'dx' : xshift, 'dy' : yshift
                }

    # finally the linked ones
    for apnam, aper in ccdaper.items():

        if aper.is_linked():
            aper.x += store[aper.link]['dx']
            aper.y += store[aper.link]['dy']

            store[apnam] = {
                'ex' : -1, 'ey' : -1,
                'fwhm' : 0, 'efwhm' : -1,
                'beta' : 0, 'ebeta' : -1,
                'dx' : store[aper.link]['dx'],
                'dy' : store[aper.link]['dy']
            }

    # update mfwhm and mbeta if we can. If not, set them to -1 as a flag down
    # the line that there is no measured value of these parameters. If things
    # work, store them in the config file for retrieval if mfwhm and mbeta go
    # wrong later.
    if wfsum > 0.:
        mfwhm = fsum / wfsum
        apsec['fit_fwhm'] = mfwhm
    else:
        mfwhm = -1

    if wbsum > 0.:
        mbeta = bsum / wbsum
        apsec['fit_beta'] = mbeta
    else:
        mbeta = -1

    return (mfwhm, mbeta)

def extractFlux(cnam, ccd, ccdaper, ccdwin, rfile, read, gain, store, mfwhm):
    """This extracts the flux of all apertures of a given CCD.

    The steps are (1) aperture resizing, (2) sky background estimation, (3)
    flux extraction.

    It returns the results as a dictionary keyed on the aperture label. Each
    entry returns a list:

    [x, ex, y, ey, fwhm, efwhm, beta, ebeta, counts, ecounts, sky, esky, nsky, nrej, flag]

    flag = bitmask. If flag = 0, all is OK.

       bit 1  : no FWHM fitted for variable extraction
       bit 2  : no sky pixels
       bit 3  : sky aperture goes off edge of window
       bit 4  : target aperture goes off edge of window

    This code::

       bset = flag & (1 << n)

    determines whether bit 'n' is set or not.
    """

    # get the control parameters
    resize, extype, r1fac, r1min, r1max, r2fac, r2min, r2max, \
        r3fac, r3min, r3max = rfile['extraction'][cnam]

    results = {}

    if resize == 'variable':
        if mfwhm <= 0:
            # return early here as there is nothing much we can do.
            print('** CCD {:s}: no measured FWHM to re-size apertures; no extraction of any aperture')
            flag = 1 << 0
            for apnam, aper in ccdaper:
                info = store[apnam]
                results[apnam] = \
                    {
                    'x' : aper.x, 'ex' : info['ex'],
                    'y' : aper.y, 'ey' : info['ey'],
                    'fwhm' : info['fwhm'], 'efwhm' : info['efwhm'],
                    'beta' : info['beta'], 'ebeta' : info['ebeta'],
                    'counts' : 0., 'ecounts' : -1, 
                    'sky' : 0., 'esky' : 0., 'nsky' : 0, 'nrej' : 0,
                    'flag' : flag
                    }
            return results

        else:
            # Re-size the apertures
            for aper in ccdaper.values():
                aper.rtarg = min(r1min, max(r1max, r1fac*mfwhm))
                aper.rsky1 = min(r2min, max(r2max, r2fac*mfwhm))
                aper.rsky2 = min(r3min, max(r3max, r3fac*mfwhm))

    elif resize == 'fixed':
        # do nothing
        pass
    else:
        raise ValueError("CCD {:s}: 'variable' and 'fixed' are the only aperture resizing options".format(
            cnam))

    # apertures are now positioned and re-sized. Finally extract something.
    for apnam, aper in ccdaper.items():

        # initialise flag
        flag = 0

        wnam = ccdwin[apnam]
        wind = ccd[wnam]

        # extract sub-window that includes all of the pixels that could
        # conceivably affect the aperture
        x1,x2,y1,y2 = aper.x-aper.rsky2-wind.xbin, aper.x+aper.rsky2+wind.xbin, \
                      aper.y-aper.rsky2-wind.ybin, aper.y+aper.rsky2+wind.ybin

        swind = wind.window(x1,x2,y1,y2)

        xlo,xhi,ylo,yhi = swind.extent()
        if xlo > aper.x-aper.rsky2 or xhi < aper.x+aper.rsky2 or \
           ylo > aper.y-aper.rsky2 or yhi < aper.y+aper.rsky2:
            # the sky aperture overlaps the edge of the window, set bit 3
            flag |= (1 << 3)

        if xlo > aper.x-aper.rtarg or xhi < aper.x+aper.rtarg or \
           ylo > aper.y-aper.rtarg or yhi < aper.y+aper.rtarg:
            # the target aperture overlaps the edge of the window, set bit 4
            flag |= (1 << 4)

        if isinstance(read, hcam.CCD):
            sread = read[wnam].window(x1,x2,y1,y2).data

        if isinstance(gain, hcam.CCD):
            sgain = gain[wnam].window(x1,x2,y1,y2).data

        # compute X, Y arrays over the sub-window relative to the centre
        # of the aperture and squared to save a little effort
        x = swind.x(np.arange(swind.nx))-aper.x
        y = swind.y(np.arange(swind.ny))-aper.y
        X, Y = np.meshgrid(x, y)
        Rsq = X**2 + Y**2

        # squared aperture radii for comparison
        R1sq, R2sq, R3sq = aper.rtarg**2, aper.rsky1**2, aper.rsky2**2

        # sky selection, accounting for masks
        sok = (Rsq > R2sq) & (Rsq < R3sq)
        for xoff, yoff, radius in aper.mask:
            sok &= (X-xoff)**2 + (Y-yoff)**2 < radius**2

        # generate sky
        dsky = swind.data[sok]
        if len(dsky):

            # we have some sky

            if rfile['sky']['method'] == 'clipped':

                # clipped mean. Take average, compute RMS,
                # reject pixels > thresh*rms from the mean.
                # repeat until no new pixels are rejected.

                thresh = rfile['sky']['thresh']
                ok = np.ones_like(dsky, dtype=bool)
                nrej = 1
                while nrej:
                    slevel = dsky[ok].mean()
                    srms = dsky[ok].std()
                    nold = len(dsky[ok])
                    ok = ok & (np.abs(dsky-slevel) < thresh*srms)
                    nrej = nold - len(dsky[ok])

                nsky = len(dsky[ok])

                # serror -- error in the sky estimate.
                serror = srms/np.sqrt(nsky)

            else:

                # 'median' goes with 'photon'
                slevel = dsky.median()
                nsky = len(dsky)
                nrej = 0

                # srms will be used to substitute use read / gain parameters
                if isinstance(read, hcam.CCD):
                    rd = sread[sok][ok]
                else:
                    rd = read

                if isinstance(gain, hcam.CCD):
                    gn = sgain[sok][ok]
                else:
                    gn = gain
                serror = np.sqrt((rd**2 + np.max(0, dsky[ok])/gn).sum()/nsky**2)

        else:
            # no sky. still return the flux in the aperture but set bit 2 in flag
            flag |= (1 << 2)
            slevel = 0
            serror = -1
            nsky = 0
            nrej = 0

        # target selection
        dok = Rsq < R1sq
        dtarg = swind.data[dok]

        if nsky and rfile['sky']['error'] == 'variance':
            rd = srms
        elif isinstance(read, hcam.CCD):
            rd = sread[dok]
        else:
            rd = read

        if isinstance(gain, hcam.CCD):
            gn = sgain[dok]
        else:
            gn = gain

        if extype == 'normal':
            counts = (dtarg-slevel).sum()
            var = (rd**2 + np.maximum(0, dtarg)/gn).sum()
            if serror > 0:
                # add in factor due to uncertainty in sky estimate
                var += (len(dtarg)*serror)**2
            ecounts = np.sqrt(var)

        elif extype == 'optimal':
            raise NotImplementedError('yet to add optimal extraction; sorry. feel free to whinge')

        info = store[apnam]

        results[apnam] = \
        {
            'x' : aper.x, 'ex' : info['ex'],
            'y' : aper.y, 'ey' : info['ey'],
            'fwhm' : info['fwhm'], 'efwhm' : info['efwhm'],
            'beta' : info['beta'], 'ebeta' : info['ebeta'],
            'counts' : counts, 'ecounts' : ecounts, 
            'sky' : slevel, 'esky' : serror, 'nsky' : nsky, 'nrej' : nrej,
            'flag' : flag
            }

    # finally, we are done
    return results

class Panel:
    """
    Keeps track of the configuration of particular panels of plots so that
    they can be easily re-plotted and selected for additional plotting if
    need be.
    """
    def __init__(self, device, xv1, xv2, yv1, yv2,
                 xlabel, ylabel, tlabel, xopt, yopt,
                 x1, x2, y1, y2):
        """
        This takes all the arguments needs to set up some axes
        at an arbitrary location, using PGPLOT commands pgsvp, pgswin,
        pgbox and pglab. 'device' is the hipercam.pgp.Device to use
        for the plot. It only stores these values. 'plot' actually
        draws the axes.
        """
        self.device = device
        self.xv1 = xv1
        self.xv2 = xv2
        self.yv1 = yv1
        self.yv2 = yv2
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.tlabel = tlabel
        self.xopt = xopt
        self.yopt = yopt
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2

        # to indicate whether this has been used yet
        self.used = False

    def set_lims(self, x1, x2, y1, y2):
        """
        Re-sets the physical limits. Needs a call to plot to actually
        perform the graphical update
        """
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2

    def plot(self):
        """
        Plot the Panel with physical scales x1 to x2, y1 to y2, using the
        properties set on instantiation. x1, x2, y1, y2 are stored for future
        use when re-selecting the panel. After running this you can plot to
        the Panel. If you plot to another Panel, you can return to this one
        using 'select'.
        """

        # select the device
        self.device.select()

        # draw the axes
        pgsci(hcam.pgp.Params['axis.ci'])
        pgsch(hcam.pgp.Params['axis.number.ch'])
        pgsvp(self.xv1, self.xv2, self.yv1, self.yv2)
        pgswin(self.x1,self.x2,self.y1,self.y2)
        pgbox(self.xopt, 0, 0, self.yopt, 0, 0)

        # plot the labels
        pgsci(hcam.pgp.Params['axis.label.ci'])
        pgsch(hcam.pgp.Params['axis.label.ch'])
        pglab(self.xlabel, self.ylabel, self.tlabel)

        self.used = True

    def select(self):
        """
        Selects this panel as the one to plot to. You can only use this if you
        have plotted the panel.
        """
        if not self.used:
            raise hcam.HipercamError('You must plot a panel before slecting it')

        # select the device
        self.device.select()

        # set the physical scales and the viewport
        pgsvp(self.xv1, self.xv2, self.yv1, self.yv2)
        pgswin(self.x1,self.x2,self.y1,self.y2)


def plotLight(lcpanel, t, results, rfile, lcbuffer):
    """Plots one set of results in the light curve panel. Handles storage of
    points for future re-plots, computing new plot limits where needed.

    It returns a bool which if True means that there is a need to re-plot.
    This is deferred since there may be other panels to be adjusted as well
    since if a plot is cleared, everything has to be re-built.
    """

    # by default, don't re-plot
    replot = False

    # shorthand
    lsect = rfile['light']

    # select the light curve panel
    lcpanel.select()

    # Get the current y-range
    ymin, ymax = lcpanel.y1, lcpanel.y2
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    # add points to the plot and buffers, tracking the minimum and maximum values
    tmax = fmin = fmax = None
    for lc in lcbuffer:
        f = lc.add_point(t, results)
        if f is not None:
            fmin = f if fmin is None else min(fmin, f)
            fmax = f if fmax is None else max(fmax, f)
            tmax = t if tmax is None else max(tmax, t)

    # now determine whether we need to re-plot, and the new plot limits.
    # first set the default
    if tmax > lcpanel.x2:
        replot = True
        x1, x2 = lcpanel.x1, lcpanel.x2
        while tmax > x2:
            x2 += lsect['extend_xrange']
    else:
        x1, x2 = lcpanel.x1, lcpanel.x2 

    if lsect['yrange_fixed']:
        y1, y2 = lcpanel.y1, lcpanel.y2

    else:

        if fmin < ymin or fmax > ymax:
            # we are going to have to replot because we have moved
            # outside the y-limits of the panel. We extend a little bit
            # more than necessary according to extend_yrange in order to
            # reduce the amount of such re-plotting
            replot = True
            extend = lsect['extend_yrange']*(ymax-ymin)
            if fmin < ymin:
                ymin = fmin - extend
            if fmax > ymax:
                ymax = fmax + extend

            if lsect['linear']:
                y1, y2 = ymin, ymax
            else:
                y1, y2 = ymax, ymin

        else:
            y1, y2 = lcpanel.y1, lcpanel.y2

    if replot:
        lcpanel.set_lims(x1, x2, y1, y2)

    return replot

class LC:
    """
    Container for light curves so they can be re-plotted as they come in.
    There should be one of these per plot line in the 'light' section.
    """
    def __init__(self, tconf, linear):
        self.cnam = tconf['ccd']
        self.targ = tconf['targ']
        self.comp = tconf['comp']
        self.off = tconf['off']
        self.fac = tconf['fac']
        self.dcol = tconf['dcol']
        self.ecol = tconf['ecol']
        self.linear = linear
        self.t  = []
        self.f  = []
        self.fe = []

    def add_point(self, t, results):
        """
        Extracts the data to be plotted on the light curve
        plot for the LC given the time and the results (for
        all CCDs, as returned by extractFlux. Assuming all is
        OK (errors > 0 for both comparison and target), it stores
        (t, f, fe) in a recarray for possible re-plotting, plots
        the point and returns the value of 'f' plotted to help with
        re-scaling or None if nothing was plotted.
        't' is the time in minutes since the start of the run
        """

        res = results[self.cnam]

        targ = res[self.targ]
        ft = targ['counts']
        fte = targ['ecounts']

        if fte > 0:

            if self.comp != '!':
                comp = res[self.comp]
                fc = comp['counts']
                fce = comp['ecounts']

                if fc > 0:
                    if fce > 0.:
                        f = ft / fc
                        fe = np.sqrt((fte/fc)**2+(t*fce/fc**2)**2)
                    else:
                        return None
                else:
                    return None
            else:
                f = ft
                fe = fte
        else:
            return None

        if not self.linear:
            if f <= 0.:
                return None

            fe = 2.5/np.log(10)*(fe/f)
            f = -2.5*np.log10(f)

        # apply scaling factor and offset
        f *= self.fac
        fe *= self.fac
        f += self.off

        # OK, we are done. Store new point
        self.t.append(t)
        self.f.append(f)
        self.fe.append(fe)

        # Plot the point in minutes from start point
        pgsci(self.ecol)
        pgmove(t, f-fe)
        pgdraw(t, f+fe)
        pgsci(self.dcol)
        pgpt1(t,f,17)

        # return f up the line
        return f

def ctrans(cname):
    """
    Translates a colour name (cname) into a PGPLOT index
    which is the return value. Defaults to index 1 and prints
    a message if cname not recognised.
    """

    if cname == 'red':
        cindex = 2
    elif cname == 'green':
        cindex = 3
    elif cname == 'blue':
        cindex = 4
    else:
        print('Failed to recognize colour = {:s}; defaulting to black'.format(
            cname))
    return cindex


def toBool(rfile, section, param):
    """
    Converts yes / no responses into True / False. This is used a few times
    in the code to read a reduce file.

    Arguments::

       rfile  : (Rfile)
         the reduce file, an Odict of Odicts

      section : (str)
         the section name

      param   : (str)
         the parameter

    Returns nothing; rfile modified on exit. A ValueError is
    raised if the initial value is neither 'yes' nor 'no.
    """

    if rfile[section][param] == 'yes':
        rfile[section][param] = True
    elif rfile[section][param] == 'no':
        rfile[section][param] = False
    else:
        raise ValueError(
            "{:s}.{:s}: 'yes' or 'no' are the only supported values".format(
                section,param)
            )
