import sys
import os
import math
from collections import OrderedDict as Odict
from multiprocessing import Pool

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

###########################################
#
# makedata -- generates fake multi-CCD data
#
###########################################

def makedata(args=None):
    """Script to generate multi-CCD test data given a set of parameters defined in
    a config file (parsed using configparser). This allows things such as a
    bias frame, flat field variations offsets, scale factor and rotations
    between CCDs, and temporal variations.

    Arguments::

       config   : (string)
          file defining the parameters.

       parallel : (bool)
          True / yes etc to run in parallel. Be warned: it does not always make things faster,
          which I assume is the result of overheads when parallelising. It will simply use whatever
          CPUs are available.

    Depending upon the setting in the config file, this could generate a large
    number of different files and so the first time you run it, you may want
    to do so in a clean directory.

    Config file format: see the documentation of configparser for the general
    format of the config files expected by this routine. Essentially there are
    a series of sections, e.g.:

    [ccd 1]
    nxtot = 2048
    nytot = 1048
    .
    .
    .

    which define all the parameters needed. There are many others to simulate
    a bias offset, a flat field; see the example file ??? for a fully-documented
    version.

    """
    import configparser
    global _gframe, _gfield

    if args is None:
        args = sys.argv[1:]

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'makedata', args) as cl:

        # Register parameters
        cl.register('config', Cline.LOCAL, Cline.PROMPT)
        cl.register('parallel', Cline.LOCAL, Cline.PROMPT)

        config = cl.get_value('config', 'configuration file',
                              cline.Fname('config'))
        parallel = cl.get_value('parallel', 'add targets in parallel?', False)

    # Read the config file
    conf = configparser.ConfigParser()
    conf.read(config)

    # Determine whether files get overwritten or not
    overwrite = conf.getboolean('general', 'overwrite') \
                if 'overwrite' in conf['general'] else False
    dtype = conf['general']['dtype'] \
            if 'dtype' in conf['general'] else None

    # Top-level header
    thead = fits.Header()
    thead.add_history('Created by makedata')

    # Store the CCD labels and parameters and their dimensions. Determine
    # maximum dimensions for later use when adding targets
    ccd_pars = Odict()
    maxnx = 0
    maxny = 0
    for key in conf:
        if key.startswith('ccd'):

            # translate parameters
            nxtot = int(conf[key]['nxtot'])
            nytot = int(conf[key]['nytot'])
            xcen = float(conf[key]['xcen'])
            ycen = float(conf[key]['ycen'])
            angle = float(conf[key]['angle'])
            scale = float(conf[key]['scale'])
            xoff = float(conf[key]['xoff'])
            yoff = float(conf[key]['yoff'])
            fscale = float(conf[key]['fscale'])
            toff = float(conf[key]['toff'])

            field = hcam.Field.rjson(conf[key]['field']) \
                if 'field' in conf[key] else None
            ndiv = int(conf[key]['ndiv']) \
                if field is not None else None

            # determine maximum total dimension
            maxnx = max(maxnx, nxtot)
            maxny = max(maxny, nytot)

            # store parameters
            ccd_pars[key[3:].strip()] = {
                'nxtot' : nxtot,
                'nytot' : nytot,
                'xcen' : xcen,
                'ycen' : ycen,
                'angle' : angle,
                'scale' : scale,
                'xoff' : xoff,
                'yoff' : yoff,
                'fscale' : fscale,
                'toff'  : toff,
                'field'  : field,
                'ndiv'  : ndiv
            }

    if not len(ccd_pars):
        raise ValueError('hipercam.makedata: no CCDs found in ' + config)

    # get the timing data
    utc_start = Time(conf['timing']['utc_start'])
    exposure = float(conf['timing']['exposure'])
    deadtime = float(conf['timing']['deadtime'])

    # Generate the CCDs, store the read / gain values
    ccds = hcam.Group()
    rgs = {}
    for cnam, pars in ccd_pars.items():
        # Generate the Windats
        winds = hcam.Group()
        rgs[cnam] = {}
        for key in conf:
            if key.startswith('window'):
                iccd, wnam = key[6:].split()
                if iccd == cnam:
                    llx = int(conf[key]['llx'])
                    lly = int(conf[key]['lly'])
                    nx = int(conf[key]['nx'])
                    ny = int(conf[key]['ny'])
                    xbin = int(conf[key]['xbin'])
                    ybin = int(conf[key]['ybin'])
                    wind = hcam.Windat(hcam.Window(llx,lly,nx,ny,xbin,ybin))

                    # Accumulate Windats in the CCD
                    winds[wnam] = wind

                    # Store read / gain value
                    rgs[cnam][wnam] = (
                        float(conf[key]['read']), float(conf[key]['gain'])
                    )

        # Generate header with timing data
        head = fits.Header()
        td = TimeDelta(pars['toff'], format='sec')
        utc = utc_start + td
        head['UTC'] = (utc.isot, 'UTC at mid exposure')
        head['MJD'] = (utc.mjd, 'MJD at mid exposure')
        head['EXPOSE'] = (exposure, 'Exposure time, seconds')
        head['TIMEOK'] = (True, 'Time status flag')

        # Accumulate CCDs
        ccds[cnam] = hcam.CCD(winds,pars['nxtot'],pars['nytot'],head)

    # Make the template MCCD
    mccd = hcam.MCCD(ccds, thead)

    # Make a flat field
    flat = mccd.copy()
    if 'flat' in conf:
        rms = float(conf['flat']['rms'])
        for ccd in flat.values():
            # Generate dust
            nspeck = int(conf['flat']['nspeck'])
            if nspeck:
                radius = float(conf['flat']['radius'])
                depth = float(conf['flat']['depth'])
                specks = []
                for n in range(nspeck):
                    x = np.random.uniform(0.5,ccd.nxtot+0.5)
                    y = np.random.uniform(0.5,ccd.nytot+0.5)
                    specks.append(Dust(x,y,radius,depth))

            # Set the flat field values
            for wind in ccd.values():
                wind.data = np.random.normal(1.,rms,(wind.ny,wind.nx))
                if nspeck:
                    wind.add_fxy(specks)

        flat.head['DATATYPE'] = ('Flat field','Artificially generated')
        fname = hcam.add_extension(conf['flat']['flat'],hcam.HCAM)
        flat.wfits(fname, overwrite)
        print('Saved flat field to ',fname)
    else:
        # Set the flat to unity
        flat.set_const(1.0)
        print('No flat field generated')

    # Make a bias frame
    bias = mccd.copy()
    if 'bias' in conf:
        mean = float(conf['bias']['mean'])
        rms = float(conf['bias']['rms'])
        for ccd in bias.values():
            for wind in ccd.values():
                wind.data = np.random.normal(mean,rms,(wind.ny,wind.nx))

        bias.head['DATATYPE'] = ('Bias frame','Artificially generated')
        fname = hcam.add_extension(conf['bias']['bias'],hcam.HCAM)
        bias.wfits(fname, overwrite)
        print('Saved bias frame to ',fname)
    else:
        # Set the bias to zero
        bias.set_const(0.)
        print('No bias frame generated')

    # Everything is set to go, so now generate data files
    nfiles = int(conf['files']['nfiles'])
    if nfiles == 0:
        out = mccd*flat + bias
        fname = hcam.add_extension(conf['files']['root'],hcam.HCAM)
        out.wfits(fname, overwrite)
        print('Written data to',fname)
    else:
        # file naming info
        root = conf['files']['root']
        ndigit = int(conf['files']['ndigit'])

        # movement
        xdrift = float(conf['movement']['xdrift'])
        ydrift = float(conf['movement']['ydrift'])
        nreset = int(conf['movement']['nreset'])
        jitter = float(conf['movement']['jitter'])

        print('Now generating data')

        tdelta = TimeDelta(exposure+deadtime, format='sec')

        for nfile in range(nfiles):

            # copy over template (into a global variable for multiprocessing speed)
            _gframe = mccd.copy()

            # get x,y offset
            xoff = np.random.normal(xdrift*(nfile % nreset), jitter)
            yoff = np.random.normal(ydrift*(nfile % nreset), jitter)

            # create target fields for each CCD
            _gfield = {}
            for cnam in _gframe.keys():
                p = ccd_pars[cnam]
                if p['field'] is not None:
                    # get field modification settings
                    transform = Transform(
                        p['nxtot'], p['nytot'], p['xcen'], p['ycen'],
                        p['angle'], p['scale'],
                        p['xoff']+xoff, p['yoff']+yoff
                    )
                    fscale = p['fscale']
                    _gfield[cnam] = p['field'].modify(transform, fscale)



            # add the targets in (slow step)
            if parallel:
                # run in parallel on whatever cores are available
                args = [(cnam,ccd_pars[cnam]['ndiv']) for cnam in _gfield]
                with Pool() as pool:
                    ccds = pool.map(worker, args)
                for cnam in _gfield:
                    _gframe[cnam] = ccds.pop(0)
            else:
                # single core
                for cnam in _gfield:
                    ccd = _gframe[cnam]
                    ndiv = ccd_pars[cnam]['ndiv']
                    field = _gfield[cnam]
                    for wind in ccd.values():
                        field.add(wind, ndiv)

            # Apply flat
            _gframe *= flat

            # Add noise
            for cnam, ccd in _gframe.items():
                for wnam, wind in ccd.items():
                    readout, gain = rgs[cnam][wnam]
                    wind.add_noise(readout, gain)

            # Apply bias
            _gframe += bias

            # data type on output
            if dtype == 'float32':
                _gframe.float32()
            elif dtype == 'uint16':
                _gframe.uint16()

            # Save
            fname = '{0:s}{1:0{2:d}d}{3:s}'.format(
                root,nfile+1,ndigit,hcam.HCAM
            )
            _gframe.wfits(fname, overwrite)
            print('Written file {0:d} to {1:s}'.format(nfile+1,fname))

            # update times in template
            for ccd in mccd.values():
                head = ccd.head
                utc = Time(head['UTC']) + tdelta
                head['UTC'] = (utc.isot, 'UTC at mid exposure')
                head['MJD'] = (utc.mjd, 'MJD at mid exposure')

#############################################
#
# makefield -- makes an artificial star field
#
#############################################

def makefield(args=None):
    """Script to generate an artificial star field which is saved to a disk
    file, a first step in generating fake data. The resulting files are needed
    by makedata if you want to add in artificial star fields. If the name
    supplied corresponds to an existing file, an attempt will be made to read
    it in first and then add to it. In this way a complex field can be
    generated. The targets are distributed at random, with random peak heights
    based on constant luminosity objects distributed throughout 3D space, and
    random angles uniform over the input range. All targets otherwise have the
    same shape thus multiple calls are needed to generate a field of objects
    of multiple shapes. Ellipsoidal "Moffat" functions [1/(1+r^2)^beta] are
    used.

    Arguments::

       fname : (string)
          file to add to. Will be created if it does not exist. [input, optional / output]

       ntarg : (int)
          The number of targets to add to the field.

       x1 : (float)
          The left-hand limit of the field [unbinned pixels]

       x2 : (float)
          The right-hand limit of the field [unbinned pixels]

       y1 : (float)
          The lowest limit of the field [unbinned pixels]

       y2 : (float)
          The highest limit of the field [unbinned pixels]

       h1 : (float)
          Minimum peak height [counts per unbinned pixel]

       h2 : (float)
          Maximum peak height [counts per unbinned pixel]

       angle1 : (float)
          Lower limit on axis 1, anti-clockwise from X-axis [degrees]

       angle2 : (float)
          Upper limit on axis 1, anti-clockwise from X-axis [degrees]

       fwhm1 : (float)
          FWHM along axis 1 [unbinned pixels]

       fwhm2 : (float)
          FWHM along axis 2 [unbinned pixels]

       beta : (float)
          Moffat function exponent

    """
    if args is None:
        args = sys.argv[1:]

    # create Cline object
    with Cline('HIPERCAM_ENV', '.hipercam', 'makefield', args) as cl:

        # register parameters
        cl.register('fname', Cline.LOCAL, Cline.PROMPT)
        cl.register('ntarg', Cline.LOCAL, Cline.PROMPT)
        cl.register('x1', Cline.LOCAL, Cline.PROMPT)
        cl.register('x2', Cline.LOCAL, Cline.PROMPT)
        cl.register('y1', Cline.LOCAL, Cline.PROMPT)
        cl.register('y2', Cline.LOCAL, Cline.PROMPT)
        cl.register('h1', Cline.LOCAL, Cline.PROMPT)
        cl.register('h2', Cline.LOCAL, Cline.PROMPT)
        cl.register('angle1', Cline.LOCAL, Cline.PROMPT)
        cl.register('angle2', Cline.LOCAL, Cline.PROMPT)
        cl.register('fwhm1', Cline.LOCAL, Cline.PROMPT)
        cl.register('fwhm2', Cline.LOCAL, Cline.PROMPT)
        cl.register('beta', Cline.LOCAL, Cline.PROMPT)
        cl.register('fmin', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        fname = cl.get_value('fname', 'file to save field to',
                             cline.Fname('field', hcam.FIELD, cline.Fname.NEW))
        if os.path.exists(fname):
            # Initialise the field from a file
            field = hcam.Field.rjson(fname)
            print('>> Loaded a field of',len(field),'objects from',fname)
        else:
            # Create an empty field
            field = hcam.Field()
            print('>> Created an empty field.')

        ntarg = cl.get_value('ntarg', 'number of targets', 100, 1)
        x1 = cl.get_value('x1', 'left-hand limit of field', -10.)
        x2 = cl.get_value('x2', 'right-hand limit of field', 2000., x1)
        y1 = cl.get_value('y1', 'lower limit of field', -10.)
        y2 = cl.get_value('y2', 'upper limit of field', 1000., y1)
        h1 = cl.get_value('h1', 'lower peak height limit', 0.1, 1.e-6)
        h2 = cl.get_value('h2', 'upper peak height limit', 1000., h1)
        angle1 = cl.get_value('angle1', 'lower limit of axis 1 angle', 0., -360., 360.)
        angle2 = cl.get_value('angle2', 'angle of major axis', 0., angle1, 360.)
        fwhm1 = cl.get_value('fwhm1', 'FWHM along axis 1', 4., 1.e-6)
        fwhm2 = cl.get_value('fwhm2', 'FWHM along axis 2', 4., 1.e-6)
        beta = cl.get_value('beta', 'Moffat exponent', 4., 1.0)
        fmin = cl.get_value('fmin', 'minimum flux level (counts/pix)',
                            0.001, 0.0)

    # add targets
    field.add_random(ntarg, x1, x2, y1, y2, h1, h2, angle1, angle2,
                     fwhm1, fwhm2, beta, fmin)

    # save result
    field.wjson(fname)

    print('>> Saved a field of',len(field),'objects to',fname)

############################################################################
#
# From this point on come helper methods and classes that are not externally
# visible
#
############################################################################

class Dust:
    """This to generate gaussian dust specks on a flat field using
    the add_fxy method of Windats"""

    def __init__(self, xcen, ycen, rms, depth):
        self.xcen = xcen
        self.ycen = ycen
        self.rms = rms
        self.depth = depth

    def __call__(self, x, y, f):
        ok = (x > self.xcen-5.*self.rms) & (x < self.xcen+5.*self.rms) & \
             (y > self.ycen-5.*self.rms) & (y < self.ycen+5.*self.rms)

        if ok.any():
            rsq = (x[ok]-self.xcen)**2 + (y[ok]-self.ycen)**2
            f[ok] -= self.depth*np.exp(-rsq/(2.*self.rms**2))

class Transform:
    """Field transformation class. This is a callable that can be sent to the
    :class:`Field` method `modify`.
    """

    def __init__(self, nxtot, nytot, xcen, ycen, angle, scale, xoff, yoff):
        self.nxtot = nxtot
        self.nytot = nytot
        self.xcen = xcen
        self.ycen = ycen
        self.angle = angle
        self.scale = scale
        self.xoff = xoff
        self.yoff = yoff

    def __call__(self, x, y):
        # Rotator centre
        xc = (self.nxtot+1)/2. + self.xcen
        yc = (self.nytot+1)/2. + self.ycen

        # Rotate around CCD centre
        cosa = np.cos(np.radians(self.angle))
        sina = np.sin(np.radians(self.angle))
        xcen =  cosa*(x-xc) + sina*(y-yc)
        ycen = -sina*(x-xc) + cosa*(y-yc)

        # Change image scale
        xcen *= self.scale
        ycen *= self.scale

        # Apply offset
        xcen += xc + self.xoff
        ycen += yc + self.yoff

        # Return change in position
        return (xcen-x,ycen-y)

# globals and a function used in the multiprocessing used in 'makedata'. The globals
# are to reduce the need for pickling / unpickling large bits of data
_gframe = None
_gfield = None

def worker(arg):
    global _gframe, _gfield
    cnam, ndiv = arg
    ccd = _gframe[cnam]
    field = _gfield[cnam]
    for wind in ccd.values():
        field.add(wind, ndiv)
    return ccd
