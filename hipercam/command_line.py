"""
command_line contains command line scripts that bundle the various routines together
to carry out various operations such as displaying an image. They are a useful resource
for developing further code.
"""
 
import sys
import os
import math
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from trm.pgplot import *
import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

###########################################################
#
# carith -- arithematic with multi-CCD images
#
###########################################################

def carith(args=None):
    """
    Carries out operations of the form output = input [op] constant
    where [op] is '+', '-', '*' or '/'. This is meant to be used as
    an entry point script. If called within a program, then the first
    argument should be 'cadd', 'csub', 'cmul' or 'cdiv' to define the
    default parameter file name.

    """

    if args is None:
        args = sys.argv
    command = os.path.split(args.pop(0))[1]

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('input', Cline.LOCAL, Cline.PROMPT)
        cl.register('const', Cline.LOCAL, Cline.PROMPT)
        cl.register('output', Cline.LOCAL, Cline.PROMPT)

        prompts = {'cadd' : 'add', 'csub' : 'subtract',
                   'cmul' : 'multiply by', 'cdiv' : 'divide by'}

        infile = cl.get_value('input', 'input file',
                              cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(infile)

        constant = cl.get_value('const', 'constant to ' + prompts[command], 0.)

        outfile = cl.get_value('output', 'output file',
                               cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW))

    # carry out operation
    if command == 'cadd':
        mccd += constant
    elif command == 'csub':
        mccd -= constant
    elif command == 'cmul':
        mccd *= constant
    elif command == 'cdiv':
        mccd /= constant

    # Add a history line
    mccd.head.add_history(
        '{:s} {:s} {:f} {:s}'.format(command,infile,constant,outfile)
    )

    # save result
    mccd.wfits(outfile, True)

###############################################################
#
# combine -- combines images
#
###############################################################

def combine(args=None):
    """Combines a series of images defined by a list using median
    combination.

    Arguments::

        list   : (string)
           list of hcm files with images to combine. The formats of the
           images should all match

        output : (string)
           output file
    """

    if args is None:
        args = sys.argv[1:]

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'combine', args) as cl:

        # register parameters
        cl.register('list', Cline.GLOBAL, Cline.PROMPT)
        cl.register('output', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        flist = cl.get_value(
            'list', 'list of files to combine',
            cline.Fname('files', hcam.LIST)
        )

        outfile = cl.get_value(
            'output', 'output file',
            cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW)
        )

    # Read files into memory
    print('Loading files from {:s}'.format(flist))
    mccds = []
    with hcam.HcamListSpool(flist) as spool:
        for mccd in spool:
            mccds.append(mccd)
    print('Loaded {:d} files'.format(len(mccds)))

    # Combine
    ofile = mccds[0].copy()
    for cnam, occd in ofile.items():
        print('Combining CCD {:s} frames'.format(cnam))
        for wnam, owin in occd.items():
            arrs = []
            for mccd in mccds:
                arrs.append(mccd[cnam][wnam].data)
            arr3d = np.stack(arrs)
            owin.data = np.median(arr3d,axis=0)

    # write out
    ofile.wfits(outfile)
    print('Written result to {:s}'.format(outfile))

###########################################################
#
# grab -- downloads a series of images from a raw data file
#
###########################################################

def grab(args=None):
    """This downloads a sequence of images from a raw data file and writes them
    out to a series CCD / MCCD files. Arguments::

      source : (string)
         's' = server, 'l' = local

      inst   : (string) [hidden]
         Instrument involved. Two choices, 'u' for ULTRCAM or ULTRASPEC, 'h'
         for HiPERCAM. This is needed because of the different formats.

      run    : (string)
         run name to access

      ndigit : (int)
         Files created will be written as 'run005_0013.fits' etc. `ndigit` is
         the number of digits used for the frame number (4 in this case)

      first  : (int)
         First frame to access

      last   : (int)
         Last frame to access, 0 for the lot

      bias   : (string)
         Name of bias frame to subtract, 'none' to ignore.

      dtype  : (string)
         Data type on output. Options::

            'r' : "raw", do nothing; this is safe but could
                  be profligate as you could end up storing 64-bit
                  double precision values if you subtract a bias.

            'f' : Convert to 32-bit floats if already in 64-bit form. Leave
                  integers untouched. This implies some loss of precision,
                  but requires half the space.

            'i' : Convert to 16-bit unsigned integers. A warning will be
                  issued if loss of precision occurs; an exception will
                  be raised if the data are outside the range 0 to 65535.
    """

    if args is None:
        args = sys.argv[1:]

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'grab', args) as cl:

        # register parameters
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('inst', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ndigit', Cline.LOCAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('last', Cline.LOCAL, Cline.PROMPT)
        cl.register('bias', Cline.GLOBAL, Cline.PROMPT)
        cl.register('dtype', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value(
            'source', 'data source [s(erver), l(ocal)]',
            'l', lvals=('s','l')
        )

        inst = cl.get_value(
            'inst', 'instrument used [h(ipercam), u(ltracam/spec)]',
            'h', lvals=('h','u')
        )

        run = cl.get_value('run', 'run name', 'run005')

        ndigit = cl.get_value(
            'ndigit', 'number of digits in frame identifier', 3, 0)

        first = cl.get_value('first', 'first frame to grab', 1, 1)
        last = cl.get_value('last', 'last frame to grab', 0)
        if last < first and last != 0:
            sys.stderr.write('last must be >= first or 0')
            sys.exit(1)

        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        if bias is not None:
            # read the bias frame
            bframe = hcam.MCCD.rfits(bias)

        dtype = cl.get_value(
            'dtype', 'data type [r(aw), f(32-bit float), i(16-bit uint)]',
            'f', lvals=('r','f','i')
        )

    # Now the actual work. First set up the arguments for data_source
    server = source == 's'
    flist = None
    if inst == 'u':
        instrument = 'ULTRA'
    elif inst == 'h':
        instrument = 'HIPER'

    # strip off extensions
    if run.endswith(hcam.HRAW):
        run = run[:run.find(hcam.HRAW)]

    # Finally, we can go
    with hcam.data_source(instrument, run, flist, server, first) as spool:
        nframe = first
        root = os.path.basename(run)
        for frame in spool:

            # subtract bias
            if bias is not None:
                frame -= bframe

            if dtype == 'i':
                frame.uint16()
            elif dtype == 'f':
                frame.float32()

            # write to disk
            fname = '{:s}_{:0{:d}}{:s}'.format(run,nframe,ndigit,hcam.HCAM)
            frame.wfits(fname,True)

            print('Written frame {:d} to {:s}'.format(nframe,fname))
            nframe += 1
            if last and nframe > last: break

###########################################################
#
# hplot -- plots a multi-CCD image.
#
###########################################################

def hplot(args=None):
    """
    Plots a multi-CCD image. Can use PGPLOT or matplotlib. The matplotlib
    version is slightly clunky in its choice of the viewing area but has some
    features that could be useful, in particular, the interactive plot
    (device='/mpl') allows one to pan and zoom and to compare the same part of
    multiple CCDs easily.

    Arguments::

      input  : (string)
         name of MCCD file

      device : (string) [hidden]
         Plot device name. Uses characters after a final trailing '/' to
         identify the type in PGPLOT style. Thus::

                   /xs : PGPLOT xserver interactive plot
                  1/xs : PGPLOT xserver interactive plot called '1'
           plot.ps/cps : PGPLOT colour postscript called 'plot.ps'
           plot.ps/vps : PGPLOT B&W portrait oriented plot
                  /mpl : matplotlib interactive plot
          plot.pdf/mpl : matplotlib PDF plot

      ccd    : (string)
         CCD(s) to plot, '0' for all. If not '0' then '1', '2' or even '3 4'
         are possible inputs (without the quotes). '3 4' will plot CCD '3' and
         CCD '4'. If you want to plot more than one CCD, then you will be
         prompted for the number of panels in the X direction. This parameter
         will not be prompted if there is only one CCD in the file.

      nx     : (int)
         number of panels across to display, prompted if more than one CCD is
         to be plotted.

      iset   : (string) [single character]
         determines how the intensities are determined. There are three
         options: 'a' for automatic simply scales from the minimum to the
         maximum value found on a per CCD basis. 'd' for direct just takes two
         numbers from the user. 'p' for percentile dtermines levels based upon
         percentiles determined from the entire CCD on a per CCD bais.

      ilo    : (float) [if iset=='d']
         lower intensity level

      ihi    : (float) [if iset=='d']
         upper intensity level

      plo    : (float) [if iset=='p']
         lower percentile level

      phi    : (float) [if iset=='p']
         upper percentile level

      xlo    : (float)
         left X-limit, for PGPLOT and matplotlib harcopy plots only as the
         matplotlib interactive plot is zoomable

      xhi    : (float)
         right X-limit

      ylo    : (float)
         bottom Y-limit

      yhi    : (float)
         top Y-limit

      width  : (float) [hidden]
         plot width (inches). Set = 0 to let the program choose.

      height : (float) [hidden]
         plot height (inches). Set = 0 to let the program choose. BOTH width
         AND height must be non-zero to have any effect
    """

    if args is None:
        args = sys.argv[1:]

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', 'hplot', args) as cl:

        # register parameters
        cl.register('input', Cline.LOCAL, Cline.PROMPT)
        cl.register('device', Cline.LOCAL, Cline.HIDE)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('nx', Cline.LOCAL, Cline.PROMPT)
        cl.register('iset', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ilo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ihi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('plo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('phi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('xlo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('xhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ylo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('yhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)

        # get inputs
        frame = cl.get_value('input', 'frame to plot',
                             cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(frame)

        device = cl.get_value('device', 'plot device name', 'term')

        # set type of plot (PGPLOT or matplotlib) and the name of the file
        # if any in the case of matplotlib
        fslash = device.rfind('/')
        if fslash > -1:
            if device[fslash+1:] == 'mpl':
                ptype = 'MPL'
                hard = device[:fslash].strip()
            else:
                ptype = 'PGP'

        else:
            raise ValueError(
                'Could not indentify plot type from device = {:s}'.format(device)
                )

        # define the panel grid
        try:
            nxdef = cl.get_default('nx')
        except:
            nxdef = 3

        max_ccd = len(mccd)
        if max_ccd > 1:
            ccd = cl.get_value('ccd', 'CCD(s) to plot [0 for all]', '0')
            if ccd == '0':
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
            if len(ccds) > 1:
                nxdef = min(len(ccds), nxdef)
                cl.set_default('nx', nxdef)
                nx = cl.get_value('nx', 'number of panels in X', 3, 1)
            else:
                nx = 1
        else:
            ccds = list(mccd.keys())
            nx = 1

        # define the display intensities
        iset = cl.get_value(
            'iset', 'set intensity a(utomatically), d(irectly) or with p(ercentiles)?',
            'a', lvals=['a','A','d','D','p','P'])
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == 'd':
            ilo = cl.get_value('ilo', 'lower intensity limit', 0.)
            ihi = cl.get_value('ihi', 'upper intensity limit', 1000.)
        elif iset == 'p':
            plo = cl.get_value('plo', 'lower intensity limit percentile', 5., 0., 100.)
            phi = cl.get_value('phi', 'upper intensity limit percentile', 95., 0., 100.)

        nxmax, nymax = 0, 0
        for cnam in ccds:
            nxmax = max(nxmax, mccd[cnam].nxtot)
            nymax = max(nymax, mccd[cnam].nytot)

        if ptype == 'PGP' or hard != '':
            xlo = cl.get_value('xlo', 'left-hand X value', 0., 0., nxmax+1)
            xhi = cl.get_value('xhi', 'right-hand X value', float(nxmax), 0., nxmax+1)
            ylo = cl.get_value('ylo', 'lower Y value', 0., 0., nymax+1)
            yhi = cl.get_value('yhi', 'upper Y value', float(nymax), 0., nymax+1)
        else:
            xlo, xhi, ylo, yhi = 0, nxmax+1, 0, nymax+1

        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

    # all inputs obtained, plot
    if ptype == 'MPL':
        if width > 0 and height > 0:
            fig = plt.figure(figsize=(width,height))
        else:
            fig = plt.figure()

        nccd = len(ccds)
        ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

        ax = None
        for n, cnam in enumerate(ccds):
            if ax is None:
                axes = ax = fig.add_subplot(ny, nx, n+1)
                axes.set_aspect('equal', adjustable='box')
            else:
                axes = fig.add_subplot(ny, nx, n+1, sharex=ax, sharey=ax)
                axes.set_aspect('equal', adjustable='datalim')
            hcam.mpl.pccd(plt,mccd[cnam],iset,plo,phi,ilo,ihi)

            plt.title('CCD {:s}'.format(cnam))
            plt.xlabel('X')
            plt.ylabel('Y')

        plt.tight_layout()
        if hard == '':
            plt.show()
        else:
            plt.savefig(hard)

    elif ptype == 'PGP':
        # open the plot
        dev = hcam.pgp.Device(device)
        if width > 0 and height > 0:
            pgpap(width,height/width)

        nccd = len(ccds)
        ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

        # set up panels and axes
        pgsubp(nx,ny)

        for cnam in ccds:
            pgsci(hcam.pgp.Params['axis.ci'])
            pgsch(hcam.pgp.Params['axis.number.ch'])
            pgenv(xlo, xhi, ylo, yhi, 1, 0)
            pglab('X','Y','CCD {:s}'.format(cnam))
            hcam.pgp.pccd(mccd[cnam],iset,plo,phi,ilo,ihi)


#############################################################
#
# makedata -- generates fake multi-CCD data
#
#############################################################

def makedata(args=None):
    """Script to generate multi-CCD test data given a set of parameters defined in
    a config file (parsed using configparser). This allows things such as a
    bias frame, flat field variations offsets, scale factor and rotations
    between CCDs, and temporal variations.

    Arguments::

       config: (string)
          file defining the parameters.

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

    if args is None:
        args = sys.argv[1:]

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'makedata', args) as cl:

        # Register parameters
        cl.register('config', Cline.LOCAL, Cline.PROMPT)

        config = cl.get_value('config', 'configuration file',
                              cline.Fname('config'))

    # Read the config file
    conf = configparser.ConfigParser()
    conf.read(config)

    # Determine whether files get overwritten or not
    if 'general':
        overwrite = conf.getboolean('general', 'overwrite')
    else:
        overwrite = False

    # Top-level header
    thead = fits.Header()
    thead.add_history('Created by makedata')

    # Store the CCD labels and parameters and their dimensions. Determine
    # maximum dimensions for later use when adding targets
    ccd_pars = {}
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

            # determine maximum total dimension
            maxnx = max(maxnx, nxtot)
            maxny = max(maxny, nytot)

            # store parameters
            ccd_pars[int(key[3:])] = {
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
            }

    if not len(ccd_pars):
        raise ValueError('hipercam.makedata: no CCDs found in ' + config)

    # Generate the CCDs, store the read / gain values
    ccds = []
    rgs = {}
    for nccd, pars in ccd_pars.items():
        # Generate the Windats
        winds = []
        rgs[nccd] = {}
        for key in conf:
            if key.startswith('window'):
                iccd, nwin = key[6:].split()
                if int(iccd) == nccd:
                    nwin = int(nwin)
                    llx = int(conf[key]['llx'])
                    lly = int(conf[key]['lly'])
                    nx = int(conf[key]['nx'])
                    ny = int(conf[key]['ny'])
                    xbin = int(conf[key]['xbin'])
                    ybin = int(conf[key]['ybin'])
                    wind = hcam.Windat(hcam.Window(llx,lly,nx,ny,xbin,ybin))

                    # Accumulate Windats in the CCD
                    winds.append((nwin, wind))

                    # Store read / gain value
                    rgs[nccd][nwin] = (
                        float(conf[key]['read']), float(conf[key]['gain'])
                    )

        # Accumulate CCDs
        ccds.append((nccd, hcam.CCD(winds,pars['nxtot'],pars['nytot'])))

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

    # Now build up Fields
    fields = []
    for key in conf:
        if key.startswith('field'):
            field = hcam.Field()
            r = conf[key]
            ntarg = int(r['ntarg'])
            height1 = float(r['height1'])
            height2 = float(r['height2'])
            angle1 = float(r['angle1'])
            angle2 = float(r['angle2'])
            fwhm1 = float(r['fwhm1'])
            fwhm2 = float(r['fwhm2'])
            beta = float(r['beta'])
            fmin = float(r['fmin'])
            ndiv = int(r['ndiv'])

            field.add_random(ntarg,-5.,maxnx+5.,-5.,maxny+5.,
                             height1, height2, angle1, angle2,
                             fwhm1, fwhm2, beta, fmin)

            fields.append((ndiv, field))

    # Everything is set to go, so now generate data files
    nfiles = int(conf['files']['nfiles'])
    if nfiles == 0:
        out = mccd*flat + bias
        fname = hcam.add_extension(conf['files']['root'],hcam.HCAM)
        out.wfits(fname, overwrite)
        print('Written data to',fname)
    else:
        root = conf['files']['root']
        ndigit = int(conf['files']['ndigit'])

        print('Now generating data')

        for nfile in range(nfiles):

            # copy over template
            frame = mccd.copy()

            # add targets
            for nccd, ccd in frame.items():
                # get field modification settings
                p = ccd_pars[nccd]
                transform = Transform(
                    p['nxtot'], p['nytot'], p['xcen'], p['ycen'],
                    p['angle'], p['scale'], p['xoff'], p['yoff'])
                fscale = p['fscale']

                # wind through each window
                for nwin, wind in ccd.items():

                    # add targets from each field
                    for ndiv, field in fields:
                        mfield = field.modify(transform, fscale)
                        mfield.add(wind, ndiv)

            # Apply flat
            frame *= flat

            # Add noise
            for nccd, ccd in frame.items():
                for nwin, wind in ccd.items():
                    readout, gain = rgs[nccd][nwin]
                    wind.add_noise(readout, gain)

            # Apply bias
            frame += bias

            # Save
            fname = '{0:s}{1:0{2:d}d}{3:s}'.format(root,nfile+1,ndigit,hcam.HCAM)
            frame.wfits(fname, overwrite)
            print('Written file {0:d} to {1:s}'.format(nfile+1,fname))

#####################################################################
#
# makefield -- makes an artificial star field
#
#####################################################################

def makefield(args=None):
    """Script to generate an artificial star field which is saved to disk file, a
    first step in generating fake data. If the name supplied corresponds to an
    existing file, an attempt will be made to read it in first and then add to
    it. In this way a complex field can be generated. The targets are
    distributed at random, with random peak heights based on constant
    luminosity objects distributed throughout 3D space, and random angles
    uniform over the input range. All targets otherwise have the same shape
    thus multiple calls are needed to generate a field of objects of multiple
    shapes. Ellisoidal "Moffat" functions [1/(1+r^2)^beta] are used.

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

       fwmin : (float)
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

        # get inputs
        fname = cl.get_value('fname', 'file to save field to',
                             cline.Fname('field', hcam.FIELD, hcam.Fname.NEW))
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

    # add targets
    field.add_random(ntarg, x1, x2, y1, y2, h1, h2, angle1, angle2, fwmax, fwmin, beta)

    # save result
    field.wjson(fname)

    print('>> Saved a field of',len(field),'objects to',fname)

###############################################################
#
# rtplot -- plots images as they come in "real time" hence 'rt'
#
###############################################################

def rtplot(args=None):
    """Plots a sequence of images as a movie in near 'real time', hence
    'rt'. Designed to be used to look at images coming in while at the
    telescope.

    Arguments::

        device : (string)
          Plot device. PGPLOT is used so this should be a PGPLOT-style name,
          e.g. '/xs', '1/xs' etc. At the moment only ones ending /xs are supported.

        source : (string)
           's' = server, 'l' = local, 'f' = file list (ucm for ULTRACAM / ULTRASPEC,
           FITS files for HiPERCAM)

        inst   : (string)
           If 's' = server, 'l' = local, name of instrument. Choices: 'u' for
           ULTRACAM / ULTRASPEC, 'h' for HiPERCAM. This is needed because of the
           different formats.

        run    : (string)
           If source == 's' or 'l', run number to access, e.g. 'run034'

        flist  : (string)
           If source == 'f', name of file list

        first  : (int)
           If source='s' or 'l', this is the exposure to start from. 1 = first frame; set
           = 0 to always try to get the most recent frame (if it has changed)

        nccd   : (int)
           CCD number to plot, 0 for all.

        nx     : (int)
           number of panels across to display.

        bias   : (string)
           Name of bias frame to subtract, 'none' to ignore.

        iset   : (string) [single character]
           determines how the intensities are determined. There are three
           options: 'a' for automatic simply scales from the minimum to the
           maximum value found on a per CCD basis. 'd' for direct just takes
           two numbers from the user. 'p' for percentile dtermines levels
           based upon percentiles determined from the entire CCD on a per CCD
           basis.

        ilo    : (float) [if iset='d']
           lower intensity level

        ihi    : (float) [if iset='d']
           upper intensity level

        plo    : (float) [if iset='p']
           lower percentile level

        phi    : (float) [if iset='p']
           upper percentile level

        width  : (float) [hidden]
           plot width (inches). Set = 0 to let the program choose.

        height : (float) [hidden]
           plot height (inches). Set = 0 to let the program choose. BOTH width
           AND height must be non-zero to have any effect
    """

    if args is None:
        args = sys.argv[1:]

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', 'rtplot', args) as cl:

        # register parameters
        cl.register('device', Cline.LOCAL, Cline.HIDE)
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('inst', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)
        cl.register('first', Cline.LOCAL, Cline.PROMPT)
        cl.register('flist', Cline.LOCAL, Cline.PROMPT)
        cl.register('ccd', Cline.LOCAL, Cline.PROMPT)
        cl.register('nx', Cline.LOCAL, Cline.PROMPT)
        cl.register('bias', Cline.GLOBAL, Cline.PROMPT)
        cl.register('iset', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ilo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ihi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('plo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('phi', Cline.LOCAL, Cline.PROMPT)
        cl.register('xlo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('xhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('ylo', Cline.GLOBAL, Cline.PROMPT)
        cl.register('yhi', Cline.GLOBAL, Cline.PROMPT)
        cl.register('width', Cline.LOCAL, Cline.HIDE)
        cl.register('height', Cline.LOCAL, Cline.HIDE)

        # get inputs
        device = cl.get_value('device', 'plot device', '1/xs')

        source = cl.get_value('source', 'data source [s(erver), l(ocal), f(ile list)]',
                              'l', lvals=('s','l','f'))

        if source == 's' or source == 'l':
            inst = cl.get_value('inst', 'instrument [h(ipercam), u(ltracam/spec)]',
                                'h', lvals=('h','u'))
            run = cl.get_value('run', 'run name', 'run005')
            first = cl.get_value('first', 'first frame to plot', 1, 1)

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

        try:
            nxdef = cl.get_default('nx')
        except:
            nxdef = 3

        if len(ccdinf) > 1:
            ccd = cl.get_value('ccd', 'CCD(s) to plot [0 for all]', '0')
            if ccd == '0':
                ccds = list(ccdinf.keys())
            else:
                ccds = ccd.split()

            if len(ccds) > 1:
                nxdef = min(len(ccds), nxdef)
                cl.set_default('nx', nxdef)
                nx = cl.get_value('nx', 'number of panels in X', 3, 1)
            else:
                nx = 1
        elif len(ccdinf) == 1:
            nx = 1
            ccds = list(ccdinf.keys())

        # bias frame (if any)
        bias = cl.get_value(
            'bias', "bias frame ['none' to ignore]",
            cline.Fname('bias', hcam.HCAM), ignore='none'
        )
        if bias is not None:
            # read the bias frame
            bframe = hcam.MCCD.rfits(bias)

        # define the display intensities
        iset = cl.get_value(
            'iset', 'set intensity a(utomatically), d(irectly) or with p(ercentiles)?',
            'a', lvals=['a','A','d','D','p','P'])
        iset = iset.lower()

        plo, phi = 5, 95
        ilo, ihi = 0, 1000
        if iset == 'd':
            ilo = cl.get_value('ilo', 'lower intensity limit', 0.)
            ihi = cl.get_value('ihi', 'upper intensity limit', 1000.)
        elif iset == 'p':
            plo = cl.get_value('plo', 'lower intensity limit percentile', 5., 0., 100.)
            phi = cl.get_value('phi', 'upper intensity limit percentile', 95., 0., 100.)

        nxmax, nymax = 0, 0
        for cnam in ccds:
            nxtot, nytot = ccdinf[cnam]
            nxmax = max(nxmax, nxtot)
            nymax = max(nymax, nytot)

        xlo = cl.get_value('xlo', 'left-hand X value', 0., 0., nxmax+1)
        xhi = cl.get_value('xhi', 'right-hand X value', float(nxmax), 0., nxmax+1)
        ylo = cl.get_value('ylo', 'lower Y value', 0., 0., nymax+1)
        yhi = cl.get_value('yhi', 'upper Y value', float(nymax), 0., nymax+1)

        width = cl.get_value('width', 'plot width (inches)', 0.)
        height = cl.get_value('height', 'plot height (inches)', 0.)

    # Open device
    imdev = hcam.pgp.Device(device)
    if width > 0 and height > 0:
        pgpap(width,height/width)

    # set up panels and axes
    nccd = len(ccds)
    ny = nccd // nx if nccd % nx == 0 else nccd // nx + 1

    pgsubp(nx,ny)

    for cnam in ccds:
        pgsci(hcam.pgp.Params['axis.ci'])
        pgsch(hcam.pgp.Params['axis.number.ch'])
        pgenv(xlo, xhi, ylo, yhi, 1, 0)
        pglab('X','Y','CCD {:s}'.format(cnam))
    pgpanl(1,1)

    # Finally plot stuff
    with hcam.data_source(instrument, run, flist, server, first, True) as spool:
        n = 1
        for frame in spool:
            print('Exposure {:d}'.format(n+1))
            n += 1
            for cnam in ccds:
                # subtract bias
                if bias is not None:
                    frame[cnam] -= bframe[cnam]
                # plot CCD
                hcam.pgp.pccd(frame[cnam],iset,plo,phi,ilo,ihi)
                pgpage()
            pgpanl(1,1)

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

        # Retunr change in position
        return (xcen-x,ycen-y)
