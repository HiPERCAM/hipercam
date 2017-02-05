import sys
import os
import math
from collections import OrderedDict

import numpy as np
from astropy.io import fits

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

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

    # create a Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', command, args)

    # register parameters
    cl.register('input', Cline.LOCAL, Cline.PROMPT)
    cl.register('const', Cline.LOCAL, Cline.PROMPT)
    cl.register('output', Cline.LOCAL, Cline.PROMPT)

    try:
        prompts = {'cadd' : 'add', 'csub' : 'subtract',
                   'cmul' : 'multiply by', 'cdiv' : 'divide by'}

        # get inputs
        infile = cl.get_value('input', 'input file',
                              cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(infile)

        constant = cl.get_value('const', 'constant to ' + prompts[command], 0.)

        outfile = cl.get_value('output', 'output file',
                               cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW))

    except cline.ClineError as err:
        sys.stderr.write('Error on parameter input:\n{0:s}\n'.format(str(err)))
        sys.exit(1)

    import matplotlib.pyplot as plt

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
    mccd.head.add_history('{0:s} {1:s} {2:f} {3:s}'.format(command,infile,constant,outfile))

    # save result
    mccd.wfits(outfile, True)

def hplot(args=None):
    """
    Plots a multi-CCD image.

    """

    if args is None:
        args = sys.argv[1:]

    # create a Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'hplot', args)

    # register parameters
    cl.register('input', Cline.LOCAL, Cline.PROMPT)
    cl.register('nccd', Cline.LOCAL, Cline.PROMPT)
    cl.register('nx', Cline.LOCAL, Cline.PROMPT)

    try:
        # get inputs
        frame = cl.get_value('input', 'frame to plot',
                             cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(frame)
        max_ccd = len(mccd)
        if max_ccd > 1:
            nccd = cl.get_value('nccd', 'CCD number to plot', 0, 0)
            if nccd == 0:
                nx = cl.get_value('nx', 'number of panels in X', 3, 1)
            else:
                nx = 1
        else:
            nccd = 1

    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        sys.exit(1)

    import matplotlib.pyplot as plt

    if nccd == 0:
        ny = int(math.ceil(max_ccd / nx))
        fig = plt.figure()
        i = 1
        for label, ccd in mccd.items():
            if i == 1:
                ax = fig.add_subplot(ny,nx,i)
                asave = ax
            else:
                ax = fig.add_subplot(ny,nx,i,sharex=asave,sharey=asave)
            i += 1
            hcam.mpl.pccd(ax, ccd)
            plt.xlim(0.5,ccd.nxtot+0.5)
            plt.ylim(0.5,ccd.nytot+0.5)
            plt.title('CCD ' + str(label))
            plt.xlabel('X pixels')
            plt.ylabel('Y pixels')
            plt.tight_layout(pad=1.01)
    else:
        try:
            ccd = mccd[nccd]
        except KeyError:
            sys.stderr.write('No CCD number {0:d} found in file = {1:s}\n'.format(nccd,frame))
            sys.exit(1)

        hcam.mpl.pccd(plt, ccd)
        plt.xlim(0.5,ccd.nxtot+0.5)
        plt.ylim(0.5,ccd.nytot+0.5)
        plt.title('CCD ' + str(nccd))
        plt.xlabel('X pixels')
        plt.ylabel('Y pixels')

    plt.show()

def makedata(args=None):
    """Script to create multi-CCD data given an artificial star field (generated
    e.g. with `makefield`), a bias and a matching flat field frame. This
    allows a sequence of images to be generated allowing for drifting with
    time, jitter and variable transparency.

    Arguments::

       field : (string)
          star field

       bias : (string)
          bias frame

       flat : (string)
          flat field frame

    """
    if args is None:
        args = sys.argv[1:]

    # create Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'makefield', args)

    # register parameters
    cl.register('field', Cline.LOCAL, Cline.PROMPT)
    cl.register('bias', Cline.LOCAL, Cline.PROMPT)
    cl.register('flat', Cline.LOCAL, Cline.PROMPT)
    cl.register('output', Cline.LOCAL, Cline.PROMPT)

    try:
        field = cl.get_value('field', 'star field file', cline.Fname('field', hcam.FIELD))
        bias = cl.get_value('bias', 'bias frame', cline.Fname('bias', hcam.HCAM))
        flat = cl.get_value('flat', 'flat field', cline.Fname('flat', hcam.FIELD))
        outfile = cl.get_value('output', 'output frame',
                                cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW))

???????????
    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    # add targets
    field.add_random(ntarg, x1, x2, y1, y2, h1, h2, fwmax, fwmin, angle, beta)

    # save result
    field.wjson(fname)

    print('>> Saved a field of',len(field),'objects to',fname)

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
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'makefield', args)

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

    try:
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

    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    # add targets
    field.add_random(ntarg, x1, x2, y1, y2, h1, h2, angle1, angle2, fwmax, fwmin, beta)

    # save result
    field.wjson(fname)

    print('>> Saved a field of',len(field),'objects to',fname)


def makehcam(args=None):
    """Script to generate a fake hipercam frame. Use this to create artificial
    bias and flat field frames.  This will split each CCD up into 2x2 windows,
    so the total dimensions must be multiples of 2.

    Arguments::

       frame : (string)
          name of frame [output]

       nccd : (int)
          number of CCDs

       nxtot : (int)
          total X dimension, a multiple of 2 [unbinned pixels]

       nytot : (int)
          total Y dimension, a multiple of 2 [unbinned pixels]

       mean : (float)
          mean level of frame [counts]

       sigma : (float)
          RMS gaussian noise level [counts]

       gradx : (float)
          change in level in X-direction from left-to-right [counts]

       grady : (float)
          change in level in Y-direction from bottom-to-top [counts]

    """
    if args is None:
        args = sys.argv[1:]

    # create Cline object
    cl = Cline('HIPERCAM_ENV', '.hipercam', 'makehcam', args)

    # register parameters
    cl.register('frame', Cline.LOCAL, Cline.PROMPT)
    cl.register('nccd', Cline.LOCAL, Cline.PROMPT)
    cl.register('nxtot', Cline.LOCAL, Cline.PROMPT)
    cl.register('nytot', Cline.LOCAL, Cline.PROMPT)
    cl.register('xbin', Cline.LOCAL, Cline.PROMPT)
    cl.register('ybin', Cline.LOCAL, Cline.PROMPT)
    cl.register('mean', Cline.LOCAL, Cline.PROMPT)
    cl.register('sigma', Cline.LOCAL, Cline.PROMPT)
    cl.register('gradx', Cline.LOCAL, Cline.PROMPT)
    cl.register('grady', Cline.LOCAL, Cline.PROMPT)

    try:
        # get cls
        frame = cl.get_value('frame', 'name of output frame',
                                cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW))
        nccd = cl.get_value('nccd', 'number of CCDs/frame', 5, 1)
        nxtot = cl.get_value('nxtot', 'number of CCDs/frame', 2048, 2, multipleof=2)
        nytot = cl.get_value('nytot', 'number of CCDs/frame', 1024, 2, multipleof=2)
        xbin = cl.get_value('xbin', 'X-binning factor', 1, 1)
        ybin = cl.get_value('ybin', 'Y-binning factor', 1, 1)
        mean = cl.get_value('mean', 'mean counts/pixel', 0.)
        sigma = cl.get_value('sigma', 'RMS counts', 0., 0.)
        gradx = cl.get_value('gradx', 'left-to-right change, counts', 0.)
        grady = cl.get_value('grady', 'bottom-to-top change, counts', 0.)

    except cline.ClineError as err:
        print('Error on parameter input:')
        print(err)
        exit(1)

    nx = nxtot // 2 // xbin
    ny = nytot // 2 // ybin

    print('Creating 4 Windows at corners of CCD(s)')
    win1  = hcam.Window(1,1,nx,ny,xbin,ybin)
    win2  = hcam.Window(nxtot-xbin*nx+1,1,nx,ny,xbin,ybin)
    win3  = hcam.Window(1,nytot-ybin*ny+1,nx,ny,xbin,ybin)
    win4  = hcam.Window(nxtot-xbin*nx+1,nytot-ybin*ny+1,nx,ny,xbin,ybin)

    print('Creating equivalent Windats')
    winds = (hcam.Windat(win1), hcam.Windat(win2), hcam.Windat(win3), hcam.Windat(win4))

    def func(x,y):
        xcen = (1+nxtot)/2.
        ycen = (1+nytot)/2.
        return gradx*(x-xcen)/(nxtot-1) + grady*(y-ycen)/(nytot-1)

    print('Adding mean, slope and noise')
    for wind in winds:
        wind.add_fxy(func)
        wind.data += np.random.normal(mean,sigma,(ny,nx))

    # Top-level header
    head = fits.Header()
    head['OBJECT'] = ('Fake star field', 'Target name')
    head['INSTRUME'] = ('HiPERCAM', 'Instrument name')
    head['RUN'] = (3,'Run number')
    head['MEAN'] = (mean, 'mean level used by makehcam')
    head['SIGMA'] = (sigma, 'RMS used by makehcam')
    head['GRADX'] = (gradx, 'Change in X used by makehcam')
    head['GRADY'] = (grady, 'Change in Y used by makehcam')
    head.add_history('Created by makehcam')

    # Prepare Windats to be made into CCDs
    winds = OrderedDict(zip(range(1,nccd+1),winds))

    ihead = fits.Header()

    print('Creating CCDs')
    ccds = OrderedDict()

    ihead['FILTER'] = ('u','Filter name')
    ccds[1] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('g','Filter name')
    ccds[2] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('r','Filter name')
    ccds[3] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('i','Filter name')
    ccds[4] = hcam.CCD(winds,nxtot,nytot,ihead)
    ihead['FILTER'] = ('z','Filter name')
    ccds[5] = hcam.CCD(winds,nxtot,nytot,ihead)

    print('Creating the MCCD')
    mccd = hcam.MCCD(ccds, head)

    print('Writing the MCCD')
    mccd.wfits(frame,True)


