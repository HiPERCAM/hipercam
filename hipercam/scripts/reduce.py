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

    # plot images
    with hcam.data_source(instrument, run, flist, server, first) as spool:

        # 'spool' is an iterable source of MCCDs
        for n, mccd in enumerate(spool):

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

            # do something else ...


# From here is support code not visible outside

class Rfile(configparser.ConfigParser):
    """Class to read and interpret reduce scripts."""

    # Data
    VERSION = '2017-10-05'

    @classmethod
    def fromFile(self, filename):
        # read it in, stripping off trailing comments
        rfile = configparser.ConfigParser(inline_comment_prefixes='#')
        with open(filename) as fp:
            rfile.read_file(fp)

        if rfile['general']['version'] != Rfile.VERSION:
            # check the version
            raise ValueError(
                'Version mismatch: file = {:s}, reduce = {:s}'.format(
                    rfile['general']['version'], Rfile.VERSION)
            )

        return rfile
