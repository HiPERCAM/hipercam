import sys
import os
import time
import glob
import re

import numpy as np
from astropy.time import Time, TimeDelta

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ = ['hlogger',]

############################################
#
# logger -- generates html logs for hipercam
#
###############################

def hlogger(args=None):
    """Generates html logs for hipercam runs.

    hlogger expects to work on directories containing the runs for each night.
    It then extract information from each run file which it writes to a file
    in JSON format and also as an html file. The JSON file acts as a short-cut
    for any re-run.

    Arguments::

        server   : bool
           True if the log is to be compiled from the server data. This option
           might be useful during observing. If server is True, the JSON and
           html files will be written in the current working directory.

        dirnam   : string [if not server]
           Directory name to work on. This is run through 'glob' which gives
           UNIX-style name matching, so e.g. '2017*' will match 2017_10_17,
           2017_10_18 etc. The JSON and html files will be written into the
           directory or directories.

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('server', Cline.GLOBAL, Cline.PROMPT)
        cl.register('dirnam', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        server = cl.get_value('server', 'access data via the server?', False)
        if not server:
            dirnam = cl.get_value('dirnam', 'directory to log', '2017-10-18')

    if not server:

        # get matching directories
        dnames = [dname for dname in glob.glob(dirnam) if os.path.isdir(dname)]

        if len(dnames) == 0:
            print('there were no directories matching the pattern = {:s}'.format(dirnam),
                  file=sys.stderr)
            return

        # match runs of specific form 'run1234.fits' (4 digits)
        rre = re.compile('^run\d\d\d\d\.fits$')

        for dname in dnames:
            runs = [run for run in os.listdir(dname) if rre.match(run)]
            print(dname, runs)

            
