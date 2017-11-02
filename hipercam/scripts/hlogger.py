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

###########
#
# Some constant stuff
#

INTRODUCTION_HEADER = """<html>
<head>
<title>HiPERCAM data logs</title>
</head>

<body>
<h1>HiPERCAM data logs</h1>

<p>
These online logs summarise HiPERCAM runs and were automatically generated from the data files and hand-written logs.

"""

INTRODUCTION_FOOTER = """
<address>Tom Marsh, Warwick</address>
</body>
</html>
"""

NIGHT_HEADER = """<html>
<head>
<title>HiPERCAM log {0:s}</title>
</head>

<body>
<h1>HiPERCAM log {0:s}</h1>

<p>
<table>
"""

NIGHT_FOOTER = """
</table>

<address>Tom Marsh, Warwick</address>
</body>
</html>
"""

MONTHS = {
    '01' : 'January',
    '02' : 'February',
    '03' : 'March',
    '04' : 'April',
    '05' : 'May',
    '06' : 'June',
    '07' : 'July',
    '08' : 'August',
    '09' : 'September',
    '10' : 'October',
    '11' : 'November',
    '12' : 'December',
}

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

        dirnam : string [if not server] Directory name to work on. This is run
           through 'glob' which gives UNIX-style name matching, so
           e.g. '2017*' will match 2017_10_17, 2017_10_18 etc. The JSON and
           html files will be written into the directory or directories. In
           this case this expects the data to be available as a link called
           'data' within the night directory. This is the standard
           configuration of my 'logs' directory as created by 'digest'. This
           script should be run from the log directory. To run it on the command
           line, to prevent the shell from globbing, you must escape or quote the
           glob as in "hlogger no '*'" or "hlogger no \*".

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('server', Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        server = cl.get_value('server', 'access data via the server?', False)

    if server:
        raise NotImplementedError('server option yet to be enabled')

    else:

        if os.path.basename(os.getcwd()) != 'logs':
            print("** hlogger must be run in a directory called 'logs'")
            print('hlogger aborted',file=sys.stderr)
            return

        # next are regular expressions to match run directories, nights, and run files
        rre = re.compile('^\d\d\d\d-\d\d$')
        nre = re.compile('^\d\d\d\d-\d\d-\d\d$')
        fre = re.compile('^run\d\d\d\d\.fits$')

        # Get list of run directories
        rnames = [rname for rname in os.listdir('.') if \
                      rre.match(rname) and \
                      os.path.isdir(rname) and \
                      os.path.isfile(os.path.join(rname, 'telescope'))
                  ]

        rnames.sort()

        if len(rnames) == 0:
            print("there were no run directories of the form YYYY-MM with a file called 'telescope' in them", file=sys.stderr)
            print('hlogger aborted',file=sys.stderr)
            return

        with open('index.html','w') as ihtml:
            # start the top level index html file
            ihtml.write(INTRODUCTION_HEADER)

            for rname in rnames:
                # write in run date, start table of nights
                rn = os.path.basename(rname)
                year, month = rn.split('-')
                with open(os.path.join(rname,'telescope')) as tel:
                    telescope = tel.read().strip()
                ihtml.write('<h2>{:s} {:s}, {:s}</h2>\n'.format(MONTHS[month],year,telescope))
                ihtml.write('\n<p>\n<table>\n')

                # get night directories
                nnames = [os.path.join(rname, ndir) for ndir in os.listdir(rname) if \
                              nre.match(ndir) and \
                              os.path.isdir(os.path.join(rname, ndir)) and \
                              os.path.isdir(os.path.join(rname, ndir, 'data'))
                          ]
                nnames.sort()

                if len(nnames) == 0:
                    print(
                        "found no night directories of the form YYYY-MM-DD with 'data' sub-directories in {:s}".format(rname),
                        file=sys.stderr
                        )
                    print('hlogger aborted',file=sys.stderr)
                    return

                for nname in nnames:
                    # Write an entry for each night linking to the log for that night.
                    night = os.path.basename(nname)
                    fname = '{0:s}/{0:s}.html'.format(night)
                    ihtml.write('<tr><td><a href="{:s}">{:s}</a><td></tr>\n'.format(fname, night))

                    # Create the html file for the night
                    date = '{:s}, {:s}'.format(night, telescope)

                    with open(fname,'w') as nhtml:
                        # write header of night file
                        nhtml.write(NIGHT_HEADER.format(date))

                        # read and store the hand written log
                        handlog = os.path.join(night,'{:s}.dat'.format(night))
                        with open(handlog) as fin:
                            hlog = {}
                            for line in fin:
                                if line.startswith('run'):
                                    arr = line.split()
                                    hlog[arr[0]] = ' '.join(arr[1:])

                        # load all the run names
                        ddir = os.path.join(night,'data')
                        runs = [run[:-5] for run in os.listdir(ddir) if \
                                    fre.match(run)]
                        runs.sort()

                        for run in runs:
                            nhtml.write('<tr><td>{:s}</td><td>{:s}</td></tr>\n'.format(run,hlog[run]))

                        # finish off the night file
                        nhtml.write(NIGHT_FOOTER)

                # finish off the run table
                ihtml.write('\n</table>\n\n')

            # finish the main index
            ihtml.write(INTRODUCTION_FOOTER)


