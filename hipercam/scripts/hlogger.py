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
############################################

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
<link rel="stylesheet" type="text/css" href="../hiper.css" />
<title>HiPERCAM log {0:s}</title>
</head>

<body>
<h1>HiPERCAM log {0:s}</h1>

<p>
The table below lists information on runs from the night starting on {0:s}. See the end of the table for some more details on the meanings of the various columns.

<p>
<table>
<tr>
<th class="left">Run<br>no.</th>
<th class="left">Target<br>name</th>
<th class="left">RA (J2000)<br>Dec (J2000)<br>(telescope)</th>
<th class="cen">Start<br>time</th>
<th class="right">Dwell<br>(sec)</th>
<th class="right">Nframe</th>
<th class="left">Ncycle</th>
<th class="cen">Read<br>mode</th>
<th class="cen">xsll,xslr,xsul,xsur,ys,nx,ny<br>Quad1<br>Quad2</th>
<th class="cen">CDOP<br>flags</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
"""

NIGHT_FOOTER = """
</table>

<p> Along with others of more obvious meaning, the columns 'Read', 'Format and
'CDOP' specify the information needed to repeat the setup at the telescope
using 'hdriver'. 'Read' is the readout mode which can be one of 4 options,
namely 'FULL' for a full frame, 1-WIN or 2-WIN for standard windows mode, and
'DRIFT' for drift-mode. 'Format' column gives the 7 parameters that appear at
the bottom of the top-right window of hdriver. In hdriver these are called
xsll, xslr, xsul, xsur, ys, nx and ny, and there are up to two sets of them.
'CDOP' are four Y/N flags which say whether or not clears, the dummy output,
overscans and prescans were enabled or not. If overscans were enabled, extra
pixels are added in the Y direction. Prescans add extra pixels in X. 'Ncycle'
refers to how often each CCD is read out. The numbers correspond to nu, ng,
nr, ni, nz listed in hdriver.  </p>

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

# HiPERCAM CSS to define the look of the log web pages
# This is written to 'hiper.css' and referred to in the
# html pages
CSS = """

body{
    background-color: #000000;
    font: 11pt sans-serif;
    color: #FFFFFF;
    margin: 10px;
    border: 0px;
    overflow: auto;
    height: 100%;
    max-height: 100%;
}

/* This for the left-hand side guide */

#guidecontent{
    position: absolute;
    margin: 10px;
    top: 0;
    bottom: 0;
    left: 0;
    width: 200px;
    height: 100%;
    overflow: auto;
}

#titlecontent{
    position: fixed;
    margin: 10px;
    top: 0;
    left: 200px;
    right: 0;
    bottom: 220px;
    overflow: auto;
}

#logcontent{
    position: fixed;
    margin: 10px;
    top: 220px;
    left: 200px;
    right: 0;
    bottom: 0;
    overflow: auto;
}

p {color: #ffffe0}

h1 {color: #ffffff}

input.text {margin-right: 20px; margin-left: 10px}

spa {margin-right: 20px; margin-left: 20px}

table {
    font: 9pt sans-serif;
    padding: 1px;
    border-top:1px solid #655500;
    border-right:1px solid #655500;
    border-bottom:2px solid #655500;
    border-left:1px solid #655500;
}

td {
    vertical-align: top;
    text-align: center;
    white-space: nowrap;
    padding-right: 10px;
}

td.left {
    vertical-align: top;
    text-align: left;
    white-space: nowrap;
    padding-right: 10px;
}

td.right {
    vertical-align: top;
    text-align: right;
    white-space: nowrap;
    padding-right: 10px;
}

td.cen {
    vertical-align: top;
    text-align: center;
    white-space: nowrap;
}

td.undef {
    background-color: #100000;
}

td.format {
    vertical-align: top;
    text-align: center;
    white-space: nowrap;
}

td.long {
    vertical-align: top;
    text-align: left;
    padding-right: 10px;
    font: 12pt sans-serif;
    white-space: normal;
}

td.bleft {
    color: #ffffa0;
    vertical-align: top;
    text-align: left;
    white-space: nowrap;
    padding-right: 10px;
    font: 12pt sans-serif;
    font-weight: bold;
}

th {
    vertical-align:top;
    text-align: center;
}

th.left {
    vertical-align: top;
    text-align: left;
}

th.cen {
    vertical-align: top;
    text-align: center;
}

/* links */

a:link {
    color: #7070ff;
    text-decoration:underline;
    font-size: 11pt;
}

a:visited {
    color: #e0b0e0;
    text-decoration:underline;
    font-size: 11pt;
}

a:hover {
    color: red;
    text-decoration:underline;
    font-size: 11pt;
}
"""

TRANSLATE_MODE = {
    'FullFrame' : 'FULL',
    'OneWindow' : '1-WIN',
    'TwoWindows' : '2-WIN',
    'DriftWindow' : 'DRIFT',
}

def correct_ra_dec(ra, dec):
    """Fixes up RAs and Decs"""
    if ra == 'UNDEF':
        ra = ''
    else:
        rah,ram,ras = ra.split(':')
        rah, ram, ras = int(rah), int(ram), float(ras)
        ra = '{:02d}:{:02d}:{:05.2f}'.format(rah,ram,ras)

    if dec == 'UNDEF':
        dec = ''
    else:
        decd,decm,decs = dec.split(':')
        if decd.find('-') > -1:
            negative = True
        else:
            negative = False
        decd, decm, decs = int(decd), int(decm), float(decs)
        decd = -abs(decd) if negative else abs(decd)
        dec = '{:+03d}:{:02d}:{:04.1f}'.format(decd,decm,decs)

    return (ra, dec)

def hlogger(args=None):
    """Generates html logs for hipercam runs.

    hlogger expects to work on directories containing the runs for each
    night. It should be run from the top-level of the hipercam/logs directory.
    It extracts information from each run file which it writes to a file
    in JSON format and also as an html file. The JSON file acts as a short-cut
    for any re-run.

    Arguments::

        server   : bool
           True if the log is to be compiled from the server data. This option
           might be useful during observing. If server is True, the JSON and
           html files will be written in the current working directory.

        dirnam : string [if not server]
           Directory name to work on. This is run through 'glob' which gives
           UNIX-style name matching, so e.g. '2017*' will match 2017_10_17,
           2017_10_18 etc. The JSON and html files will be written into the
           directory or directories. In this case this expects the data to be
           available as a link called 'data' within the night directory. This
           is the standard configuration of my 'logs' directory as created by
           'digest'. This script should be run from the log directory. To run
           it on the command line, to prevent the shell from globbing, you
           must escape or quote the glob as in "hlogger no '*'" or "hlogger no
           \*".

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
            print("there were no run directories of the form "
                  "YYYY-MM with a file called 'telescope' in them", file=sys.stderr)
            print('hlogger aborted',file=sys.stderr)
            return

        with open('index.html','w') as ihtml:
            # start the top level index html file
            ihtml.write(INTRODUCTION_HEADER)

            for rname in rnames:
                print('\nProcessing run {:s}'.format(rname))

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
                        "found no night directories of the form YYYY-MM-DD with"
                        " 'data' sub-directories in {:s}".format(rname),
                        file=sys.stderr
                        )
                    print('hlogger aborted',file=sys.stderr)
                    return

                for nname in nnames:
                    print('  night {:s}'.format(nname))

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

                        # now wind through the runs getting basic info and
                        # writing a row of info to the html file for the night
                        # in question.
                        for run in runs:

                            # open the run file as an Rtime
                            rname = os.path.join(night,'data',run)
                            rtime = hcam.hcam.Rtime(rname)
                            hd = rtime.header

                            # start the row
                            nhtml.write('<tr>\n')

                            # run number
                            nhtml.write('<td class="left">{:s}</td>'.format(run[3:]))

                            # object name
                            nhtml.write(
                                '<td class="left">{:s}</td>'.format(
                                    hd['OBJECT'])
                            )

                            # RA, Dec
                            ra, dec = correct_ra_dec(hd['RA'], hd['Dec'])
                            nhtml.write(
                                '<td class="left">{:s}<br>{:s}</td>'.format(
                                    ra, dec)
                            )

                            # First timestamp
                            try:
                                tstamp = rtime(1).isot
                                nhtml.write(
                                    '<td class="cen">{:s}</td>'.format(
                                        tstamp[tstamp.find('T')+1:tstamp.rfind(
                                            '.')])
                                )
                            except:
                                nhtml.write('<td class="cen">--:--:--</td>')

                            # timing info
                            texps, toffs, ncycs, tdead = rtime.tinfo()
                            ttotal = 0.
                            ntotal = rtime.ntotal()
                            for texp, ncyc in zip(texps, ncycs):
                                ttotal = max(
                                    ttotal, (texp+tdead)*(ntotal // ncyc)
                                )

                            # total exposure time
                            nhtml.write('<td class="right">{:d}</td>'.format(
                                int(round(ttotal)))
                            )

                            # number of frames
                            nhtml.write(
                                '<td class="right">{:d}</td>'.format(ntotal)
                            )

                            # cycle nums
                            nhtml.write('<td class="left">{:s}</td>'.format(
                                '|'.join([str(ncyc) for ncyc in ncycs]))
                            )

                            # readout mode
                            nhtml.write('<td class="cen">{:s}</td>'.format(
                                    TRANSLATE_MODE[rtime.mode])
                                        )

                            # Format
                            nhtml.write(
                                '<td class="cen">{:s}</td>'.format(
                                    '<br>'.join(rtime.wforms))
                                )

                            # flags
                            nhtml.write(
                                '<td class="cen">{:s}{:s}{:s}{:s}</td>'.format(
                                    'Y' if rtime.clear else 'N',
                                    'Y' if rtime.dummy else 'N',
                                    'Y' if rtime.oscan else 'N',
                                    'Y' if rtime.pscan else 'N')
                                )

                            # run number again
                            nhtml.write(
                                '<td class="left">{:s}</td>'.format(run[3:])
                            )

                            # hand log comment
                            nhtml.write(
                                '<td class="left">{:s}</td>'.format(hlog[run])
                            )

                            # end the row
                            nhtml.write('\n</tr>\n')

                        # finish off the night file
                        nhtml.write(NIGHT_FOOTER)

                # finish off the run table
                ihtml.write('\n</table>\n\n')

            # finish the main index
            ihtml.write(INTRODUCTION_FOOTER)


    # write out the css file
    with open('hiper.css','w') as fout:
        fout.write(CSS)
