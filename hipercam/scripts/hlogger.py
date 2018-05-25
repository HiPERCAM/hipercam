import sys
import traceback
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
These online logs summarise HiPERCAM runs and were automatically generated from
the data files and hand-written logs.

"""

INTRODUCTION_FOOTER = """
<address>Tom Marsh, Warwick</address>
</body>
</html>
"""

NIGHT_HEADER1 = """<html>
<head>

<!-- script to hide/reveal selected columns -->
<script type="text/javascript">

function hide ( column ) {
var tbl = document.getElementById( "tbl" );
var i;
for ( i = 0; i < tbl.rows.length; i++ )
        tbl.rows[ i ].cells[ column ].style.display = "none";
}

function restore () {
var tbl = document.getElementById( "tbl" );
var i;
var j;
for ( i = 0; i < tbl.rows.length; i++ )
        for ( j = 0; j < tbl.rows[ i ].cells.length; j++ )
                tbl.rows[ i ].cells[ j ].style.display = "table-cell";
}

function shrink () {
var hcols = [7, 16, 17, 18, 20, 21, 22];
var tbl = document.getElementById( "tbl" );
var i;
var j;
for ( i = 0; i < tbl.rows.length; i++ )
        for ( j = 0; j < hcols.length; j++ )
                tbl.rows[ i ].cells[ hcols[j] ].style.display = "none";
}

</script>

<!-- script to toggle auto refresh -->

<script>
var reloading;
function checkReloading() {
    if (window.location.hash=="#autoreload") {
        reloading=setTimeout("window.location.reload();",10000);
        document.getElementById("reloadCB").checked=true;
    }
}
function toggleAutoRefresh(cb) {
    if (cb.checked) {
        window.location.replace("#autoreload");
        reloading=setTimeout("window.location.reload();",10000);
    } else {
        window.location.replace("#");
        clearTimeout(reloading);
    }
}
window.onload=checkReloading;
</script>

<link rel="stylesheet" type="text/css" href="../hiper.css" />
"""

NIGHT_HEADER2 = """<title>HiPERCAM log {0:s}</title>
</head>

<body>
<h1>HiPERCAM log {0:s}</h1>

<p>
The table below lists information on runs from the night starting on {0:s}.
See the end of the table for some more details on the meanings of the various columns.
You can hide any column using the buttons in the top line and "shrink" will a standard
set of below average interest.
"""

TABLE_TOP = """<p>
<button onclick="restore();">Restore columns</button>
<button onclick="shrink();">Shrink table</button>

<p>
<table border=1 cellspacing="0" cellpadding="4" id="tbl">
"""

TABLE_HEADER = """
<tr>
<td><button id="hidden0" onclick="hide(0)"></button></td>
<td><button id="hidden1" onclick="hide(1)"></button></td>
<td><button id="hidden2" onclick="hide(2)"></button></td>
<td><button id="hidden3" onclick="hide(3)"></button></td>
<td><button id="hidden4" onclick="hide(4)"></button></td>
<td><button id="hidden5" onclick="hide(5)"></button></td>
<td><button id="hidden6" onclick="hide(6)"></button></td>
<td><button id="hidden7" onclick="hide(7)"></button></td>
<td><button id="hidden8" onclick="hide(8)"></button></td>
<td><button id="hidden9" onclick="hide(9)"></button></td>
<td><button id="hidden10" onclick="hide(10)"></button></td>
<td><button id="hidden11" onclick="hide(11)"></button></td>
<td><button id="hidden12" onclick="hide(12)"></button></td>
<td><button id="hidden13" onclick="hide(13)"></button></td>
<td><button id="hidden14" onclick="hide(14)"></button></td>
<td><button id="hidden15" onclick="hide(15)"></button></td>
<td><button id="hidden16" onclick="hide(16)"></button></td>
<td><button id="hidden17" onclick="hide(17)"></button></td>
<td><button id="hidden18" onclick="hide(18)"></button></td>
<td><button id="hidden19" onclick="hide(19)"></button></td>
<td><button id="hidden20" onclick="hide(20)"></button></td>
<td><button id="hidden21" onclick="hide(21)"></button></td>
<td><button id="hidden22" onclick="hide(22)"></button></td>
<td><button id="hidden23" onclick="hide(23)"></button></td>
<td><button id="hidden24" onclick="hide(24)"></button></td>
<td><button id="hidden25" onclick="hide(25)"></button></td>
<td align="left"><button id="hidden26" onclick="hide(26)"></button></td>
</tr>

<tr>
<th class="left">Run<br>no.</th>
<th class="left">Target<br>name</th>
<th class="left">RA (J2000)<br>(telescope)</th>
<th class="left">Dec&nbsp;(J2000)<br>(telescope)</th>
<th class="cen">Start</th>
<th class="cen">End</th>
<th class="right">Cadence<br>(sec)</th>
<th class="right">Duty<br>cycle, %</th>
<th class="right">Nframe</th>
<th class="right">Dwell<br>(sec)</th>
<th class="left">Type</th>
<th class="cen">Read<br>mode</th>
<th class="cen">Nskips</th>
<th class="cen">
xll,xlr,xul,xur,ys,nx,ny<br>
xsl,xsr,ys,nx,ny [DRIFT]<br>
Quad1</th>
<th class="cen">
xll,xlr,xul,xur,ys,nx,ny<br>
xsl,xsr,ys,nx,ny [DRIFT]<br>
Quad2</th>
<th class="cen">XxY<br>bin</th>
<th class="cen">Clr</th>
<th class="cen">Dum</th>
<th class="cen">Ov-/Pre-<br>scan</th>
<th class="cen">Readout<br>speed</th>
<th class="cen">Fast<br>clocks</th>
<th class="cen">Tbytes</th>
<th class="cen">FPslide</th>
<th class="left">PI</th>
<th class="left">PID</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
"""

NIGHT_FOOTER = """

<p> Along with others of more obvious meaning, the columns 'Read', 'xsll, ..',
and 'CDOP' specify the information needed to repeat the setup at the telescope
using 'hdriver'. 'Read' is the readout mode which can be one of 4 options,
namely 'FULL' for a full frame, 1-WIN or 2-WIN for standard windows mode, and
'DRIFT' for drift-mode. The 'xsll,xslr,..' column gives the 7 parameters that appear at
the bottom of the top-right window of hdriver (or 5 parameters in DRIFT mode).
In hdriver these are called xsll, xslr, xsul, xsur, ys, nx and ny (or xsl, xsr, ys, nx, ny
in drift mode), and there are up to two sets of them. In fullframe mode, these parameters
do not need to be specified. 'CDOP' are four Y/N flags which say whether or not clears, 
the dummy output, overscans and prescans were enabled or not. If overscans were enabled, extra
pixels are added in the Y direction. Prescans add extra pixels in X. 'Ncycle'
refers to how often each CCD is read out. The numbers correspond to nu, ng,
nr, ni, nz listed in hdriver.</p>

<p> 'Duty cycle' shows the fractional time collecting photons; the 'cadence'
is the (minimum) time between exposures, minimised over the different arms.
</p>

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
    font: 9pt sans-serif;
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

button {font-family: Arial,Helvetica,sans-serif;font-size: 10px;width: 10px; height: 10px; border-radius: 0%;}

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
    padding-right: 5px;
}

td.left {
    vertical-align: top;
    text-align: left;
    white-space: nowrap;
    padding-right: 5px;
}

td.lalert {
    color: #FFAAAA;
    vertical-align: top;
    text-align: left;
    white-space: nowrap;
    padding-right: 5px;
}

td.right {
    vertical-align: top;
    text-align: right;
    white-space: nowrap;
    padding-right: 5px;
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
    padding-right: 5px;
    font: 9pt sans-serif;
    white-space: normal;
}

td.bleft {
    color: #ffffa0;
    vertical-align: top;
    text-align: left;
    white-space: nowrap;
    padding-right: 5px;
    font: 9pt sans-serif;
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
    font-size: 9pt;
}

a:visited {
    color: #e0b0e0;
    text-decoration:underline;
    font-size: 9pt;
}

a:hover {
    color: red;
    text-decoration:underline;
    font-size: 9pt;
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

def gethead(head, key, default=''):
    """Extracts the value corresponding to key from head, returning default
    if the key is not present"""
    if key in head:
        return head[key]
    else:
        return default

def hlogger(args=None):
    """``hlogger server (dirnam)``

    Generates html logs for hipercam runs.

    hlogger expects to work on directories containing the runs for each
    night. It should be run from the top-level of the hipercam/logs directory.
    It extracts information from each run file which it writes to a file in
    JSON format and also as an html file. The JSON file acts as a short-cut
    for any re-run.

    Parameters:

        server   : bool
           True if the log is to be compiled from the server data. This option
           might be useful during observing. If server is True, the JSON and
           html files will be written in the current working directory.

        dirnam : string [if not server]
           Directory name to work on. This is run through 'glob' which gives
           UNIX-style name matching, so e.g. '2017*' will match 2017-10-17,
           2017-10-18 etc. The JSON and html files will be written into the
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

        # next are regular expressions to match run directories, nights, and
        # run files
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

                for nn, nname in enumerate(nnames):

                    print('  night {:s}'.format(nname))

                    if nn > 0:
                        bnight = os.path.basename(nnames[nn-1])
                        links = '\n<p><a href=../{0:s}/{0:s}.html>Previous night</a>'.format(bnight)

                    if nn < len(nnames)-1:
                        anight = os.path.basename(nnames[nn+1])
                        if nn > 0:
                            links += ', '
                        else:
                            links = '\n<p>'
                        links += '<a href=../{0:s}/{0:s}.html>Next night</a>\n</p>\n'.format(anight)
                    elif nn > 0:
                        links += '\n</p>\n'

                    # Write an entry for each night linking to the log for that night.
                    night = os.path.basename(nname)
                    fname = '{0:s}/{0:s}.html'.format(night)
                    ihtml.write('<tr><td><a href="{:s}">{:s}</a><td></tr>\n'.format(fname, night))

                    # Create the html file for the night
                    date = '{:s}, {:s}'.format(night, telescope)

                    with open(fname,'w') as nhtml:
                        # write header of night file
                        nhtml.write(NIGHT_HEADER1)
                        nhtml.write(NIGHT_HEADER2.format(date))
                        nhtml.write(links)
                        nhtml.write(TABLE_TOP)

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
                        for nrun, run in enumerate(runs):

                            if nrun % 20 == 0:
                                # write table header
                                nhtml.write(TABLE_HEADER)

                            # open the run file as an Rtime
                            rname = os.path.join(night,'data',run)
                            try:
                                rtime = hcam.hcam.Rtime(rname)
                            except:
                                exc_type, exc_value, exc_traceback = sys.exc_info()
                                traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
                                traceback.print_exc(file=sys.stdout)
                                print('Problem on run = ',rname)

                                # dummy info line just to allow us to proceed
                                nhtml.write('<tr>\n')
                                # run number
                                nhtml.write('<td class="lalert">{:s}</td>'.format(run[3:]))
                                nhtml.write('</tr>\n')
                                continue

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
                                '<td class="left">{:s}</td><td class="left">{:s}</td>'.format(
                                    ra, dec)
                            )


                            # timing info
                            ntotal = rtime.ntotal()
                            texps, toffs, nskips, tdead = rtime.tinfo()
                            # total = total time on target
                            # duty = worst duty cycle, percent
                            # tsamp = shortest sample time
                            ttotal = 0.
                            duty = 100
                            tsamp = 99000.
                            for texp, nskip in zip(texps, nskips):
                                ttotal = max(ttotal, (texp+tdead)*(ntotal // nskip))
                                duty = min(duty, 100.*texp/(texp+tdead))
                                tsamp = min(tsamp, texp+tdead)

                            # First & last timestamp
                            try:
                                tstart = rtime(1)[0].isot
                                tend = rtime(ntotal)[0].isot
                                nhtml.write(
                                    '<td class="cen">{:s}</td><td class="cen">{:s}</td>'.format(
                                        tstart[tstart.find('T')+1:tstart.rfind('.')],
                                        tend[tend.find('T')+1:tend.rfind('.')]
                                        )
                                    )
                            except:
                                exc_type, exc_value, exc_traceback = sys.exc_info()
                                traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
                                traceback.print_exc(file=sys.stdout)
                                print('Run =',rname)
                                nhtml.write('<td class="cen">--:--:--</td><td class="cen">--:--:--</td>')

                            # sample time
                            nhtml.write('<td class="right">{:.3f}</td>'.format(tsamp))

                            # duty cycle
                            nhtml.write('<td class="right">{:.1f}</td>'.format(duty))

                            # number of frames
                            nhtml.write('<td class="right">{:d}</td>'.format(ntotal))

                            # total exposure time
                            nhtml.write('<td class="right">{:d}</td>'.format(
                                int(round(ttotal)))
                            )

                            # run type
                            nhtml.write(
                                '<td class="left">{:s}</td>'.format(
                                    gethead(hd, 'IMAGETYP', '---'))
                            )

                            # readout mode
                            nhtml.write('<td class="cen">{:s}</td>'.format(
                                    TRANSLATE_MODE[rtime.mode])
                                        )

                            # cycle nums
                            nhtml.write('<td class="cen">{:s}</td>'.format(
                                ','.join([str(nskip) for nskip in nskips]))
                            )

                            # window formats
                            nhtml.write(
                                '<td class="cen">{:s}</td>'.format(rtime.wforms[0])
                                )
                            nhtml.write(
                                '<td class="cen">{:s}</td>'.format(
                                    rtime.wforms[1] if len(rtime.wforms) > 1 else '')
                                )

                            # binning
                            nhtml.write(
                                '<td class="cen">{:d}x{:d}</td>'.format(
                                    rtime.xbin,rtime.ybin)
                                )

                            # clear
                            nhtml.write(
                                '<td class="cen">{:s}</td>'.format(
                                    'On' if rtime.clear else 'Off')
                                )

                            # dummy
                            nhtml.write(
                                '<td class="cen">{:s}</td>'.format(
                                    'On' if rtime.dummy else 'Off')
                                )

                            # overscan/prescan
                            nhtml.write(
                                '<td class="cen">{:s},{:s}</td>'.format(
                                    'On' if rtime.oscan else 'Off',
                                    'On' if rtime.pscan else 'Off')
                                )

                            # CCD speed
                            if 'ESO DET SPEED' in hd:
                                speed = gethead(hd,'ESO DET SPEED')
                                if speed == 0:
                                    speed = 'Slow'
                                elif speed == 1:
                                    speed = 'Fast'
                                else:
                                    speed = '--'
                            else:
                                speed = '--'
                            nhtml.write(
                                '<td class="cen">{:s}</td>'.format(speed)
                                )

                            # Fast clocks
                            if 'ESO DET FASTCLK' in hd:
                                fclock = gethead(hd,'ESO DET FASTCLK')
                                if fclock == 0:
                                    fclock = 'No'
                                elif fclock == 1:
                                    fclock = 'Yes'
                                else:
                                    fclock = '--'
                            else:
                                fclock = '--'
                            nhtml.write(
                                '<td class="cen">{!s}</td>'.format(fclock)
                                )

                            # Tbytes problem
                            nhtml.write(
                                '<td class="cen">{:s}</td>'.format(
                                    'OK' if rtime.ntbytes == 36 else 'BAD'
                                    )
                                )

                            # Focal plane slide
                            nhtml.write(
                                '<td class="cen">{!s}</td>'.format(
                                    gethead(hd,'FPslide','----')
                                    )
                                )

                            # PI
                            nhtml.write(
                                '<td class="left">{:s}</td>'.format(
                                    gethead(hd, 'PI', '---'))
                            )

                            # Program ID
                            nhtml.write(
                                '<td class="left">{:s}</td>'.format(
                                    gethead(hd, 'PROGRM', '---'))
                            )

                            # run number again
                            nhtml.write(
                                '<td class="left">{:s}</td>'.format(run[3:])
                            )

                            # comments
                            pcomm = gethead(hd, 'RUNCOM', '').strip()
                            if pcomm == 'None':
                                pcomm = ''
                            if pcomm != '':
                                if not pcomm.endswith('.'):
                                    pcomm += '. '
                                else:
                                    pcomm += '.'
                            nhtml.write(
                                '<td class="left">{:s}{:s}</td>'.format(pcomm,hlog[run])
                            )

                            # end the row
                            nhtml.write('\n</tr>\n')

                        # finish off the night file
                        nhtml.write('</table>\n{:s}'.format(links))
                        nhtml.write(NIGHT_FOOTER)

                # finish off the run table
                ihtml.write('\n</table>\n\n')

            # finish the main index
            ihtml.write(INTRODUCTION_FOOTER)


    # write out the css file
    with open('hiper.css','w') as fout:
        fout.write(CSS)
