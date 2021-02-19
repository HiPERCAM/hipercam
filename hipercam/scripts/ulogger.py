import sys
import traceback
import os
import time
import glob
import re
import math
import warnings

import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta, TimeISO
from astropy.coordinates import get_sun, get_moon, EarthLocation, SkyCoord, AltAz
import astropy.units as u

import hipercam as hcam
from hipercam.utils import format_ulogger_table

__all__ = [
    "ulogger",
]

############################################################################
#
# ulogger -- generates html logs and spreadsheets for ultracam and ultraspec
#
############################################################################

###########
#
# Some constant stuff
#

# Header of the main file

INDEX_HEADER = """<html>
<head>

<title>{instrument} data logs</title>
</head>

<body>
<h1>{instrument} data logs</h1>

<p> These on-line logs summarise all {instrument} runs. They were
automatically generated from the data files and hand-written logs.

<p>For a searchable file summarising the same information, along with
calculated data on the Sun and Moon for each run, see this <a
href="{linstrument}-log.xlsx">spreadsheet</a>. This also contains a column
'Nlink' with hyperlinks to the on-line log of whichever night contains
a given run. In addition after the comment column the spreadsheet has
a number of angles (altitudes and azimuths etc, all in degrees) and
the lunar phase. 'Start' and 'end' refer to the mid-time of the first
and last frame of a run, while 'mid' is half-way between them. These
parameters could allow one to select dark runs for instance. If there is
anything else that you think could be added to this file let me know.

"""

# Footer of the main file

INDEX_FOOTER = """
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
var hcols = [4, 9];

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

<link rel="stylesheet" type="text/css" href="ultra.css" />
"""

NIGHT_HEADER2 = """<title>{instrument} log {date}</title>
</head>

<body>

<h1>{instrument} log {date}</h1>

<p>
The table below lists information on runs from the night starting on {date}.
See the end of the table for some more details on the meanings of the various columns.
You can hide any column using the buttons in the top line and "shrink" will reduce clutter
to focus on those of most interest.
"""

TABLE_TOP = """<p>
<button style="width: 120px; height: 30px;" onclick="shrink();">Shrink table</button>
<button style="width: 140px; height: 30px;" onclick="restore();">Restore columns</button>


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
<td><button id="hidden26" onclick="hide(26)"></button></td>
<td><button id="hidden27" onclick="hide(27)"></button></td>
<td align="left"><button id="hidden28" onclick="hide(28)"></button></td>
</tr>

<tr>
<th class="left">Run<br>no.</th>
<th class="left">Target<br>name</th>
<th class="left">Auto<br>ID</th>
<th class="left">RA (J2000)</th>
<th class="left">Dec&nbsp;(J2000)</th>
<th class="cen">Date<br>(start)</th>
<th class="cen">Start</th>
<th class="cen">End</th>
<th class="right">Total<br>(sec)</th>
<th class="right">Cadence<br>(sec)</th>
<th class="right">Exposure<br>(sec)</th>
<th class="right">Nframe</th>
<th class="right">Nok</th>
<th class="cen">Filters</th>
<th class="left">Type</th>
<th class="cen">Read<br>mode</th>
<th class="left">Nb</th>
<th class="cen">
xl1,xr1,ys1,nx1,ny1<br>
xl2,xr2,ys2,nx2,ny2<br>
xl3,xr3,ys3,nx3,ny3</th>
<th class="cen">XxY<br>bin</th>
<th class="cen">Clr</th>
<th class="cen">Read<br>speed</th>
<th class="cen">FPslide</th>
<th class="cen">Observers</th>
<th class="left">PI</th>
<th class="left">PID</th>
<th class="left">Tel</th>
<th class="cen">Size<br>(MB)</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
"""

NIGHT_FOOTER = """

<p> 'Instr. PA' is the instrumental PA. 'Clr' indicates whether clears were enabled;
'Read mode' is the readout mode which can be one of several options:
'FFCLR' for full frames with clear; 'FFNCLR' full frames with no clear;
'1-PAIR', '2-PAIR', '3-PAIR', for standard windowed modes, etc. The 'xl1,xr1,..' column gives the 5
parameters defining the window pair (ULTRACAM). In fullframe mode,
these parameters do not need to be specified. The 'cadence' is the
time between exposures, the 'exposure' the actual time spent integrating, although
note that these numbers are not always entirely accurate.
</p>

<address>Tom Marsh, Warwick</address>
</body>
</html>
"""

MONTHS = {
    "01": "January",
    "02": "February",
    "03": "March",
    "04": "April",
    "05": "May",
    "06": "June",
    "07": "July",
    "08": "August",
    "09": "September",
    "10": "October",
    "11": "November",
    "12": "December",
}

# CSS to define the look of the log web pages
# This is written to 'ultra.css' and referred to in the
# html pages
CSS = """

body{
    background-color: #000000;
    font: 10pt sans-serif;
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
    font: 10pt sans-serif;
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
    font-size: 10pt;
}

a:visited {
    color: #e0b0e0;
    text-decoration:underline;
    font-size: 10pt;
}

a:hover {
    color: red;
    text-decoration:underline;
    font-size: 10pt;
}
"""

# end of CSS

def correct_ra_dec(ra, dec):
    """Fixes up RAs and Decs"""
    if ra == "UNDEF":
        ra = ""
    else:
        rah, ram, ras = ra.split(":")
        rah, ram, ras = int(rah), int(ram), float(ras)
        ra = "{:02d}:{:02d}:{:05.2f}".format(rah, ram, ras)

    if dec == "UNDEF":
        dec = ""
    else:
        decd, decm, decs = dec.split(":")
        if decd.find("-") > -1:
            negative = True
        else:
            negative = False
        decd, decm, decs = int(decd), int(decm), float(decs)
        decd = -abs(decd) if negative else abs(decd)
        dec = "{:+03d}:{:02d}:{:04.1f}".format(decd, decm, decs)

    return (ra, dec)

def ulogger(args=None):
    """``ulogger``

    Generates html logs for ULTRACAM and ULTRASPEC runs, and also generates
    a searchable spreadsheet.

    ulogger expects to work on directories containing the runs for
    each night. It should be run from the ultra(cam|spec)/raw_data
    directory.  It extracts information from each run file. It runs
    without any parameters.  ulogger also writes out a spreadsheet
    with useful info line-by-line for each run.

    If run at Warwick, it writes data to the web pages. Otherwise it
    writes them to a sub-directory of raw_data called "logs".

    Bit of a specialist routine this one; if you have access to the
    on-line logs at Warwick, it should be unnecessary. It requires the
    installation of the python module xlsxwriter in order to write an
    excel spreadsheet; this is not a required module for the whole
    pipeline and may not exist. The script will fail early on if it is
    not installed.

    """
    from astroplan import moon_phase_angle
    warnings.filterwarnings('ignore')

    barr = []
    cwd = os.getcwd()
    if os.path.basename(cwd) != "raw_data":
        print("** ulogger must be run in a directory called 'raw_data'")
        print("ulogger aborted", file=sys.stderr)
        return

    if cwd.find("ultracam") > -1:
        instrument = "ULTRACAM"
    elif cwd.find("ultraspec") > -1:
        instrument = "ULTRASPEC"
    else:
        print("** ulogger: cannot find either ultracam or ultraspec in path")
        print("ulogger aborted", file=sys.stderr)
        return
    linstrument = instrument.lower()

    # location to write files
    if os.path.exists(f'/storage/astro2/www/phsaap/{linstrument}/logs'):
        root = f'/storage/astro2/www/phsaap/{linstrument}/logs/logs'
    else:
        root = 'logs'

    # Try an import now to get an early warning of problems
    import xlsxwriter

    # next are regular expressions to match run directories, nights, and
    # run files
    rre = re.compile("^\d\d\d\d-(\d\d|P\d\d\d)$")
    nre = re.compile("^\d\d\d\d-\d\d-\d\d$")
    fre = re.compile("^run\d\d\d\.xml$")

    # Get list of run directories
    rnames = [
        rname
        for rname in os.listdir(".")
        if rre.match(rname)
        and os.path.isdir(rname)
        and os.path.isfile(os.path.join(rname, "telescope"))
    ]
    rnames.sort()

    if len(rnames) == 0:
        print(
            "there were no run directories of the form "
            "YYYY-MM with a file called 'telescope' in them",
            file=sys.stderr,
        )
        print("ulogger aborted", file=sys.stderr)
        return

    # Get started. First make sure all files are created with the
    # right permissions
    os.umask(0o022)

    # Ensure the root directory exists.
    os.makedirs(root, exist_ok=True)

    print(f'Will write to directory = "{root}".')

    # Load target positional information
    targets = Targets('TARGETS')
    auto_targets = Targets('AUTO_TARGETS')

    # Load targets to skip and failed targets
    skip_targets, failed_targets = load_skip_fail()

    # Index file
    index = os.path.join(root, 'index.html')
    with open(index, "w") as ihtml:
        # start the top level index html file
        ihtml.write(
            INDEX_HEADER.format(
                instrument=instrument,
                linstrument=linstrument
            )
        )

        for rname in rnames:
            print(f"\nProcessing run {rname}")

            # write in run date, start table of nights
            rn = os.path.basename(rname)
            year, month = rn.split("-")
            with open(os.path.join(rname, "telescope")) as tel:
                telescope = tel.read().strip()
            ihtml.write(
                f"<h2>{MONTHS[month]} {year}, {telescope}</h2>\n"
            )
            ihtml.write("\n<p>\n")

            # set site
            if telescope == 'WHT':
                observatory = EarthLocation.of_site('Roque de los Muchachos')
            elif telescope == 'VLT':
                observatory = EarthLocation.of_site('Cerro Paranal')
            elif telescope == 'NTT':
                observatory = EarthLocation.of_site('Cerro La Silla')
            else:
                raise ValueError('did not recognise telescope =',telescope)

            # get night directories
            nnames = [
                os.path.join(rname, ndir)
                for ndir in os.listdir(rname)
                if nre.match(ndir)
                and os.path.isdir(os.path.join(rname, ndir))
            ]
            nnames.sort()

            if len(nnames) == 0:
                print(
                    "found no night directories of the form YYYY-MM-DD"
                    f" in {rname}", file=sys.stderr,
                )
                print("ulogger aborted", file=sys.stderr)
                return

            for nn, nname in enumerate(nnames):

                print(f"  night {nname}")

                # create directory for any meta info such as the times
                meta = os.path.join(nname, 'meta')
                os.makedirs(meta, exist_ok=True)

                links = '\n<p><a href="index.html">Run index</a>'
                if nn > 0:
                    bnight = os.path.basename(nnames[nn - 1])
                    links += f', <a href="{bnight}.html">Previous night</a>'
                else:
                    links += f', Previous night'

                if nn < len(nnames) - 1:
                    anight = os.path.basename(nnames[nn + 1])
                    links += f', <a href="{anight}.html">Next night</a>\n</p>\n'
                else:
                    links += f', Next night\n</p>\n'

                links += "\n</p>\n"

                # Write an entry for each night linking to the log for that night.
                night = os.path.basename(nname)
                ihtml.write(
                    f'<a href="{night}.html">{night}</a> \n'
                )

                # Create the html file for the night
                date = f"{night}, {telescope}"
                fname = os.path.join(root, f"{night}.html")

                with open(fname, "w") as nhtml:

                    # write header of night file
                    nhtml.write(NIGHT_HEADER1)
                    nhtml.write(NIGHT_HEADER2.format(date=date,instrument=instrument))
                    nhtml.write(links)
                    nhtml.write(TABLE_TOP)

                    # read and store the hand written log
                    handlog = os.path.join(night, f"{night}.dat")
                    hlog = Log(handlog)

                    # load all the run names
                    runs = [run[:-4] for run in os.listdir(night) if fre.match(run)]
                    runs.sort()

                    ####################################
                    #
                    # Get or create timing info

                    times = os.path.join(meta, 'times')
                    tdata = {}
                    if os.path.exists(times):
                        # pre-existing file found
                        with open(times) as tin:
                            for line in tin:
                                arr = line.split()
                                tdata[arr[0]] = [
                                    '' if val == 'UNDEF' else val for val in arr[1:]
                                ]
                        print('Read timing data from',times)

                    else:
                        # need to generate timing data, which can take a while so
                        # we store the results to a disk file for fast lookup later.
                        with open(times,'w') as tout:
                            for run in runs:
                                dfile = os.path.join(night, run)
                                try:
                                    rtime = hcam.ucam.Rtime(dfile)

                                    # Find first good time
                                    for n, tdat in enumerate(rtime):
                                        time, tinfo = tdat[:2]
                                        if time.good:
                                            mjd_start = time.mjd
                                            ts = Time(mjd_start, format="mjd", precision=2)
                                            ut_start = ts.hms_custom
                                            expose = round(time.expose,3)
                                            n_start = n+1
                                            break
                                    else:
                                        raise hcam.HipercamError(f'No good time found in {dfile}')

                                    # Find last good time. First we just go for times near the
                                    # end of the run. Failing that, we try again from the start,
                                    # to account for runs with time stamp issues.
                                    if rtime.header['MODE'] == 'DRIFT':
                                        # ultracam
                                        win = rhead.win[0]
                                        nyu = win.ny*rhead.ybin
                                        nback = int((1033/nyu + 1) / 2) + 3
                                    elif rtime.header['MODE'] == 'UDRIFT':
                                        # ultraspec
                                        win = rhead.win[0]
                                        nyu = win.ny*rhead.ybin
                                        nback = int((1037/nyu + 1) / 2) + 3
                                    else:
                                        # non drift mode
                                        nback = 4

                                    nbytes = os.stat(dfile + '.dat').st_size
                                    ntotal = nbytes // rtime.framesize
                                    nreset = max(1, ntotal - nback)
                                    rtime.set(nreset)

                                    flast = False
                                    for n, tdat in enumerate(rtime):
                                        time, tinfo = tdat[:2]
                                        if time.good:
                                            mjd_end = time.mjd
                                            ts = Time(mjd_start, format="mjd", precision=2)
                                            ut_end = ts.hms_custom
                                            n_end = nreset + n
                                            expose = max(expose, round(time.expose,3))
                                            flast = True

                                    if not flast:
                                        # no time found near end,
                                        # grind it out by going
                                        # through the whole run
                                        rtime.set()
                                        for n, tdat in enumerate(rtime):
                                            time, tinfo = tdat[:2]
                                            if time.good:
                                                mjd_end = time.mjd
                                                ts = Time(mjd_start, format="mjd", precision=2)
                                                ut_end = ts.hms_custom
                                                n_end = n + 1
                                                expose = max(expose, round(time.expose,3))

                                    nok = n_end-n_start+1
                                    if n_end > n_start:
                                        cadence = round(86400*(mjd_end-mjd_start)/(n_end-n_start),3)
                                        tdata[run] = [ut_start,mjd_start,ut_end,mjd_end,cadence,expose,nok,ntotal]
                                    else:
                                        cadence = 'UNDEF'
                                        tdata[run] = [ut_start,mjd_start,ut_end,mjd_end,'',expose,nok,ntotal]
                                    tout.write(f'{run} {ut_start} {mjd_start} {ut_end} {mjd_end} {cadence} {expose} {nok} {ntotal}\n')

                                except hcam.ucam.PowerOnOffError:
                                    # Power on/off
                                    tdata[run] = ['power-on-off',]
                                    tout.write(f'{run} power-on-off\n')

                                except hcam.HipercamError:
                                    # No good times
                                    tdata[run] = ['','','','','','',0,ntotal]
                                    tout.write(f'{run} UNDEF UNDEF UNDEF UNDEF UNDEF UNDEF 0 {ntotal}\n')

                                except:
                                    # some other failure
                                    exc_type, exc_value, exc_traceback = sys.exc_info()
                                    traceback.print_tb(
                                        exc_traceback, limit=1, file=sys.stderr
                                    )
                                    traceback.print_exc(file=sys.stderr)
                                    print("Problem on run = ", dfile)

                                    # Load of undefined
                                    tdata[run] = 8*['']
                                    tout.write(f'{run} {" ".join(8*["UNDEF"])}\n')

                        print('Written timing data to',times)


                    ###############################
                    # Get or create positional info

                    posdata = os.path.join(meta, 'posdata')
                    pdata = {}
                    if os.path.exists(posdata):
                        # pre-existing file found
                        with open(posdata) as pin:
                            for line in pin:
                                arr = line.split()
                                arr[3] = arr[3].replace(' ','~')
                                pdata[arr[0]] = [
                                    val if val != 'UNDEF' else '' for val in arr[1:]
                                ]
                        print('Read position data from',posdata)

                    else:
                        # need to generate sun / moon / airmass data
                        # which takes a while so store the results
                        # to a disk file for fast lookup in the future. Only
                        # bother with runs for which we can find a position
                        with open(posdata,'w') as pout:
                            for run in runs:

                                if len(tdata[run]) == 1:
                                    # means its a power on/off
                                    continue

                                # open the run file as an R
                                runname = os.path.join(night, run)
                                try:
                                    rhead = hcam.ucam.Rhead(runname)
                                except:
                                    continue

                                # object name
                                if hlog.format == 1:
                                    target = hlog.target[run]
                                else:
                                    target = rhead.header.get("TARGET",'')

                                # RA, Dec lookup
                                if target == '':
                                    ra, dec, autoid = 'UNDEF', 'UNDEF', 'UNDEF'
                                elif target in skip_targets:
                                    ra, dec, autoid = 'UNDEF', 'UNDEF', 'UNDEF'
                                elif target in targets.lnames:
                                    autoid, dct = targets(target)
                                    ra, dec = dct['ra'], dct['dec']
                                elif target in auto_targets.lnames:
                                    autoid, dct = auto_targets(target)
                                    ra, dec = dct['ra'], dct['dec']
                                else:
                                    # attempt simbad lookup here
                                    ra, dec, autoid = 'UNDEF', 'UNDEF', 'UNDEF'

                                # start accumulating stuff to write out
                                arr = [ra, dec, autoid]

                                # time-dependent info
                                ut_start, mjd_start, ut_end, mjd_end, cadence, expose, nok, ntotal = tdata[run]
                                try:

                                    mjd_start = float(mjd_start)
                                    mjd_end = float(mjd_end)
                                    tstart = Time(mjd_start, format='mjd')
                                    tmid = Time((mjd_start+mjd_end)/2, format='mjd')
                                    tend = Time(mjd_end, format='mjd')

                                    if ra != 'UNDEF' and dec != 'UNDEF':
                                        # Calculate the Alt, Az at
                                        # start, middle, end
                                        frames = AltAz(obstime=[tstart,tmid,tend], location=observatory)
                                        points = SkyCoord(f'{ra} {dec}',unit=(u.hourangle, u.deg)).transform_to(frames)
                                        alts = [round(alt,1) for alt in points.alt.degree]
                                        azs = [round(az,1) for az in points.az.degree]
                                        arr += alts + azs

                                        # Now calculate the angular
                                        # distance from the Sun and
                                        # Moon at the mid-time
                                        sun_mid = get_sun(tmid).transform_to(frames[1])
                                        moon_mid = get_moon(tmid).transform_to(frames[1])
                                        point_mid = points[1]
                                        sun_dist = point_mid.separation(sun_mid).degree
                                        moon_dist = point_mid.separation(moon_mid).degree
                                        arr += [round(sun_dist,1),round(moon_dist,1)]

                                    else:
                                        arr += 8*['UNDEF']

                                    # Now some data on the altitude of the Sun & Moon
                                    frame = AltAz(obstime=tstart, location=observatory)
                                    sun_start = get_sun(tstart).transform_to(frame)
                                    moon_start = get_moon(tstart).transform_to(frame)

                                    # end
                                    frame = AltAz(obstime=tend, location=observatory)
                                    sun_end = get_sun(tend).transform_to(frame)
                                    moon_end = get_moon(tend).transform_to(frame)

                                    # Lunar phase at the mid-point only.
                                    moon_phase = moon_phase_angle(tmid).value / np.pi

                                    arr += [
                                        round(sun_start.alt.degree,1), round(sun_end.alt.degree,1),
                                        round(moon_start.alt.degree,1), round(moon_end.alt.degree,1),
                                        round(moon_phase,2),
                                    ]

                                except:
                                    # write out info
                                    arr += 13*['UNDEF']

                                autoid_nospace = arr[2].replace(' ','~')
                                pout.write(
                                    f'{run} {arr[0]} {arr[1]} {autoid_nospace} {arr[3]} ' +
                                    f'{arr[4]} {arr[5]} {arr[6]} {arr[7]} {arr[8]} {arr[9]} ' +
                                    f'{arr[10]} {arr[11]} {arr[12]} {arr[13]} {arr[14]} {arr[15]}\n'
                                )

                                pdata[run] = [
                                    '' if val == 'UNDEF' else val for val in arr
                                ]

                        print('Written positional data to',posdata)

                    # Right, finally!
                    #
                    # Have now either read or created (and dumped)
                    # basic timing and positional data for all runs
                    # for this night.  Now wind through the runs
                    # getting basic info and writing a row of info to
                    # the html file for the night in question, and accccumulating
                    # data for the spreadsheet

                    for nrun, run in enumerate(runs):

                        if nrun % 20 == 0:
                            # write table header
                            nhtml.write(TABLE_HEADER)

                        if len(tdata[run]) == 1:
                            # means its a power on/off
                            continue

                        # open the run file as an Rhead
                        runname = os.path.join(night, run)
                        try:
                            rhead = hcam.ucam.Rhead(runname)
                        except:
                            exc_type, exc_value, exc_traceback = sys.exc_info()
                            traceback.print_tb(
                                exc_traceback, limit=1, file=sys.stderr
                            )
                            traceback.print_exc(file=sys.stderr)
                            print("Problem on run = ", runname)

                            # dummy info line just to allow us to proceed
                            nhtml.write("<tr>\n")
                            # run number
                            nhtml.write(f'<td class="lalert">{run[3:]}</td>')
                            nhtml.write("</tr>\n")
                            brow = [night, run[3:]] + 45*[None]
                            continue

                        hd = rhead.header

                        # start the row
                        nhtml.write("<tr>\n")

                        # run number
                        runno = run[3:]
                        nhtml.write('<td class="left">{:s}</td>'.format(runno))

                        # start list to append to array for spreadsheet
                        brow = [night, runno,]

                        # object name
                        if hlog.format == 1:
                            target = hlog.target[run]
                        else:
                            target = hd.get("TARGET",'')

                        # position
                        pdat = pdata[run]
                        ra, dec, autoid = pdat[:3]

                        try:
                            ra, dec = float(ra), float(dec)
                            rastr = dec2sexg(ra, False, 2)
                            decstr = dec2sexg(dec, True, 1)
                            ra *= 15
                            ra = round(ra,5)
                            dec = round(dec,4)
                        except:
                            rastr, decstr = '', ''

                        nhtml.write(f'<td class="left">{target}</td><td class="left">{autoid}</td>')
                        nhtml.write(f'<td class="left">{rastr}</td><td class="left">{decstr}</td>')
                        brow += [target,autoid, rastr, decstr, ra, dec, rname]

                        # Timing info
                        ut_start, mjd_start, ut_end, mjd_end, cadence, expose, nok, ntotal = tdata[run]

                        try:
                            # Start time, date
                            mjd_start = float(mjd_start)
                            tstart = Time(mjd_start, format='mjd').isot
                            date_start = tstart[:tstart.find("T")]
                            nhtml.write(
                                f'<td class="cen">{date_start}</td> <td class="cen">{ut_start}</td>'
                            )
                            brow += [date_start, ut_start]
                        except:
                            nhtml.write(2*'<td></td>')
                            brow += 2*['']

                        try:
                            # End time, total time, cadence (assumes we have a valid start time)
                            mjd_end = float(mjd_end)
                            ttime = round(86400.*(mjd_end - mjd_start))
                            nhtml.write(
                                f'<td class="cen">{ut_end}</td> <td class="right">{ttime}</td> <td class="right">{cadence}</td>'
                            )
                            brow += [ut_end, ttime, cadence]
                        except:
                            nhtml.write(3*'<td></td>')
                            brow += 3*['']

                        # Actual exposure time per frame
                        nhtml.write(f'<td class="right">{expose}</td>')
                        brow.append(expose)

                        # number of frames, number ok
                        nhtml.write(f'<td class="right">{ntotal}</td><td class="right">{nok}</td>')
                        brow += [ntotal,nok]

                        # filters used
                        if hlog.format == 1:
                            filters = hlog.filters.get(run, '')
                        else:
                            filters = hd.get("filters", '')
                        nhtml.write(f'<td class="cen">{filters}</td>')
                        brow.append(None if filters == '' else filters)

                        # run type
                        itype = hd.get("DTYPE", '')
                        nhtml.write(f'<td class="left">{itype}</td>')
                        brow.append(None if itype == '' else itype)

                        # readout mode
                        nhtml.write(f'<td class="cen">{rhead.mode}</td>')
                        brow.append(rhead.mode)

                        # nblue
                        nblue = hd['NBLUE']
                        nhtml.write(f'<td class="cen">{nblue}</td>')
                        brow.append(nblue)

                        # window formats
                        win1 = rhead.wforms[0]
                        nhtml.write(f'<td class="cen">{win1}')
                        win1 = win1.replace('&nbsp;',' ')
                        brow.append(win1)

                        win2 = rhead.wforms[1] if len(rhead.wforms) > 1 else ""
                        if win2 != '':
                            nhtml.write(f'<br>{win2}')
                        win2 = win2.replace('&nbsp;',' ')
                        brow.append(win2)

                        win3 = rhead.wforms[1] if len(rhead.wforms) > 2 else ""
                        if win3 != '':
                            nhtml.write(f'<br>{win3}</td>')
                        else:
                            nhtml.write('</td>')
                        win3 = win3.replace('&nbsp;',' ')
                        brow.append(win3)

                        # binning
                        binning = f'{rhead.xbin}x{rhead.ybin}'
                        nhtml.write(f'<td class="cen">{binning}</td>')
                        brow.append(binning)

                        # clear
                        clear = "Y" if hd['CLEAR'] else "N"
                        nhtml.write(f'<td class="cen">{clear}</td>')
                        brow.append(clear)

                        # CCD speed
                        speed = hd['GAINSPED']
                        nhtml.write(f'<td class="cen">{speed}</td>')
                        brow.append(speed)

                        # Focal plane slide
                        fpslide = hd.get('SLIDEPOS','UNKNOWN')
                        try:
                            fpslide = float(fpslide)
                            fpslide = round(fpslide,1)
                        except:
                            fpslide = ''
                        nhtml.write(f'<td class="cen">{fpslide}</td>')
                        brow.append(None if fpslide == '' else fpslide)

                        if instrument == 'ULTRASPEC':
                            # instr PA
                            instpa = hd.get("INSTRPA", "----")
                            nhtml.write(f'<td class="cen">{instpa}</td>')
                            brow.append(instpa)

                        # Observers
                        observers = hd.get("OBSERVRS", '')
                        nhtml.write(f'<td class="cen">{observers}</td>')
                        brow.append(None if observers == '' else observers)

                        # PI
                        pi = hd.get("PI", '')
                        nhtml.write(f'<td class="left">{pi}</td>')
                        brow.append(None if pi == '' else pi)

                        # Program ID
                        pid = hd.get("ID", "")
                        nhtml.write(f'<td class="left">{pid}</td>')
                        brow.append(None if pid == '' else pid)

                        # Telescope name
                        nhtml.write(f'<td class="left">{telescope}</td>')
                        brow.append(telescope)

                        try:
                            # file size
                            mbytes = round(os.stat(runname + '.dat').st_size / (1024*1024),1)
                            nhtml.write(f'<td class="cen">{mbytes}</td>')
                            brow.append(mbytes)
                        except:
                            nhtml.write(f'<td class="cen"></td>')
                            brow.append(None)

                        # run number again. Add null to the
                        # spreadsheet to allow a formula
                        nhtml.write(f'<td class="left">{runno}</td>')
                        brow.append('NULL')

                        # comments
                        pcomm = hd.get("RUNCOM", "").strip()
                        if pcomm == "None" or pcomm == "UNDEF":
                            pcomm = ""
                        if pcomm != "":
                            if not pcomm.endswith("."):
                                pcomm += " "
                            else:
                                pcomm += ". "

                        lcomm = hlog.comment.get(run,'No comment in log')
                        comments = f'{pcomm}{lcomm}'
                        nhtml.write(f'<td class="left">{comments}</td>')
                        brow.append(comments)

                        # Finally tack on extra positional stuff for the spreadsheet only
                        brow += pdat[3:]

                        # at last: end the row
                        nhtml.write("\n</tr>\n")

                        barr.append(brow)

                    # finish off the night file
                    nhtml.write("</table>\n{:s}".format(links))
                    nhtml.write(NIGHT_FOOTER)

        # finish the main index
        ihtml.write(INDEX_FOOTER)

    # write out the css file
    css = os.path.join(root, "ultra.css")
    with open(css, "w") as fout:
        fout.write(CSS)

    # Create and write out spreadsheet
    colnames = (
        'Night', 'Run no.', 'Target name', 'RA (J2000)', 'Dec (J2000)',
        'RA (deg)', 'Dec (deg)', 'Obs run', 'Date\n(start)',
        'UTC\nStart', 'UTC\nEnd', 'Total\n(sec)', 'Cadence\n(sec)',
        'Exposure\n(sec)', 'Nframe', 'Nok', 'Filters', 'Run\ntype',
        'Read\nmode', 'Nb', 'Win1', 'Win2' , 'Win3', 'XxY\nbin',
        'Clr', 'Read\nspeed', 'FPslide', 'Observers', 'PI', 'PID',
        'Tel', 'Size\n(MB)', 'Nlink', 'Comment', 'Alt\nstart',
        'Alt\nmiddle', 'Alt\nend', 'Az\nstart', 'Az\nmiddle',
        'Az\nend', 'Sun sep\nmid', 'Moon sep\nmid', 'Sun alt\nstart',
        'Sun alt\nend', 'Moon alt\nstart', 'Moon alt\nend',
        'Moon\nphase',
    )

    ptable = pd.DataFrame(barr, columns=colnames)

    spreadsheet = os.path.join(root, f"{linstrument}-log.xlsx")

    format_ulogger_table(spreadsheet, ptable, linstrument)

    print(f'\nAll done. Look in {root} for the outputs.')



class Log(object):
    """
    Class to read and store log file data. These come in two formats:

    1) Old style: run, target name, filters, comment
    2) New style: run, comment (target names are in the xml files)

    The class just stores the data in a couple of dictionaries
    'comment' and 'target'; 'format' is an integer specifying
    the format as above. 'target' is blank in the case of format == 2.
    """

    def __init__(self, fname):
        """
        Constructs a new Log given a file. Makes empty
        dictionaries if none found and reports an error
        """
        self.format  = 2
        self.target  = {}
        self.filters = {}
        self.comment = {}

        try:
            rec    = re.compile('file\s+object\s+filter', re.I)
            old    = re.compile('\s*(\S+)\s+(\S+)\s+(.*)$')
            oldii  = re.compile('\s*(\S+)\s*$')
            with open(fname) as f:
                for line in f:
                    m = rec.search(line)
                    if m:
                        self.format = 1
                        if len(self.comment):
                            raise Exception('Error in night log = ' + fname + ', line = ' + line)

                    if line.startswith('run'):
                        run = line[:6]
                        if self.format == 2:
                            self.comment[run] = line[6:].strip()
                        else:
                            m = old.search(line[6:])
                            if m:
                                self.target[run]  = m.group(1)
                                self.filters[run] = m.group(2)
                                self.comment[run] = m.group(3)
                            else:
                                m = oldii.search(line[6:])
                                if m:
                                    self.target[run]  = m.group(1)

        except Exception as err:
           sys.stderr.write('Night log problem:' + str(err) + '\n')


class Targets(dict):
    """Class to read and store the target positions and regular expressions.

    It is a dictionary keyed on target names. Each item is in turn a
    dictionary with the following keys:

    'ra'     -- RA (decimal hours)
    'dec'    -- Dec (decimal degrees)
    'names'  -- A list of matching names. These are names from the logs with
                their typos etc that will be identified with the object.

    """

    def __init__(self, *fnames):
        """Constructs a new Targets object. Makes empty dictionaries if none
        found and reports an error. Multiple file names can be
        specified; any which do not exist will be skipped.  The files
        must have the format

        name hh:mm:ss.ss [+-]dd:mm:ss.s lname1 lname2 .. lnameN

        where name is the name that will be attached to the target,
        and lname1, lname2 etc, are the names in the logs that will be
        matched to give name. The set of 'name's must be unique as
        must the set of 'lname's

        The RA and Dec are assumed to be ICRS. Any '~' in either name
        or lname will be replaced by single spaces.

        An attribute called lnames is maintained which is a dictionary
        keyed on the lnames and translating to the name
        """
        self.lnames = {}
        for fname in fnames:
            if os.path.isfile(fname):
                with open(fname) as f:
                    for line in f:
                        if not line.startswith('#') and not line.isspace():
                            tokens = line.strip().split()
                            if len(tokens) < 4:
                                # must have at least one log name
                                raise Exception('Targets: invalid target, file = ' + fname + ', line = ' + line)

                            # prepare data
                            target,ra,dec = tokens[:3]
                            target = target.replace('~',' ')
                            ra,dec,system = str2radec(ra + ' ' + dec)
                            names = [token.replace('~',' ') for token in tokens[3:]]

                            # check that, if the target has been
                            # entered before, as is possible, that it
                            # is self-consistent
                            if target in self:
                                entry = self[target]
                                if entry['ra'] != ra or entry['dec'] != dec:
                                    raise Exception(
                                        'Targets: file = ' + fname + ', line = ' + line + \
                                        '\nTarget =' + target + ' already has an entry but with a different position.'
                                    )

                            # add names to the dictionary maintained to check for uniqueness
                            for name in names:
                                if name in self.lnames:
                                    raise Exception(
                                        'Targets: file = ' + fname + ', line = ' + line + \
                                        '\nName = ' + name + ' already exists.'
                                    )
                                self.lnames[name] = target

                            self[target] = {'ra' : ra, 'dec' : dec, 'names' : names}

                print(len(self),'targets after loading',fname)
            else:
                print('No targets loaded from',fname,'as it does not exist.')

    def __call__(self, name):
        """
        Call with name to lookup. Returns standardised name plus
        dictionary with 'ra', 'dec' and 'names'
        """
        target = self.lnames[name]
        return target, self[target]

    def write(self, fname):
        """
        Write targets out to disk file fname.
        """

        # write in RA order
        ras = dict([(targ,entry['ra']) for targ, entry in self.items()])
        targs = sorted(ras, key=ras.get)

        with open(fname,'w') as f:
            f.write("""
#
# File of targets written by ulogger.Targets.write
#

""")

            for targ in targs:
                entry = self[targ]
                pos = dec2sexg(entry['ra'],False,2) + '   ' + dec2sexg(entry['dec'],True,1)
                f.write('%-32s %s' % (targ.replace(' ','~'),pos))
                lnames = ' '.join([name.replace(' ','~') for name in entry['names']])
                f.write(' ' + lnames + '\n')

    def tohtml(self, fname):
        """
        Write targets out to an html file (give full name)
        """

        # write in RA order
        ras   = dict([(targ,entry['ra']) for targ, entry in self.items()])
        targs = sorted(ras, key=ras.get)

        with open(fname,'w') as f:
            f.write("""
<html>
<head>
<title>Complete list of ULTRACAM targets</title>
<link rel="stylesheet" type="text/css" href="ultracam_logs.css" />
</head>
<body>
<h1>RA-ordered list of ULTRACAM targets</h1>

<p> 
This table shows the identifier name, position and matching strings used
to attach coordinates to objects in the ULTRACAM database. Since ULTRACAM does
not talk to telescopes, this is pretty much all we have to go on. If you spot
problems please let me (trm) know about them. Matching is exact and case-sensitive.
The positions are ICRS. Where you see "~" in the matching strings, they actually 
count as blanks. If you are searching for a name to use for a previously observed
object while observing which will be recognised (good for you), use one of the 
match strings rather than the ID string.

<p>
You can search for runs on particular positions <a href="ulogs.php">here</a>. Clicking
on the IDs should take you to a list of runs. If it returns no runs, try reducing the
minimum run length.

<p>
<table>
<tr><th class="left">ID</th><th>RA</th><th>Dec</th><th class="left">Matching strings</th></tr>
""")

            for targ in targs:
                entry = self[targ]
                rad = entry['ra']
                decd = entry['dec']
                ra = dec2sexg(rad,False,2)
                dec = dec2sexg(decd,True,1)
                req = urllib.urlencode(
                    {'slimits' : 'manual', 'target' : '', 'delta' : 2., 'RA1' : rad-0.01, \
                     'RA2' : rad + 0.01, 'Dec1' : decd - 0.01, 'Dec2' : decd + 0.01, 'emin' : 10.}
                )
                f.write(
                    f'<tr><td class="left"><a href="ulogs.php?{req}">{targ}<td>{ra}</td><td>{dec}</td><td class="left">'
                )
                lnames = ' '.join([name.replace(' ','~') for name in entry['names']])
                f.write(lnames + '</td></tr>\n')
                f.write('</table>\n</body>\n</html>\n')

def dec2sexg(value, sign, dp):
    v = abs(value)
    v1 = int(v)
    v2 = int(60*(v-v1))
    v3 = 3600*(v-v1-v2/60)
    if sign:
        if v >= 0.:
            return f'+{v1:02d} {v2:02d} {v3:0{3+dp}.{dp}f}'
        else:
            return f'-{v1:02d} {v2:02d} {v3:0{3+dp}.{dp}f}'
    else:
        return f'{v1:02d} {v2:02d} {v3:0{3+dp}.{dp}f}'

def load_skip_fail():
    """Looks for a file called SKIP_TARGETS containing a list of target
    names to skip as far as target position lookup, and another called
    FAILED_TARGETS of targets that have failed target lokup. Returns
    a list of names to skip and a dictionary keyed by target name that
    returns a tuple containing run name, date and run number of the failed
    name.
    """
    if os.path.exists('SKIP_TARGETS'):
        with open('SKIP_TARGETS') as fp:
            skip_targets = fp.readlines()
            skip_targets = [name.strip() for name in skip_targets if not name.startswith('#')]
        print(f'Loaded {len(skip_targets)} target names to skip from SKIP_TARGETS')
    else:
        print('No file called SKIP_TARGETS')
        skip_targets = []

    failed_targets = {}
    if os.path.exists('FAILED_TARGETS'):
        with open('FAILED_TARGETS') as fp:
            for line in fp:
                if not line.startswith('#'):
                    target,rdir,ndir,run = line.split()
                    failed_targets[target] = (rdir,ndir,run)
                    skip_targets.append(target.replace('~',' '))

        print(f'Loaded {len(failed_targets)} target names from FAILED_TARGETS')
    else:
        print('No file called FAILED_TARGETS')

    return (skip_targets, failed_targets)


def str2radec(position):
    """ra,dec,system = str2radec(position) -- translates an astronomical
    coordinate string to double precision RA and Dec.

    'ra' is the RA in decimal hours; 'dec' is the declination in
    degrees; 'system' is one of 'ICRS', 'B1950', 'J2000'. Entering
    coordinates is an error-prone and often irritating chore.  This
    routine is designed to make it as safe as possible while
    supporting a couple of common formats.

    Here are example input formats, both good and bad:

     12 34 56.1 -23 12 12.1     -- OK. 'system' will default to ICRS
     234.5 34.2  B1950          -- OK. RA, Dec entered in decimal degrees, B1950 to override default ICRS
     11 02 12.1 -00 00 02       -- OK. - sign will be picked up
     12:32:02.4 -12:11 10.2     -- OK. Colon separators allowed.

     11 02 12.1 -23.2 J2000     -- NOK. Cannot mix HMS/DMS with decimals
     11 02 12.1 -23 01 12 J4000 -- NOK. Only 'ICRS', 'B1950', 'J2000' allowed at end.
     1.02323e2 -32.5            -- NOK. Scientific notation not supported.
     11 02 12.1 23 01 12        -- NOK. In HMS mode, the sign of the declination must be supplied
     25 01 61.2 +90 61 78       -- NOK. various items out of range.

    A ValueError is raised on failure

    """

    # Try out three types of match
    m = re.search(r'^\s*(\d{1,2})(?:\:|\s+)(\d{1,2})(?:\:|\s+)(\d{1,2}(?:\.\d*)?)\s+([\+-])(\d{1,2})(?:\:|\s+)(\d{1,2})(?:\:|\s+)(\d{1,2}(?:\.\d*)?)(?:\s+(\w+))?\s*$', position)
    if m:
        rah,ram,ras,decsign,decd,decm,decs,system = m.groups()
        rah  = int(rah)
        ram  = int(ram)
        ras  = float(ras)
        decd = int(decd)
        decm = int(decm)
        decs = float(decs)
        if (rah > 23 and ram > 0 and ras > 0.) or ram > 60 or ras > 60. \
                or (decd > 89 and decm > 0 and decs > 0.) or decm > 60 or decs > 60.:
            raise ValueError('one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

        if not system: system = 'ICRS'
        if system != 'ICRS' and system != 'J2000' and system != 'B1950':
            raise ValueError('astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + ' is not recognised.')

        ra  = rah + ram/60. + ras/3600.
        dec = decd + decm/60. + decs/3600.
        if decsign == '-': dec = -dec
        return (ra,dec,system)

    # No arcseconds of dec as sometimes is the case with coords from simbad
    m = re.search(r'^\s*(\d{1,2})(?:\:|\s+)(\d{1,2})(?:\:|\s+)(\d{1,2}(?:\.\d*)?)\s+([\+-])(\d{1,2})(?:\:|\s+)(\d{1,2}\.\d*)(?:\s+(\w+))?\s*$', position)
    if m:
        rah,ram,ras,decsign,decd,decm,system = m.groups()
        rah  = int(rah)
        ram  = int(ram)
        ras  = float(ras)
        decd = int(decd)
        decm = float(decm)
        if (rah > 23 and ram > 0 and ras > 0.) or ram > 60 or ras > 60. \
                or (decd > 89 and decm > 0.) or decm > 60.:
            raise ValueError('one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

        if not system: system = 'ICRS'
        if system != 'ICRS' and system != 'J2000' and system != 'B1950':
            raise ValueError('astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + ' is not recognised.')

        ra  = rah + ram/60. + ras/3600.
        dec = decd + decm/60.
        if decsign == '-': dec = -dec
        return (ra,dec,system)

    # Decimal entries
    m = re.search(r'^\s*(\d{1,3}(?:\.\d*)?)\s+([+-]?\d{1,2}(?:\.\d*)?)\s+(\w+)?\s*$', position)
    if m:
        print ('matched decimal entries')
        (rad,dec,system) = m.groups()
        ra   = float(rad)/15.
        dec  = float(dec)
        if ra >= 24. or dec < -90. or dec > 90.:
            raise ValueError('one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

        if not system: system = 'ICRS'
        if system != 'ICRS' and system != 'J2000' and system != 'B1950':
            raise ValueError('astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + ' is not recognised.')

        return (ra,dec,system)

    raise ValueError('could not interpret "' + position + '" as astronomical coordinates')

class TimeHMSCustom(TimeISO):
    """
    Just the HMS part of a Time as "<HH>:<MM>:<SS.sss...>".
    """

    name = "hms_custom"  # Unique format name
    subfmts = (("date_hms", "%H%M%S", "{hour:02d}:{min:02d}:{sec:02d}"),)

