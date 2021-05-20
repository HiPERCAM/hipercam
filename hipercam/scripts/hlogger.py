import sys
import traceback
import os
import time
import glob
import re
import math

import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
from astropy.coordinates import get_sun, get_moon, EarthLocation, SkyCoord, AltAz
import astropy.units as u

import hipercam as hcam
from hipercam.utils import format_hlogger_table

__all__ = [
    "hlogger",
]

############################################
#
# logger -- generates html logs for hipercam
#
############################################

###########
#
# Some constant stuff
#

# Header of the main file

INDEX_HEADER = """<html>
<head>

<title>HiPERCAM data logs</title>
</head>

<body>
<h1>HiPERCAM data logs</h1>

<p> These on-line logs summarise all HiPERCAM runs. They were
automatically generated from the FITS headers of the data files
and hand-written logs.

<p>For a searchable file summarising the same information and more,
see this <a href="hipercam-log.xlsx">spreadsheet</a>. This contains a
column 'Nlink' with hyperlinks to the on-line log of whichever night
contains a given run. In addition after the comment column the
spreadsheet has a number of angles (altitudes and azimuths etc, all in
degrees) and the lunar phase. 'Start' and 'end' refer to the mid-time
of the first and last frame of a run, while 'mid' is half-way between
them. These parameters could allow one to select dark runs for
instance.  """

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
var hcols = [4, 6, 7, 9, 15, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34];

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

<link rel="stylesheet" type="text/css" href="hiper.css" />
"""

NIGHT_HEADER2 = """<title>HiPERCAM log {0:s}</title>
</head>

<body>
<h1>HiPERCAM log {0:s}</h1>

<p>
The table below lists information on runs from the night starting on {0:s}.
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
<td><button id="hidden28" onclick="hide(28)"></button></td>
<td><button id="hidden29" onclick="hide(29)"></button></td>
<td><button id="hidden30" onclick="hide(30)"></button></td>
<td><button id="hidden31" onclick="hide(31)"></button></td>
<td><button id="hidden32" onclick="hide(32)"></button></td>
<td><button id="hidden33" onclick="hide(33)"></button></td>
<td><button id="hidden34" onclick="hide(34)"></button></td>
<td align="left"><button id="hidden35" onclick="hide(35)"></button></td>
</tr>

<tr>
<th class="left">Run<br>no.</th>
<th class="left">Target<br>name</th>
<th class="left">RA (J2000)<br>(telescope)</th>
<th class="left">Dec&nbsp;(J2000)<br>(telescope)</th>
<th class="cen">Date<br>(start)</th>
<th class="cen">Start</th>
<th class="cen">End</th>
<th class="cen">Time<br>flag</th>
<th class="right">Cadence<br>(sec)</th>
<th class="right">Duty<br>cycle, %</th>
<th class="right">Nframe</th>
<th class="right">Dwell<br>(sec)</th>
<th class="cen">Filters</th>
<th class="left">Type</th>
<th class="cen">Read<br>mode</th>
<th class="cen">Nskips</th>
<th class="cen">
xll,xlr,xul,xur,ys,nx,ny<br>
xsl,xsr,ys,nx,ny&nbsp;[DRIFT]<br>
Quad1</th>
<th class="cen">
xll,xlr,xul,xur,ys,nx,ny<br>
xsl,xsr,ys,nx,ny&nbsp;[DRIFT]<br>
Quad2</th>
<th class="cen">XxY<br>bin</th>
<th class="cen">Clr</th>
<th class="cen">Dum</th>
<th class="cen">LED</th>
<th class="cen">Over-<br>scan</th>
<th class="cen">Pre-<br>scan</th>
<th class="cen">Nod</th>
<th class="cen">Readout<br>speed</th>
<th class="cen">Fast<br>clocks</th>
<th class="cen">Tbytes</th>
<th class="cen">FPslide</th>
<th class="cen">Instr.<br>PA</th>
<th class="cen">CCD<br>temps</th>
<th class="cen">Observers</th>
<th class="left">PI</th>
<th class="left">PID</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
"""

NIGHT_FOOTER = """

<p> 'Instr. PA' is the instrumental PA. On the GTC this is measured
East of North in degrees, but has an offset which may vary somewhat
from night-to-night. 'Clr' indicates whether clears were enabled;
'Dum' whether the dummy output was being used to reduce pickup noise;
'Over-scan' and 'Pre-scan' whether the over- and pre-scans were on or
off; 'Nod' whether each exposure was paused to allow the telescope to
be moved. 'Time flag' is set to NOK (not OK) if the number of
satellites == -1 and the GPS was not synced or if the time that came
back could not be understood. 'Read mode' is the readout mode which
can be one of 4 options, namely 'FULL' for a full frame, 1-WIN or
2-WIN for standard windows mode, and 'DRIFT' for drift-mode. The
'xsll,xslr,..' column gives the 7 parameters that appear at the bottom
of the top-right window of hdriver (or 5 parameters in DRIFT mode).
In hdriver these are called xsll, xslr, xsul, xsur, ys, nx and ny (or
xsl, xsr, ys, nx, ny in drift mode), and there are up to two sets of
them. In fullframe mode, these parameters do not need to be
specified. 'Nskips' refers to how often each CCD is read out. The
numbers correspond to nu, ng, nr, ni, nz listed in hdriver. 'Duty
cycle' shows the fractional time collecting photons; the 'cadence' is
the (minimum) time between exposures, minimised over the different
arms. 'Tbytes' says whether there are the correct number of timing
bytes (36).  </p>

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

# HiPERCAM CSS to define the look of the log web pages
# This is written to 'hiper.css' and referred to in the
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

TRANSLATE_MODE = {
    "FullFrame": "FULL",
    "OneWindow": "1-WIN",
    "TwoWindows": "2-WIN",
    "DriftWindow": "DRIFT",
}


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

observatory = EarthLocation.of_site('Roque de los Muchachos')

def hlogger(args=None):
    """``hlogger server (dirnam)``

    Generates html logs for hipercam runs.

    hlogger expects to work on directories containing the runs for
    each night. It should be run from the hipercam/raw_data directory.
    It extracts information from each run file. It runs without any
    parameters.  hlogger also writes out a spreadsheet with useful
    info line-by-line for each run.

    If run at Warwick, it writes data to the web pages. Otherwise it
    writes them to a sub-directory of raw_data called "logs".

    Bit of a specialist routine this one; if you have access to the
    on-line logs at Warwick, it should be unnecessary. It requires the
    installation of the python module xlsxwriter in order to write an
    excel spreadsheet; this is not a required module for the whole
    pipeline and may not exist. The script will fail early on if it is
    not installed.

    """

    barr = []

    if os.path.basename(os.getcwd()) != "raw_data":
        print("** hlogger must be run in a directory called 'raw_data'")
        print("hlogger aborted", file=sys.stderr)
        return

    # location to write files
    if os.path.exists('/storage/astro2/www/phsaap/hipercam/logs'):
        root = '/storage/astro2/www/phsaap/hipercam/logs'
    else:
        root = 'logs'

    # Try an import now to get an early warning of problems
    import xlsxwriter

    # next are regular expressions to match run directories, nights, and
    # run files
    rre = re.compile("^\d\d\d\d(-\d\d|[AB])$")
    nre = re.compile("^\d\d\d\d-\d\d-\d\d$")
    fre = re.compile("^run\d\d\d\d\.fits$")

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
        print("hlogger aborted", file=sys.stderr)
        return

    # Get started. First make sure all files are created with the
    # right permissions
    os.umask(0o022)

    # Ensure the root directory exists.
    os.makedirs(root, exist_ok=True)

    print(f'Will write to directory = "{root}".')

    # Index file
    index = os.path.join(root, 'index.html')
    with open(index, "w") as ihtml:
        # start the top level index html file
        ihtml.write(INDEX_HEADER)

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
            ihtml.write("\n<p>\n<table>\n")

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
                print("hlogger aborted", file=sys.stderr)
                return

            for nn, nname in enumerate(nnames):

                print(f"  night {nname}")

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
                    f'<tr><td><a href="{night}.html">{night}</a><td></tr>\n'
                )

                # Create the html file for the night
                date = f"{night}, {telescope}"
                fname = os.path.join(root, f"{night}.html")

                with open(fname, "w") as nhtml:

                    # write header of night file
                    nhtml.write(NIGHT_HEADER1)
                    nhtml.write(NIGHT_HEADER2.format(date))
                    nhtml.write(links)
                    nhtml.write(TABLE_TOP)

                    # read and store the hand written log
                    handlog = os.path.join(night, f"{night}.dat")
                    with open(handlog) as fin:
                        hlog = {}
                        for line in fin:
                            if line.startswith("run"):
                                arr = line.split()
                                hlog[arr[0]] = " ".join(arr[1:])

                    # load all the run names
                    runs = [run[:-5] for run in os.listdir(night) if fre.match(run)]
                    runs.sort()

                    # now wind through the runs getting basic info and
                    # writing a row of info to the html file for the night
                    # in question.
                    for nrun, run in enumerate(runs):

                        if nrun % 20 == 0:
                            # write table header
                            nhtml.write(TABLE_HEADER)

                        # open the run file as an Rtime
                        runname = os.path.join(night, run)
                        try:
                            rtime = hcam.hcam.Rtime(runname)
                        except:
                            exc_type, exc_value, exc_traceback = sys.exc_info()
                            traceback.print_tb(
                                exc_traceback, limit=1, file=sys.stdout
                            )
                            traceback.print_exc(file=sys.stdout)
                            print("Problem on run = ", runname)

                            # dummy info line just to allow us to proceed
                            nhtml.write("<tr>\n")
                            # run number
                            nhtml.write(f'<td class="lalert">{run[3:]}</td>')
                            nhtml.write("</tr>\n")
                            brow = [run[3:]] + 23*[None]
                            continue

                        hd = rtime.header

                        # start the row
                        nhtml.write("<tr>\n")

                        # run number
                        runno = run[3:]
                        nhtml.write('<td class="left">{:s}</td>'.format(runno))

                        # this is the start of a list to be appended
                        # to the array sent to a pandas.DataFrame at
                        # the end of the script in order to write out
                        # a spreadsheet
                        #link = f'=HYPERLINK("http://deneb.astro.warwick.ac.uk/phsaap/hipercam/logs/{night}.html", "{runno}")'
                        #brow = [link,]
                        brow = [runno,]

                        # object name
                        nhtml.write(f'<td class="left">{hd["OBJECT"]}</td>')
                        brow.append(hd["OBJECT"])

                        # RA, Dec
                        ra, dec = correct_ra_dec(hd["RA"], hd["Dec"])
                        nhtml.write(f'<td class="left">{ra}</td><td class="left">{dec}</td>')
                        brow += [ra, dec, rname, night]

                        # timing info
                        ntotal = rtime.ntotal()
                        texps, toffs, nskips, tdead = rtime.tinfo()
                        # total = total time on target
                        # duty = worst duty cycle, percent
                        # tsamp = shortest sample time
                        ttotal = 0.0
                        duty = 100
                        tsamp = 99000.0
                        for texp, nskip in zip(texps, nskips):
                            ttotal = max(ttotal, (texp + tdead) * (ntotal // nskip))
                            duty = min(duty, 100.0 * texp / (texp + tdead))
                            tsamp = min(tsamp, texp + tdead)

                        # First & last timestamp
                        try:

                            # get mid-exposure time of first and last frame.
                            tstamp_start, tinfo, tflag1 = rtime(1)
                            tstamp_end, tinfo, tflag2 = rtime(ntotal)

                            # compute mid-time
                            tstamp_mid = tstamp_start + (tstamp_end-tstamp_start)/2

                            # Extract understandable times
                            tstart = tstamp_start.isot
                            tend = tstamp_end.isot

                            datestart = tstart[:tstart.find("T")]
                            utcstart = tstart[tstart.find("T")+1:tstart.rfind(".")]
                            utcend = tend[tend.find("T")+1:tend.rfind(".")]
                            tflag = "OK" if tflag1 and tflag2 else "NOK"
                            nhtml.write(
                                f'<td class="cen">{datestart}</td> <td class="cen">{utcstart}</td>' +
                                f'<td class="cen">{utcend}</td> <td class="cen">{tflag}</td>'
                            )
                            brow += [datestart, utcstart, utcend, tflag]

                        except:
                            exc_type, exc_value, exc_traceback = sys.exc_info()
                            traceback.print_tb(
                                exc_traceback, limit=1, file=sys.stdout
                            )
                            traceback.print_exc(file=sys.stdout)
                            print("Run =", runname)
                            nhtml.write(
                                '<td class="cen">----</td><td class="cen">----</td><td class="cen">----</td><td>NOK</td>'
                            )
                            brow += 4*[None]

                            # set start to a bad value for later
                            tstamp_start = None

                        # sample time
                        nhtml.write(f'<td class="right">{tsamp:.3f}</td>')
                        brow.append(round(tsamp,4))

                        # duty cycle
                        nhtml.write(f'<td class="right">{duty:.1f}</td>')
                        brow.append(round(duty,1))

                        # number of frames
                        nhtml.write(f'<td class="right">{ntotal:d}</td>')
                        brow.append(ntotal)

                        # total exposure time
                        ttime = int(round(ttotal))
                        nhtml.write(f'<td class="right">{ttime:d}</td>')
                        brow.append(ttime)

                        # filters used
                        filters = hd.get("filters", "----")
                        nhtml.write(f'<td class="cen">{filters}</td>')
                        brow.append(filters)

                        # run type
                        itype = hd.get("IMAGETYP", "----")
                        nhtml.write(f'<td class="left">{itype}</td>')
                        brow.append(itype)

                        # readout mode
                        nhtml.write(f'<td class="cen">{TRANSLATE_MODE[rtime.mode]}</td>')
                        brow.append(TRANSLATE_MODE[rtime.mode])

                        # cycle nums
                        skips = ",".join([str(nskip) for nskip in nskips])
                        nhtml.write(f'<td class="cen">{skips}</td>')
                        brow.append(skips)

                        # window formats
                        win1 = rtime.wforms[0]
                        nhtml.write(f'<td class="cen">{win1}</td>')
                        win1 = win1.replace('&nbsp;',' ')
                        brow.append(win1)

                        win2 = rtime.wforms[1] if len(rtime.wforms) > 1 else ""
                        nhtml.write(f'<td class="cen">{win2}</td>')
                        win2 = win2.replace('&nbsp;',' ')
                        brow.append(win2)

                        # binning
                        binning = f'{rtime.xbin:d}x{rtime.ybin:d}'
                        nhtml.write(f'<td class="cen">{binning}</td>')
                        brow.append(binning)

                        # clear
                        clear = "On" if rtime.clear else "Off"
                        nhtml.write(f'<td class="cen">{clear}</td>')
                        brow.append(clear)

                        # dummy output in use
                        dummy = "On" if rtime.dummy else "Off"
                        nhtml.write(f'<td class="cen">{dummy}</td>')
                        brow.append(dummy)

                        # LED on
                        led = hd.get("ESO DET EXPLED", "----")
                        led = "On" if led == 1 else "Off"

                        nhtml.write(f'<td class="cen">{led}</td>')
                        brow.append(led)

                        # over-scan
                        oscan = "On" if rtime.oscan else "Off"
                        nhtml.write(f'<td class="cen">{oscan}</td>')
                        brow.append(oscan)

                        # pre-scan
                        pscan = "On" if rtime.pscan else "Off"
                        nhtml.write(f'<td class="cen">{pscan}</td>')
                        brow.append(pscan)

                        # Nodding
                        nod = hd.get("ESO DET SEQ1 TRIGGER", 0)
                        nod = "On" if nod == 1 else "Off"
                        nhtml.write(f'<td class="cen">{nod}</td>')
                        brow.append(nod)

                        # CCD speed
                        speed = hd.get("ESO DET SPEED", "----")
                        if speed == 0:
                            speed = "Slow"
                        elif speed == 1:
                            speed = "Fast"
                        nhtml.write(f'<td class="cen">{speed}</td>')
                        brow.append(speed)

                        # Fast clocks
                        fclock = hd.get("ESO DET FASTCLK", "----")
                        if fclock == 0:
                            fclock = "No"
                        elif fclock == 1:
                            fclock = "Yes"
                        nhtml.write(f'<td class="cen">{fclock}</td>')
                        brow.append(fclock)

                        # Tbytes problem
                        tbytes = "OK" if rtime.ntbytes == 36 else "NOK"
                        nhtml.write(f'<td class="cen">{tbytes}</td>')
                        brow.append(tbytes)

                        # Focal plane slide
                        fpslide = hd.get("FPslide", "----")
                        nhtml.write(f'<td class="cen">{fpslide}</td>')
                        brow.append(fpslide)

                        # instr PA [GTC]
                        instpa = hd.get("INSTRPA", "----")
                        nhtml.write(f'<td class="cen">{instpa}</td>')
                        brow.append(instpa)

                        # CCD temps
                        t1,t2,t3,t4,t5 = hd.get("CCD1TEMP", 0.0), hd.get("CCD2TEMP", 0.0), \
                            hd.get("CCD3TEMP", 0.0), hd.get("CCD4TEMP", 0.0), hd.get("CCD5TEMP", 0.0),

                        ccdtemps = f'{t1:.1f},{t2:.1f},{t3:.1f},{t4:.1f},{t5:.1f}'
                        nhtml.write(f'<td class="cen">{ccdtemps}</td>')
                        brow.append(ccdtemps)

                        # Observers
                        observers = hd.get("OBSERVER", "----")
                        nhtml.write(f'<td class="cen">{observers}</td>')
                        brow.append(observers)

                        # PI
                        pi = hd.get("PI", "----")
                        nhtml.write(f'<td class="left">{pi}</td>')
                        brow.append(pi)

                        # Program ID
                        pid = hd.get("PROGRM", "----")
                        nhtml.write(f'<td class="left">{pid}</td>')
                        brow.append(pid)

                        # run number again. Add null to the spreadsheet
                        # to allow a formula
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

                        comments = f'{pcomm}{hlog[run]}'
                        nhtml.write(f'<td class="left">{comments}</td>')
                        brow.append(comments)

                        # Finally tack on extra positional stuff for the spreadsheet only
                        if tstamp_start:

                            if ra != '' and dec != '':
                                # Target stuff. Calculate the Alt & Az at start, middle, end
                                frames = AltAz(obstime=[tstamp_start,tstamp_mid,tstamp_end], location=observatory)
                                points = SkyCoord(f'{ra} {dec}',unit=(u.hourangle, u.deg)).transform_to(frames)
                                alts = [round(alt,1) for alt in points.alt.degree]
                                azs = [round(az,1) for az in points.az.degree]
                                brow += alts + azs

                                # Now calculate the angular distance from the Sun and Moon at the mid-time
                                sun_mid = get_sun(tstamp_mid).transform_to(frames[1])
                                moon_mid = get_moon(tstamp_mid).transform_to(frames[1])
                                point_mid = points[1]
                                sun_dist = point_mid.separation(sun_mid).degree
                                moon_dist = point_mid.separation(moon_mid).degree
                                brow += [round(sun_dist,1),round(moon_dist,1)]
                            else:
                                brow += 8*[None]


                            # Now some data on the altitude of the Sun & Moon
                            frame = AltAz(obstime=tstamp_start, location=observatory)
                            sun_start = get_sun(tstamp_start).transform_to(frame)
                            moon_start = get_moon(tstamp_start).transform_to(frame)

                            # end
                            frame = AltAz(obstime=tstamp_end, location=observatory)
                            sun_end = get_sun(tstamp_end).transform_to(frame)
                            moon_end = get_moon(tstamp_end).transform_to(frame)

                            # Lunar phase at the mid-point only. Need re-doing.
                            #moon_phase = moon_phase_angle(tstamp_mid).value / np.pi
                            moon_phase = 0
                            
                            # write out info
                            brow += [
                                round(sun_start.alt.degree,1), round(sun_end.alt.degree,1),
                                round(moon_start.alt.degree,1), round(moon_end.alt.degree,1),
                                round(moon_phase,2),
                            ]

                        else:
                            brow += 13*[None]

                        # at last: end the row
                        nhtml.write("\n</tr>\n")
                        barr.append(brow)

                    # finish off the night file
                    nhtml.write("</table>\n{:s}".format(links))
                    nhtml.write(NIGHT_FOOTER)

            # finish off the run table
            ihtml.write("\n</table>\n\n")

        # finish the main index
        ihtml.write(INDEX_FOOTER)

    # write out the css file
    css = os.path.join(root, "hiper.css")
    with open(css, "w") as fout:
        fout.write(CSS)

    # Create and write out spreadsheet
    colnames = (
        'Run no.', 'Target name', 'RA (J2000)', 'Dec (J2000)', 'Obs run',
        'Night', 'Date\n(start)', 'UTC\nStart', 'UTC\nEnd', 'Tflag',
        'Cadence\n(sec)', 'Duty\ncycle(%)', 'Nframe', 'Dwell\n(sec)', 'Filters', 'Run\ntype',
        'Readout\nmode', 'Nskips', 'Win1', 'Win2' , 'XxY\nbin', 'Clr', 'Dum',
        'LED', 'Over-\nscan', 'Pre-\nscan', 'Nod', 'Readout\nspeed', 'Fast\nclocks', 'Tbytes',
        'FPslide', 'Instr\nPA', 'CCD temperatures', 'Observers', 'PI', 'PID', 'Nlink', 'Comment',
        'Alt\nstart', 'Alt\nmiddle', 'Alt\nend', 'Az\nstart', 'Az\nmiddle', 'Az\nend',
        'Sun sep\nmid', 'Moon sep\nmid',
        'Sun alt\nstart', 'Sun alt\nend', 'Moon alt\nstart', 'Moon alt\nend', 'Moon\nphase',
    )
    ptable = pd.DataFrame(barr, columns=colnames)

    spreadsheet = os.path.join(root, "hipercam-log.xlsx")

    format_hlogger_table(spreadsheet, ptable)

    print(f'\nAll done. Look in {root} for the outputs.')

