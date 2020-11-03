import sys
import traceback
import os
import time
import glob
import re

import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

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

<p>
These on-line logs summarise all HiPERCAM runs. They were automatically generated from
the data files and hand-written logs. For a searchable file summarising the same
information, see this <a href="hipercam-log.xlsx">spreadsheet</a>.
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
var hcols = [4, 6, 7, 9, 15, 19, 20, 21, 22, 24, 25, 26, 27, 28];

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
You can hide any column using the buttons in the top line and "shrink" will reduce clutter
to focus on those of most interest.
"""

TABLE_TOP = """<p>
<button style="width: 100px; height: 25px;" onclick="shrink();">Shrink table</button>
<button style="width: 125px; height: 25px;" onclick="restore();">Restore columns</button>


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
<td align="left"><button id="hidden32" onclick="hide(32)"></button></td>
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
<th class="cen">Ov-/Pre-<br>scan</th>
<th class="cen">Readout<br>speed</th>
<th class="cen">Fast<br>clocks</th>
<th class="cen">Tbytes</th>
<th class="cen">FPslide</th>
<th class="cen">Instr.<br>PA</th>
<th class="cen">CCD<br>temps</th>
<th class="left">PI</th>
<th class="left">PID</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
"""

NIGHT_FOOTER = """

<p> 'Instr. PA' is the instrumental PA. On the GTC this is measured East of North in degrees, but has
an offset which may vary somewhat from night-to-night. 'Clr' indicates whether clears were enabled; 'Dum' whether the dummy
output was being used to reduce pickup noise; 'Ov-Pre-scan' whether the over-
and pre-scans were on or off. 'Time flag' is set to NOK (not OK) if the number
of satellites == -1 and the GPS was not synced or if the time that came back
could not be understood. 'Read mode' is the readout mode which can be one of 4
options, namely 'FULL' for a full frame, 1-WIN or 2-WIN for standard windows
mode, and 'DRIFT' for drift-mode. The 'xsll,xslr,..' column gives the 7
parameters that appear at the bottom of the top-right window of hdriver (or 5
parameters in DRIFT mode).  In hdriver these are called xsll, xslr, xsul,
xsur, ys, nx and ny (or xsl, xsr, ys, nx, ny in drift mode), and there are up
to two sets of them. In fullframe mode, these parameters do not need to be
specified. 'Nskips' refers to how often each CCD is read out. The numbers
correspond to nu, ng, nr, ni, nz listed in hdriver. 'Duty cycle' shows the
fractional time collecting photons; the 'cadence' is the (minimum) time
between exposures, minimised over the different arms. 'Tbytes' says whether
there are the correct number of timing bytes (36).  </p>

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

    # next are regular expressions to match run directories, nights, and
    # run files
    rre = re.compile("^\d\d\d\d-\d\d$")
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

    # Get started. First make sure all files are created with the right permissions
    os.umask(0o022)

    # Ensure the root directory exists.
    os.makedirs(root, exist_ok=True)

    print(f'Will write to {root}.')

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

                links = '\n<p><a href="../">Run index</a>'
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
                        nhtml.write('<td class="left">{:s}</td>'.format(run[3:]))
                        link = f'=HYPERLINK("http://deneb.astro.warwick.ac.uk/phsaap/hipercam/logs/{night}/{night}.html", {run[3:]})'
                        brow = [link,]

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
                            tstamp, tinfo, tflag1 = rtime(1)
                            tstart = tstamp.isot
                            tstamp, tinfo, tflag2 = rtime(ntotal)
                            tend = tstamp.isot
                            nhtml.write(
                                f'<td class="cen">{tstart[:tstart.find("T")]}</td>' +
                                f'<td class="cen">{tstart[tstart.find("T")+1:tstart.rfind(".")]}</td>' +
                                f'<td class="cen">{tend[tend.find("T") + 1 : tend.rfind(".")]}</td>' +
                                f'<td class="cen">{"OK" if tflag1 and tflag2 else "NOK"}</td>'
                            )
                            brow += [
                                tstart[:tstart.find("T")],
                                tstart[tstart.find("T")+1:tstart.rfind(".")],
                                tend[tend.find("T")+1:tend.rfind(".")],
                                "OK" if tflag1 and tflag2 else "NOK",
                            ]

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

                        # sample time
                        nhtml.write(f'<td class="right">{tsamp:.3f}</td>')
                        brow.append(f'{tsamp:.3f}')

                        # duty cycle
                        nhtml.write(f'<td class="right">{duty:.1f}</td>')

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
                        itype = hd.get("IMAGETYP", "---")
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
                        led = hd.get("ESO DET EXPLED", "--")
                        if led == 0:
                            led = "Off"
                        elif led == 1:
                            led = "On"

                        nhtml.write(f'<td class="cen">{led}</td>')
                        brow.append(led)

                        # overscan/prescan
                        opscan = f'{"On" if rtime.oscan else "Off"},{"On" if rtime.pscan else "Off"}'
                        nhtml.write(f'<td class="cen">{opscan}</td>')
                        brow.append(opscan)

                        # CCD speed
                        speed = hd.get("ESO DET SPEED", "--")
                        if speed == 0:
                            speed = "Slow"
                        elif speed == 1:
                            speed = "Fast"
                        nhtml.write(f'<td class="cen">{speed}</td>')
                        brow.append(speed)

                        # Fast clocks
                        fclock = hd.get("ESO DET FASTCLK", "--")
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
                        instpa = hd.get("INSTRPA", "UNDEF")
                        nhtml.write(f'<td class="cen">{instpa}</td>')
                        brow.append(instpa)

                        # CCD temps
                        t1,t2,t3,t4,t5 = hd.get("CCD1TEMP", 0.0), hd.get("CCD2TEMP", 0.0), \
                            hd.get("CCD3TEMP", 0.0), hd.get("CCD4TEMP", 0.0), hd.get("CCD5TEMP", 0.0),

                        ccdtemps = f'{t1:.1f},{t2:.1f},{t3:.1f},{t4:.1f},{t5:.1f}'
                        nhtml.write(f'<td class="cen">{ccdtemps}</td>')
                        brow.append(ccdtemps)

                        # PI
                        pi = hd.get("PI", "---")
                        nhtml.write(f'<td class="left">{pi}</td>')
                        brow.append(pi)

                        # Program ID
                        pid = hd.get("PROGRM", "---")
                        nhtml.write(f'<td class="left">{pid}</td>')
                        brow.append(pid)

                        # run number again
                        nhtml.write(f'<td class="left">{run[3:]}</td>')

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
        'Cadence\n(sec)', 'Nframe', 'Dwell\n(sec)', 'Filters', 'Run\ntype',
        'Readout\nmode', 'Nskips', 'Win1', 'Win2' , 'XxY\nbin', 'Clr', 'Dum',
        'LED', 'Ov-/Pre-\nscan', 'Readout\nspeed', 'Fast\nclocks', 'Tbytes',
        'FPslide', 'Instr\nPA', 'CCD temperatures', 'PI', 'PID', 'Comment'
    )

    spreadsheet = os.path.join(root, "hipercam-log.xlsx")
    writer = pd.ExcelWriter(spreadsheet, engine='xlsxwriter')
    ptable = pd.DataFrame(barr, columns=colnames)
    ptable.to_excel(writer, sheet_name='Hipercam logs', index=False)
    workbook = writer.book
    worksheet = writer.sheets["Hipercam logs"]

    # Column widths. Should probably automate.
    worksheet.set_column('B:B', 30)
    worksheet.set_column('C:D', 11)
    worksheet.set_column('F:G', 10)
    worksheet.set_column('H:I', 9)
    worksheet.set_column('J:J', 7)
    worksheet.set_column('L:M', 8)
    worksheet.set_column('N:O', 12)
    worksheet.set_column('P:P', 8)
    worksheet.set_column('Q:Q', 10)
    worksheet.set_column('R:S', 16)
    worksheet.set_column('T:W', 4)
    worksheet.set_column('X:Y', 8)
    worksheet.set_column('Z:AC', 7)
    worksheet.set_column('AD:AD', 26)
    worksheet.set_column('AE:AF', 14)
    worksheet.set_column('AG:AG', 60)

    worksheet.set_row(0,28)
    worksheet.freeze_panes(1, 0)
    worksheet.set_zoom(200)
    writer.save()

    print(f'\nAll done. Look in {root} for the outputs.')

