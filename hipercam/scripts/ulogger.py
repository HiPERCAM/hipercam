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
<h1>HiPERCAM data logs</h1>

<p> These on-line logs summarise all {instrument} runs. They were
automatically generated from the data files and hand-written logs.

<p>For a searchable file summarising the same information and more,
see this <a href="{linstrument}-log.xlsx">spreadsheet</a>. This contains a
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
        print("** ulogger: cannot finf either ultracam or ultraspec in path")
        print("ulogger aborted", file=sys.stderr)
        return
    linstrument = instrument.lower()

    # location to write files
    if os.path.exists(f'/storage/astro2/www/phsaap/{linstrument}/logs'):
        root = f'/storage/astro2/www/phsaap/{linstrument}/logs'
    else:
        root = 'logs'

    # Try an import now to get an early warning of problems
    import xlsxwriter

    # next are regular expressions to match run directories, nights, and
    # run files
    rre = re.compile("^\d\d\d\d-\d\d$")
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
                print("ulogger aborted", file=sys.stderr)
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
                    runs = [run[:-4] for run in os.listdir(night) if fre.match(run)]
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
                            rtime = hcam.ucam.Rtime(runname)
                        except hcam.ucam.PowerOnOffError:
                            # not interesting
                            continue
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
                        target = hd["TARGET"]
                        nhtml.write(f'<td class="left">{target}</td>')
                        brow.append(target)

                        # RA, Dec lookup
                        if target in skip_targets:
                            # don't even try
                            ra, dec = '', ''
                        elif target in targets.lnames:
                            dct = targets.lnames[target]
                            ra, dec = dct['ra'], dct['dec']
                        elif target in auto_targets.lnames:
                            dct = auto_targets.lnames[target]
                            ra, dec = dct['ra'], dct['dec']
                        else:
                            # add simbad lookup here
                            ra, dec = '', ''

                        nhtml.write(f'<td class="left">{ra}</td><td class="left">{dec}</td>')
                        brow += [ra, dec, rname, night]

                        # timing info
                        utime1, info1 = rtime(1)[:2]
                        utime2, info2 = rtime(0)[:2]
                        mjd1 = utime1.mjd
                        mjd2 = utime1.mjd
                        ntotal = info2['fnum']
                        cadence = 86400.*(mjd2-mjd1)/max(1,ntotal-1)
                        ttime = 1440.*(mjd2-mjd1)

                        tstart = Time(mjd1, format='mjd').isot
                        tend = Time(mjd2, format='mjd').isot

                        datestart = tstart[:tstart.find("T")]
                        utcstart = tstart[tstart.find("T")+1:tstart.rfind(".")]
                        utcend = tend[tend.find("T")+1:tend.rfind(".")]

                        nhtml.write(
                            f'<td class="cen">{datestart}</td> <td class="cen">{utcstart}</td>' +
                            f'<td class="cen">{utcend}</td>'
                        )
                        brow += [datestart, utcstart, utcend]

                        # cadence time
                        nhtml.write(f'<td class="right">{cadence:.3f}</td>')
                        brow.append(round(cadence,4))

                        # number of frames
                        nhtml.write(f'<td class="right">{ntotal:d}</td>')
                        brow.append(ntotal)

                        # total exposure time
                        ttime = int(round(ttime))
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
                        nhtml.write(f'<td class="cen">{rtime.mode}</td>')
                        brow.append(rtime.mode)

                        # nblue
                        nblue = hd['NBLUE']
                        nhtml.write(f'<td class="cen">{nblue}</td>')
                        brow.append(nblue)

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
                        binning = f'{rtime.xbin}x{rtime.ybin}'
                        nhtml.write(f'<td class="cen">{binning}</td>')
                        brow.append(binning)

                        # clear
                        clear = "On" if hd['CLEAR'] else "Off"
                        nhtml.write(f'<td class="cen">{clear}</td>')
                        brow.append(clear)

                        # CCD speed
                        speed = hd['GAINSPED']
                        nhtml.write(f'<td class="cen">{speed}</td>')
                        brow.append(speed)

                        # Focal plane slide
                        fpslide = hd.get('SLIDEPOS','UNKNOWN')
                        nhtml.write(f'<td class="cen">{fpslide}</td>')
                        brow.append(fpslide)

                        if instrument == 'ULTRASPEC':
                            # instr PA [GTC]
                            instpa = hd.get("INSTRPA", "----")
                            nhtml.write(f'<td class="cen">{instpa}</td>')
                            brow.append(instpa)

                        # Observers
                        observers = hd.get("OBSERVERS", "----")
                        nhtml.write(f'<td class="cen">{observers}</td>')
                        brow.append(observers)

                        # PI
                        pi = hd.get("PI", "----")
                        nhtml.write(f'<td class="left">{pi}</td>')
                        brow.append(pi)

                        # Program ID
                        pid = hd.get("ID", "----")
                        nhtml.write(f'<td class="left">{pid}</td>')
                        brow.append(pid)

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

                        comments = f'{pcomm}{hlog[run]}'
                        nhtml.write(f'<td class="left">{comments}</td>')
                        brow.append(comments)

                        # Finally tack on extra positional stuff for the spreadsheet only
                        if tstart:

                            if ra != '' and dec != '':
                                # Target stuff. Calculate the Alt & Az at start, middle, end
                                frames = AltAz(obstime=[tstart,tmid,tend], location=observatory)
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
                            frame = AltAz(obstime=tstart, location=observatory)
                            sun_start = get_sun(tstart).transform_to(frame)
                            moon_start = get_moon(tstart).transform_to(frame)

                            # end
                            frame = AltAz(obstime=tend, location=observatory)
                            sun_end = get_sun(tend).transform_to(frame)
                            moon_end = get_moon(tend).transform_to(frame)

                            # Lunar phase at the mid-point only.
                            moon_phase = moon_phase_angle(tmid).value / np.pi

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
    css = os.path.join(root, "ultra.css")
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

    spreadsheet = os.path.join(root, f"{linstrument}-log.xlsx")

    format_hlogger_table(spreadsheet, ptable)

    print(f'\nAll done. Look in {root} for the outputs.')


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
                            ra,dec,system = subs.str2radec(ra + ' ' + dec)
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

    def write(self, fname):
        """
        Write targets out to disk file fname.
        """

        # write in RA order
        ras   = dict([(targ,entry['ra']) for targ, entry in self.items()])
        targs = sorted(ras, key=ras.get)

        with open(fname,'w') as f:
            f.write("""
#
# File of targets written by ulogger.Targets.write
#

""")

            for targ in targs:
                entry = self[targ]
                pos = subs.d2hms(entry['ra'],dp=2) + '   ' + subs.d2hms(entry['dec'],dp=1,sign=True)
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
                ra = subs.d2hms(rad,sep=' ')
                dec = subs.d2hms(decd,sep=' ',sign=True)
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
            skip_targets = [name.strip() for name in sskip if not name.startswith('#')]
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

