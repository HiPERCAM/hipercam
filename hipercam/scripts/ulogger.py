import sys
import traceback
import os
import shutil
import re
import warnings
import argparse
import sqlite3

import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta, TimeISO
from astropy.coordinates import get_sun, get_moon, EarthLocation, SkyCoord, AltAz
import astropy.units as u

import hipercam as hcam
from hipercam.utils import format_ulogger_table, target_lookup, dec2sexg, str2radec, LOG_CSS, LOG_MONTHS

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

def ulogger(args=None):
    description = \
    """``ulogger``

    Generates html logs for ULTRACAM and ULTRASPEC runs, a searchable
    spreadsheet and an sqlite3 database for programmatic SQL
    enquiries.

    ulogger expects to work in a directory containing the runs for
    each night. Specifically, it should be run from the
    ultra(cam|spec)/raw_data directory.  It extracts information from
    each run file. It usually runs without any parameters, although
    there are some optional unix-style command-line switches. It tries
    to read target data from four files: TARGETS, AUTO_TARGETS,
    FAILED_TARGETS and SKIP_TARGETS.

    If run at Warwick, it writes data to the web pages. Otherwise it
    writes them to a sub-directory of raw_data called "logs".  Bit of
    a specialist routine this one; if you have access to the on-line
    logs at Warwick, it should be unnecessary. It requires the
    installation of the python module xlsxwriter in order to write an
    excel spreadsheet; this is not a required module for the whole
    pipeline and may not exist. The script will fail early on if it is
    not installed.

    It also writes information to a sub-directory "meta" of each night
    directory, and can pick up information stored there by the related
    scripts |redplt| and |hmeta|. There is a slightly circular
    relationship between these scripts since |redplt| can only be run
    once ulogger has created a timing file in "meta". Thus two runs of
    |ulogger| are needed to create complete logs with images. A full
    sequence would be |ulogger|, followed by |redplt| and |hmeta|, and
    only when these two have completed, then run |ulogger| a second
    time. |hmeta| can be run idependently however, but must be run
    before the final |ulogger| run for a complete spreadsheet and SQL
    database to be created.

    """
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        dest="full",
        action="store_true",
        help="carry out full re-computation of times and positional data",
    )
    parser.add_argument(
        "-p",
        dest="positions",
        action="store_true",
        help="re-compute positional data but not times",
    )
    parser.add_argument(
        "-n",
        dest="night", default=None,
        help="use this with a YYYY-MM-DD date to update times and positions for a specific night (lots of diagnostic ouput; no html created)",
    )

    # Get command line options
    args = parser.parse_args()
    do_full = args.full
    do_positions = args.positions
    do_night = args.night

    if (do_full or do_positions) and do_night is not None:
        print('-n switch is not compatible with either -f or -p; please check help')
        return

    # start with defining a few regular expressions to match run directories, nights,
    # and run files
    rre = re.compile("^(Others|\d\d\d\d-(\d\d|P\d\d\d))$")
    nre = re.compile("^\d\d\d\d-\d\d-\d\d$")
    fre = re.compile("^run\d\d\d\.xml$")


    if do_night:

        # identify observatory
        with open(os.path.join(do_night, "telescope")) as tel:
            telescope = tel.read().strip()
        if telescope == 'WHT':
            observatory = EarthLocation.of_site('Roque de los Muchachos')
        elif telescope == 'VLT':
            observatory = EarthLocation.of_site('Cerro Paranal')
        elif telescope == 'NTT':
            observatory = EarthLocation.of_site('La Silla Observatory')
        elif telescope == 'TNT':
            observatory = EarthLocation.from_geodetic(
                '98 29 12','18 35 26',2457
            )
        else:
            raise ValueError('did not recognise telescope =',telescope)

        # read and store the hand written log
        handlog = os.path.join(do_night, f"{do_night}.dat")
        hlog = Log(handlog)

        # read target info from standard locations
        targets = Targets('TARGETS', 'AUTO_TARGETS')
        skip_targets, failed_targets = load_skip_fail()

        # Just re-do timing and positions for a particular night.
        runs = [run[:-4] for run in os.listdir(do_night) if fre.match(run)]
        runs.sort()
        if len(runs) == 0:
            print(f'Found no runs in night = {do_night}')
            return
        else:
            print(f'Found {len(runs)} runs in night = {do_night}\n')

        # create directory for any meta info such as the times
        meta = os.path.join(do_night, 'meta')
        os.makedirs(meta, exist_ok=True)

        # make the times
        times = os.path.join(meta, 'times')
        tdata = make_times(do_night, runs, observatory, times, True)
        print(f'Created & wrote timing data for {do_night} to {times}\n')

        # make the positions
        posdata = os.path.join(meta, 'posdata')
        make_positions(
            do_night, runs, observatory, hlog, targets,
            skip_targets, failed_targets, tdata, posdata, True
        )
        print(f'Created & wrote positional data for {do_night} to {posdata}')

        print(f'Finished creating time & position data for {do_night}')
        return

    barr = []
    cwd = os.getcwd()
    if os.path.basename(cwd) != "raw_data":
        print("** ulogger must be run in a directory called 'raw_data'")
        print("ulogger aborted")
        return

    if cwd.find("ultracam") > -1:
        instrument = "ULTRACAM"
        from hipercam.scripts.hmeta import ULTRACAM_META_COLNAMES
        COLNAMES = ULTRACAM_COLNAMES + ULTRACAM_META_COLNAMES[1:]
        nextra = len(ULTRACAM_META_COLNAMES)-1
    elif cwd.find("ultraspec") > -1:
        instrument = "ULTRASPEC"
        from hipercam.scripts.hmeta import ULTRASPEC_META_COLNAMES
        COLNAMES = ULTRASPEC_COLNAMES + ULTRASPEC_META_COLNAMES[1:]
        nextra = len(ULTRASPEC_META_COLNAMES)-1
    else:
        print("** ulogger: cannot find either ultracam or ultraspec in path")
        print("ulogger aborted")
        return
    linstrument = instrument.lower()

    # location to write files
    if os.path.exists(f'/storage/astro2/www/phsaap/{linstrument}/logs'):
        root = f'/storage/astro2/www/phsaap/{linstrument}/logs'
    else:
        root = 'logs'

    # Try an import now to get an early warning of problems
    import xlsxwriter

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
            "YYYY-MM with a file called 'telescope' in them"
        )
        print("ulogger aborted")
        return

    # If 'Others' exists, ensure it comes last.
    if 'Others' in rnames:
        rnames.remove('Others')
        rnames += ['Others']

    # Get started. First make sure all files are created with the
    # right permissions
    os.umask(0o022)

    # Ensure the root directory exists.
    os.makedirs(root, exist_ok=True)

    print(f'Will write to directory = "{root}".')

    # Load target positional information. Hand written, automatically
    # looked up, target names to skip and failed ones potentially
    # recoverable with some work
    targets = Targets('TARGETS', 'AUTO_TARGETS')
    skip_targets, failed_targets = load_skip_fail()

    # Index file
    index_tmp = os.path.join(root, 'index.html.tmp')
    index = os.path.join(root, 'index.html')
    with open(index_tmp, "w") as ihtml:
        # start the top level index html file
        ihtml.write(
            INDEX_HEADER.format(
                instrument=instrument,
                linstrument=linstrument
            )
        )

        for rname in rnames:
            print(f"\nProcessing run {rname}")

            with open(os.path.join(rname, "telescope")) as tel:
                telescope = tel.read().strip()

            # write in run date, start table of nights
            rn = os.path.basename(rname)
            if rn == 'Others':
                # Write "Others" to separate section at the bottom to avoid overwhelming the table
                ihtml.write(
                    f"""
</table>

<hr>
<h2>{rn}, {telescope}</h2>

<p>
""")
            else:
                try:
                    year, month = rn.split("-")
                    if month in LOG_MONTHS:
                        ihtml.write(
                            f"<tr><td>{LOG_MONTHS[month]} {year}</td><td>{telescope}</td><td>"
                        )
                    else:
                        ihtml.write(
                            f"<tr><td>{rn}</td><td>{telescope}</td><td>"
                        )
                except:
                    ihtml.write(
                        f"<tr><td>{rn}</td><td>{telescope}</td><td>"
                    )

            # set site
            if telescope == 'WHT':
                observatory = EarthLocation.of_site('Roque de los Muchachos')
            elif telescope == 'VLT':
                observatory = EarthLocation.of_site('Cerro Paranal')
            elif telescope == 'NTT':
                observatory = EarthLocation.of_site('La Silla Observatory')
            elif telescope == 'TNT':
                observatory = EarthLocation.from_geodetic(
                    '98 29 12','18 35 26',2457
                )
            else:
                raise ValueError('did not recognise telescope =',telescope)

            # get night directories which must have at least one run to be included
            nnames = []
            for ndir in os.listdir(rname):
                if nre.match(ndir):
                    if os.path.isdir(ndir):
                        nrun = len([run for run in os.listdir(ndir) if fre.match(run)])
                        if nrun == 0:
                            print(f'{ndir} contains no runs and will be ignored')
                        else:
                            nnames.append(ndir)
            nnames.sort()

            if len(nnames) == 0:
                print(
                    "found no night directories of the form YYYY-MM-DD"
                    f" in run = {rname}"
                )
                continue

            for nn, night in enumerate(nnames):

                # load all the run names
                runs = [run[:-4] for run in os.listdir(night) if fre.match(run)]
                runs.sort()
                if len(runs) == 0:
                    continue

                print(f"  night {night}")

                # create directory for any meta info such as the times
                meta = os.path.join(night, 'meta')

                os.makedirs(meta, exist_ok=True)
                links = '\n<p><a href="../index.html">Run index</a>'
                if nn > 0:
                    links += f', <a href="../{nnames[nn-1]}/">Previous night</a>'
                else:
                    links += f', Previous night'

                if nn < len(nnames) - 1:
                    links += f', <a href="../{nnames[nn+1]}/">Next night</a>\n</p>\n'
                else:
                    links += f', Next night\n</p>\n'

                links += "\n</p>\n"

                # Write an entry for each night linking to the log for
                # that night.
                if nn == 0:
                    ihtml.write(
                        f'<a href="{night}/">{night}</a>'
                    )
                    old_year = night[:4]
                else:
                    if rname == 'Others' and night[:4] != old_year:
                        ihtml.write(
                            f'<br><br>\n<a href="{night}/">{night}</a>'
                        )
                        old_year = night[:4]
                    else:
                        ihtml.write(
                            f', <a href="{night}/">{night}</a>'
                        )

                # Create the directory for the night
                date = f"{night}, {telescope}"
                ndir = os.path.join(root, night)
                os.makedirs(ndir, exist_ok=True)

                # and now the index file
                fname = os.path.join(ndir, 'index.html')
                fname_tmp = os.path.join(ndir, 'index.html.tmp')

                # see if there is a file of stats from hmeta avaialable.
                fstats = os.path.join(meta, 'statistics.csv')
                if os.path.exists(fstats):
                    stats = pd.read_csv(fstats,index_col='run_no')
                else:
                    stats = None

                with open(fname_tmp, "w") as nhtml:

                    # write header of night file
                    nhtml.write(NIGHT_HEADER1)
                    nhtml.write(NIGHT_HEADER2.format(date=date,instrument=instrument))
                    nhtml.write(links)
                    nhtml.write(TABLE_TOP)

                    # read and store the hand written log
                    handlog = os.path.join(night, f"{night}.dat")
                    hlog = Log(handlog)

                    ####################################
                    #
                    # Get or create timing info

                    times = os.path.join(meta, 'times')
                    if not do_full and os.path.exists(times):
                        # pre-existing file found
                        tdata = {}
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
                        tdata = make_times(night, runs, observatory, times, False)

                    ##################################
                    # Get or create positional info

                    posdata = os.path.join(meta, 'posdata')
                    pdata = {}
                    if not do_full and not do_positions and os.path.exists(posdata):
                        # pre-existing file found
                        with open(posdata) as pin:
                            for line in pin:
                                arr = line.split()
                                arr[3] = arr[3].replace('~',' ')
                                pdata[arr[0]] = [
                                    '' if val == 'UNDEF' else val for val in arr[1:]
                                ]
                        print('Read position data from',posdata)

                    else:
                        # create it
                        pdata = make_positions(
                            night, runs, observatory, hlog, targets, skip_targets,
                            failed_targets, tdata, posdata, False, rname
                        )

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
                            nhtml.write(
                                ULTRACAM_TABLE_HEADER if instrument == 'ULTRACAM' else
                                ULTRASPEC_TABLE_HEADER
                            )

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
                                exc_traceback, limit=1, file=sys.stdout
                            )
                            traceback.print_exc(file=sys.stdout)
                            print("Problem on run = ", runname)

                            # dummy info line just to allow us to proceed
                            nhtml.write("<tr>\n")
                            # run number
                            nhtml.write(f'<td class="lalert">{run}</td>')
                            nhtml.write("</tr>\n")
                            if instrument == 'ULTRACAM':
                                brow = [night, run[3:]] + (49+nextra)*[None]
                            else:
                                brow = [night, run[3:]] + (55+nextra)*[None]
                            continue

                        hd = rhead.header

                        # start the row
                        nhtml.write("<tr>\n")

                        # run number
                        png = os.path.join(meta,f'{run}.png')
                        if os.path.exists(png):
                            npng = os.path.join(ndir, f'{run}.png')
                            nhtml.write(f'<td class="left"><a href="{run}.png">{run}</a></td>')
                            shutil.copyfile(png, npng)
                        else:
                            nhtml.write(f'<td class="left">{run}</td>')

                        # start list to append to array for spreadsheet
                        brow = [night, run,]

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
                            rastr, decstr, ra, dec = ['','',None,None]

                        nhtml.write(f'<td class="left">{target}</td><td class="left">{autoid}</td>')
                        nhtml.write(f'<td class="left">{rastr}</td><td class="left">{decstr}</td>')

                        if instrument == 'ULTRASPEC':
                            # instr PA
                            tel_ra = hd.get("RA", "")
                            tel_dec = hd.get("Dec", "")
                            tel_pa = hd.get("PA", "")
                            nhtml.write(f'<td class="cen">{tel_ra}</td><td class="cen">{tel_dec}</td><td class="cen">{tel_pa}</td>')
                            if tel_ra != '' and tel_dec != '':
                                try:
                                    tel_ra_deg, tel_dec_deg, syst = str2radec(tel_ra + ' ' + tel_dec)
                                    tel_ra_deg = round(15*tel_ra_deg,5)
                                    tel_dec_deg = round(tel_dec_deg,4)
                                except:
                                    exc_type, exc_value, exc_traceback = sys.exc_info()
                                    traceback.print_tb(
                                        exc_traceback, limit=1, file=sys.stdout
                                    )
                                    traceback.print_exc(file=sys.stdout)
                                    tel_ra_deg, tel_dec_deg = None, None
                                    print("Position problem on run = ", runname)
                            else:
                                tel_ra_deg, tel_dec_deg = None, None
                            brow += [target, autoid, rastr, decstr, ra, dec, tel_ra, tel_dec, tel_ra_deg, tel_dec_deg, noval(tel_pa), rname]
                        else:
                            brow += [target, autoid, rastr, decstr, ra, dec, rname]

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
                            brow += ['','']
                            mjd_start = None

                        try:
                            # End time, total time, cadence (assumes we have a valid start time)
                            mjd_end = float(mjd_end)
                            ttime = round(86400.*(mjd_end - mjd_start))
                            nhtml.write(
                                f'<td class="cen">{ut_end}</td> <td class="right">{ttime}</td> <td class="right">{cadence}</td>'
                            )
                            brow += [ut_end, mjd_start, mjd_end, ttime, noval(cadence)]
                        except:
                            nhtml.write(3*'<td></td>')
                            brow += ['',mjd_start,None,None,None]

                        # Actual exposure time per frame
                        nhtml.write(f'<td class="right">{expose}</td>')
                        brow.append(noval(expose))

                        # number of frames, number ok
                        nhtml.write(f'<td class="right">{ntotal}</td><td class="right">{nok}</td>')
                        brow += [noval(ntotal), noval(nok)]

                        # filters used
                        if hlog.format == 1:
                            filters = hlog.filters.get(run, '')
                        else:
                            filters = hd.get("filters", '')
                        nhtml.write(f'<td class="cen">{filters}</td>')
                        brow.append(filters)

                        # run type
                        itype = hd.get("DTYPE", '')
                        nhtml.write(f'<td class="left">{itype}</td>')
                        brow.append(itype)

                        # readout mode
                        nhtml.write(f'<td class="cen">{rhead.mode}</td>')
                        brow.append(rhead.mode)

                        # nblue
                        if instrument == 'ULTRACAM':
                            nblue = hd['NBLUE']
                            nhtml.write(f'<td class="cen">{nblue}</td>')
                            brow.append(nblue)

                        # window formats
                        nwmax = 3 if instrument == 'ULTRACAM' else 4
                        for n in range(nwmax):
                            if n < len(rhead.wforms):
                                win = rhead.wforms[n]
                                nhtml.write(f'<td class="cen">{win}</td>')
                                brow.append(win)
                            else:
                                nhtml.write(f'<td></td>')
                                brow.append('')

                        # binning
                        binning = f'{rhead.xbin}x{rhead.ybin}'
                        nhtml.write(f'<td class="cen">{binning}</td>')
                        brow.append(binning)

                        # clear
                        clear = "Y" if hd['CLEAR'] else "N"
                        nhtml.write(f'<td class="cen">{clear}</td>')
                        brow.append(clear)

                        # CCD speed
                        speed = hd['GAINSPED'] if instrument == 'ULTRACAM' else \
                            hd['SPEED']
                        nhtml.write(f'<td class="cen">{speed}</td>')
                        brow.append(speed)

                        if instrument == 'ULTRASPEC':
                            output = hd.get('OUTPUT','')
                            nhtml.write(f'<td class="cen">{output}</td>')
                            brow.append(output)

                            hv_gain = hd.get('HV_GAIN','')
                            nhtml.write(f'<td class="cen">{hv_gain}</td>')
                            brow.append(hv_gain)

                        # Focal plane slide
                        fpslide = hd.get('SLIDEPOS','UNKNOWN')
                        try:
                            fpslide = float(fpslide)
                            fpslide = round(fpslide,1)
                        except:
                            fpslide = ''
                        nhtml.write(f'<td class="cen">{fpslide}</td>')
                        brow.append(noval(fpslide))

                        # Observers
                        observers = hd.get("OBSERVRS", '')
                        nhtml.write(f'<td class="cen">{observers}</td>')
                        brow.append(observers)

                        # PI
                        pi = hd.get("PI", '')
                        nhtml.write(f'<td class="left">{pi}</td>')
                        brow.append(pi)

                        # Program ID
                        pid = hd.get("ID", "")
                        nhtml.write(f'<td class="left">{pid}</td>')
                        brow.append(pid)

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
                        nhtml.write(f'<td class="left">{run}</td>')
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
                        brow += [noval(v) for v in pdat[3:]]

                        # at last: end the row
                        nhtml.write("\n</tr>\n")

                        try:
                            # add in statistics
                            brow += stats.loc[run].values.tolist()
                        except:
                            brow += nextra*[None]

                        if len(brow) != len(COLNAMES):
                            print(f'{runname}: data items vs colnames = {len(brow)} vs {len(COLNAMES)}')

                        barr.append(brow)

                    # finish off the night file
                    nhtml.write("</table>\n{:s}".format(links))
                    nhtml.write(NIGHT_FOOTER)

                # rename night file
                shutil.move(fname_tmp, fname)

            ihtml.write('</td></tr>\n')

        if 'Others' not in rnames:
            ihtml.write('</table>\n')

        # finish the main index
        ihtml.write(INDEX_FOOTER)

    # rename index file (means that any old index
    # file still exists while new one is being written)
    shutil.move(index_tmp, index)

    # write out the css file
    css = os.path.join(root, "ultra.css")
    with open(css, "w") as fout:
        fout.write(LOG_CSS)

    # write out page of info about sql database
    sqdb = os.path.join(root, 'sqldb.html')
    with open(sqdb, "w") as sqout:
        sqout.write(f'''<html>
<head>

<title>{instrument} sqlite3 database</title>
</head>

<body>
<h1>{instrument} sqlite3 database</h1>

<p>
The sqlite3 database file <a href="{linstrument}.db">{linstrument}.db</a> is designed to allow SQL queries
to be applied using Python's sqlite3 module. This has a single table called "{linstrument}". The columns of
the database are defined below. Note although some are naturally integers, they are converted to floats to
allow null values due to details of pandas.

<p>
Here is a simple example of carrying out a search for ULTRACAM. The position selected matches AR Sco and in this case
only runs longer than 200 seconds are returned (37 runs in the ULTRACAM database as of Feb 2021).
<pre>
import sqlite3
import pandas as pd

# Connect to the database
cnx = sqlite3.connect('{linstrument}.db')

# Query string
query = 'select * from ultracam WHERE ra_deg > 245 and ra_deg < 246 and dec_deg > -23 and dec_deg < -22 and total > 200'

# return the result as a pandas DataFrame
res = pd.read_sql_query(query, cnx)
print(res)
cnx.close()
</pre>

<p>
Here is a more complex example in which we carry out the same simple search as before
to derive a sub-table labelled t2, but then search through the original table (t1) for
all runs taken within 2 days of the t2 ones which match them in terms of format. This returns
matching calibration frames, checks that the filters are right for flats and the red speed
is OK for biases. It also writes the results to disk.

<pre>
import sqlite3
import pandas as pd
from hipercam.utils import format_ulogger_table

# Connect to the database
cnx = sqlite3.connect('ultracam.db')

query = """
SELECT t1.*
FROM ultracam AS t1, (
   SELECT *
   FROM ultracam
   WHERE ra_deg > 245 AND ra_deg < 246 AND dec_deg > -23 AND dec_deg < -22 AND total > 200
) as t2
WHERE ABS(t1.mjd_start-t2.mjd_start) < 2
AND (t1.run_type != 'flat' OR t1.filters == t2.filters)
AND (t1.run_type != 'bias' OR t1.read_speed == t2.read_speed)
AND (
     (t1.win1 == t2.win1 AND t1.binning == t2.binning AND ((t1.ra_hms == t2.ra_hms) OR t1.ra_deg is NULL))
     OR
     (t1.win1 == '1,513,1,512,1024' AND t1.binning == '1x1') AND ((t1.ra_hms == t2.ra_hms) OR t1.ra_deg is NULL)
    )
GROUP BY t1.night, t1.run_no
"""

res = pd.read_sql_query(query, cnx)
cnx.close()
print(res)
format_ulogger_table('results.xlsx', res, 'ultracam')
</pre>

<p>
This search comes up with 283 runs as of 22 Feb; there tend to be a fair few matching calibrations for
each bit of real data, but this does cut down significantly on what to look through.

<h2>Column names, definitions and data types</h2>

<p>
The database is called {linstrument}.db and contains a single table called '{linstrument}' with the following columns:

<p>
<table>
<tr><th class="left">Column name</th><th class="left">Data type</th><th class="left">Definition</th></tr>
''')

        dtypes = {}
        cnames = []
        for cname, dtype, definition in COLNAMES:
            sqout.write(
                f'<tr><td class="left">{cname}</td><td>{dtype}</td><td class="left">{definition}</td></tr>\n'
            )
            cnames.append(cname)
            dtypes[cname] = dtype

        sqout.write("""</table>

<hr>
<address>Tom Marsh, Warwick</address>
</body>
</html>
""")

    print('\nFinished generation of web pages; now generating a spreadsheet and an SQL database')

    # create pd.DataFrame containing all info
    ptable = pd.DataFrame(data=barr,columns=cnames)

    # enforce data types
    ptable = ptable.astype(dtypes)

    spreadsheet = os.path.join(root, f"{linstrument}-log.xlsx")
    format_ulogger_table(spreadsheet, ptable, linstrument)
    print(f'Written spreadsheet to {linstrument}-log.xlsx')

    # write out sqlite database
    sqldb = os.path.join(root, f'{linstrument}.db')
    cnx = sqlite3.connect(sqldb)
    ptable.to_sql(name=f'{linstrument}', con=cnx, if_exists='replace')
    cnx.commit()
    cnx.close()
    print(f'Written sqlite database to {linstrument}.db')
    print(f'Table dimensions (rows,columns) = {ptable.shape}')
    print(f'\nAll done. Look in {root} for the outputs.')

# End of main section

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
        except FileNotFoundError:
            sys.stdout.write(f'Night log = {fname} does not exist\n')
        except Exception as err:
           sys.stdout.write(f'Problem on night log = {fname}:' + str(err) + '\n')


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

    def add_target(self, target, ra, dec, uid, dmax=10):
        """Add a target at position ra, dec decimal hours and degrees. dmax
        is maximum separation in arcsec from previous entry with
        matching lookup name. 'target' is a name that will be matched
        to a unique identifier uid

        """

        if target in self.lnames:
            raise hcam.HipercamError(
                'Target = {target} name = {name} already has an entry'
            )

        if uid in self:
            entry = self[uid]
            ora, odec = entry['ra'], entry['dec']
            dist = np.sqrt((15*(ra-ora))**2+(np.cos(np.radians(dec))*(dec-odec))**2)
            if 3600*dist > dmax:
                raise hcam.HipercamError(
                    f'Target id = {uid} already entered but the new position is {3600*dist}" (>{dmax}) way from first position'
                )

            # add target to the list of names associated with ID
            self[uid]['names'].append(target)

        else:
            # new entry
            self[uid] = {'ra' : ra, 'dec' : dec, 'names' : [target]}

        # ensure target maps to the id
        self.lnames[target] = uid

    def __call__(self, target):
        """
        Call with name to lookup. Returns standardised name plus
        dictionary with 'ra', 'dec' and 'names'
        """
        if target in self:
            uid = target
        else:
            uid = self.lnames[target]
        dct = self[uid]
        return uid, dct['ra'], dct['dec']

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
<table border=1 cellspacing="1" cellpadding="4">
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


class TimeHMSCustom(TimeISO):
    """
    Just the HMS part of a Time as "<HH>:<MM>:<SS.sss...>".
    """

    name = "hms_custom"  # Unique format name
    subfmts = (("date_hms", "%H%M%S", "{hour:02d}:{min:02d}:{sec:02d}"),)


def make_times(night, runs, observatory, times, full):
    """
    Generates timing data for a set of runs from a particular night.
    Results are stored in a file called "times", and returned as
    dictionary keyed on run numbers.
    """

    # use this to check times are vaguely right. time of runs
    # must lie between 06.00 local time on date corresponding to
    # start of night date and 1.5 days later
    mjd_ref = Time(night).mjd - observatory.lon.degree/360 + 0.25

    tdata = {}
    with open(times,'w') as tout:
        for run in runs:
            if full:
                print(f'Analysing times for run {run}')
            dfile = os.path.join(night, run)
            try:
                ntotal = 0
                rtime = hcam.ucam.Rtime(dfile)

                # Find first good time, has to roughly match the start
                # date of the night because some times can just be
                # junk
                not_alerted = True
                for n, tdat in enumerate(rtime):
                    time, tinfo = tdat[:2]
                    if time.good:
                        mjd_start = time.mjd
                        tdelta = mjd_start-mjd_ref
                        if tdelta > 0 and tdelta < 1.5:
                            ts = Time(mjd_start, format="mjd", precision=2)
                            ut_start = ts.hms_custom
                            expose = round(time.expose,3)
                            n_start = n+1
                            if expose > 0 and expose < 500:
                                break
                        elif not_alerted and (tdelta < 0 or tdelta > 1.5):
                            # maximum one warning per run
                            not_alerted = False
                            print(f'  Bad time: tdelta = {tdelta} < 0 or > 1.5 on time {n} of {dfile}')
                else:
                    ntotal = 0
                    raise hcam.HipercamError(f'No good times found in {dfile}')

                # Find last good time. First we just go for times near the
                # end of the run. Failing that, we try again from the start,
                # to account for runs with time stamp issues.
                if rtime.header['MODE'] == 'DRIFT':
                    # ultracam
                    win = rtime.win[0]
                    nyu = win.ny*rtime.ybin
                    nback = int((1033/nyu + 1) / 2) + 3
                elif rtime.header['MODE'] == 'UDRIFT':
                    # ultraspec
                    win = rtime.win[0]
                    nyu = win.ny*rtime.ybin
                    nback = int((1037/nyu + 1) / 2) + 3
                else:
                    # non drift mode
                    nback = 4

                nbytes = os.stat(dfile + '.dat').st_size
                ntotal = nbytes // rtime.framesize

                if ntotal > 20000:
                    # this is a risk-reducing
                    # strategy in case the end of
                    # a long run is
                    # corrupt. Better to look at
                    # more than the necessary
                    # number of frames if it
                    # prevents us from having to
                    # wind through the whole lot.
                    nback = max(nback, 500)

                nreset = max(1, ntotal - nback)
                rtime.set(nreset)

                flast = False
                for n, tdat in enumerate(rtime):
                    time, tinfo = tdat[:2]
                    if time.good:
                        mjd = time.mjd
                        if mjd >= mjd_start and mjd < mjd_start + 0.4:
                            mjd_end = mjd
                            ts = Time(mjd_end, format="mjd", precision=2)
                            ut_end = ts.hms_custom
                            n_end = nreset + n
                            nexpose = round(time.expose,3)
                            if nexpose < 500:
                                expose = max(expose, nexpose)
                                flast = True

                if not flast:
                    # no good time found near end. There must be
                    # one or we wouldn't get to this point, so
                    # grind it out the hard way by going through
                    # the whole run, which can be slow.
                    rtime.set()
                    for n, tdat in enumerate(rtime):
                        time, tinfo = tdat[:2]
                        if time.good:
                            mjd = time.mjd
                            if mjd >= mjd_start and mjd < mjd_start + 0.4:
                                mjd_end = mjd
                                ts = Time(mjd_end, format="mjd", precision=2)
                                ut_end = ts.hms_custom
                                n_end = n + 1
                                nexpose = round(time.expose,3)
                                if nexpose < 500:
                                    expose = max(expose, nexpose)

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
                if full: print(f'{run} was a power-on or -off')

            except hcam.HipercamError:
                # No good times
                tdata[run] = ['','','','','','',0,ntotal]
                tout.write(f'{run} UNDEF UNDEF UNDEF UNDEF UNDEF UNDEF 0 {ntotal}\n')
                if full:
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_tb(exc_traceback, limit=1)
                    traceback.print_exc()
                    print(f'No good times found for {run}; ntotal = {ntotal}')

            except:
                # some other failure
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_tb(exc_traceback, limit=1)
                traceback.print_exc()
                print("Problem on run = ", dfile)

                # Load of undefined
                tdata[run] = 8*['']
                tout.write(f'{run} {" ".join(8*["UNDEF"])}\n')

    print('Written timing data to',times)
    return tdata

def make_positions(
        night, runs, observatory, hlog, targets,
        skip_targets, failed_targets, tdata, posdata, full,
        rname=None
):
    """
    Determine positional info, write to podata,
    return as dictionary keyed on the runs. Uses pre-determined
    timing data from make_times
    """
    from astroplan import moon_phase_angle

    pdata = {}
    with open(posdata,'w') as pout:
        for run in runs:
            if len(tdata[run]) == 1:
                # means its a power on/off
                continue

            if full:
                print(f'Analysing positions for run {run}')

            # open the run file as an Rhead
            runname = os.path.join(night, run)
            try:
                rhead = hcam.ucam.Rhead(runname)
            except:
                if full:
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_tb(exc_traceback, limit=1)
                    traceback.print_exc()
                    print(f"Failed top open {runname} as an Rhead")
                continue

            # object name
            if hlog.format == 1:
                target = hlog.target[run]
            else:
                target = rhead.header.get("TARGET",'')
            target = target.strip()

            # RA, Dec lookup
            if target == '' or target in skip_targets or target in failed_targets:
                autoid, ra, dec = 'UNDEF', 'UNDEF', 'UNDEF'
            else:
                try:
                    autoid, ra, dec = targets(target)
                except:
                    try:
                        # attempt simbad lookup here
                        autoid, ra, dec = target_lookup(target)
                        targets.add_target(target, ra, dec, autoid)
                        print(f'  Added {target} to targets')
                    except:
                        print(f'  No position found for {runname}, target = "{target}"')
                        if rname is not None:
                            failed_targets[target] = (rname,night,run)
                        autoid, ra, dec = 'UNDEF', 'UNDEF', 'UNDEF'

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
                    arr = arr[:3] + 8*['UNDEF']

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
                arr = arr[:3] + 13*['UNDEF']

            autoid_nospace = arr[2].replace(' ','~')
            pout.write(
                f'{run} {arr[0]} {arr[1]} {autoid_nospace} {arr[3]} ' +
                f'{arr[4]} {arr[5]} {arr[6]} {arr[7]} {arr[8]} {arr[9]} ' +
                f'{arr[10]} {arr[11]} {arr[12]} {arr[13]} {arr[14]} {arr[15]}\n'
            )

            arr[2] = arr[2].replace('~',' ')
            pdata[run] = [
                '' if val == 'UNDEF' else val for val in arr
            ]

    print('Written positional data to',posdata)
    return pdata

# Create and write out spreadsheet
ULTRACAM_COLNAMES = (
    ('night' , 'str', 'Date at the start of the night'),
    ('run_no', 'str', 'Run number'),
    ('target', 'str', 'Target name'),
    ('auto_id', 'str', 'Lookup name, from disk files or SIMBAD'),
    ('ra_hms', 'str', 'Target RA (J2000), HMS'),
    ('dec_dms', 'str', 'Target Dec (J2000), DMS'),
    ('ra_deg', 'float64', 'Target RA (J2000), degrees'),
    ('dec_deg', 'float64', 'Target Dec (J2000), degrees'),
    ('obs_run', 'str', 'Run group title'),
    ('date_start', 'str', 'Date at the start of the run'),
    ('utc_start', 'str', 'UTC at the start of the run'),
    ('utc_end', 'str', 'UTC at the end of the run'),
    ('mjd_start', 'float64', 'MJD (UTC) at the start of the run'),
    ('mjd_end', 'float64', 'MJD (UTC) at the end of the run'),
    ('total', 'float32', 'Total exposure time (seconds)'),
    ('cadence', 'float32', 'Sampling time between exposures (seconds)'),
    ('exposure', 'float32', 'Actual exposure time (seconds)'),
    ('nframe', 'float32', 'Total number of frames in file'),
    ('nok', 'float32', 'Number of frames from first to last with OK times'),
    ('filters', 'str', 'Colour filters used'),
    ('run_type', 'str', 'Provisional type of the run, indicative only'),
    ('read_mode', 'str', 'ULTRACAM readout mode'),
    ('nb', 'int', 'Nblue, CCD 1 readout skip'),
    ('win1', 'str', 'Format of window pair 1'),
    ('win2', 'str', 'Format of window pair 2'),
    ('win3', 'str', 'Format of window pair 3'),
    ('binning', 'str', 'X by Y binning'),
    ('clr', 'str', 'Clear enabled or not'),
    ('read_speed', 'str', 'Readout speed'),
    ('fpslide', 'float32', 'Focal plane slide'),
    ('observers', 'str', 'Observers'),
    ('pi', 'str', 'PI of data'),
    ('pid', 'str', 'Proposal ID'),
    ('tel', 'str', 'Telescope'),
    ('size', 'float32', 'Data file size, MB'),
    ('nlink', 'str', 'Link to html log of night containing run (spreadsheet only)'),
    ('comment', 'str', 'Hand comments from night'),
    ('alt_start', 'float32', 'Target altitude at start of run, degrees'),
    ('alt_middle', 'float32', 'Target altitude in middle of run, degrees'),
    ('alt_end', 'float32', 'Target altitude at end of run, degrees'),
    ('az_start', 'float32', 'Target azimuth at start of run, degrees'),
    ('az_middle', 'float32', 'Target azimuth in middle of run, degrees'),
    ('az_end', 'float32', 'Target azimuth at end of run, degrees'),
    ('sun_dist', 'float32', 'Distance from Sun, degrees, middle of run'),
    ('moon_dist', 'float32', 'Distance from Moon, degrees, middle of run'),
    ('sun_alt_start', 'float32', 'Altitude of Sun at start of run'),
    ('sun_alt_end', 'float32', 'Altitude of Sun at end of run'),
    ('moon_alt_start', 'float32', 'Altitude of Moon at start of run'),
    ('moon_alt_end', 'float32', 'Altitude of Moon at end of run'),
    ('moon_phase', 'float32', 'Lunar phase (fraction illuminated')
)

# Create and write out spreadsheet
ULTRASPEC_COLNAMES = (
    ('night' , 'str', 'Date at the start of the night'),
    ('run_no', 'str', 'Run number'),
    ('target', 'str', 'Target name'),
    ('auto_id', 'str', 'Lookup name, from disk files or SIMBAD'),
    ('ra_hms', 'str', 'Target RA (J2000), HMS'),
    ('dec_dms', 'str', 'Target Dec (J2000), DMS'),
    ('ra_deg', 'float64', 'Target RA (J2000), degrees'),
    ('dec_deg', 'float64', 'Target Dec (J2000), degrees'),
    ('ra_tel', 'str', 'Telescope pointing RA (J2000), HMS'),
    ('dec_tel', 'str', 'Telescope pointing Dec (J2000), DMS'),
    ('ra_tel_deg', 'float64', 'Telescope pointing RA (J2000), degrees'),
    ('dec_tel_deg', 'float64', 'Telescope pointing Dec (J2000), degrees'),
    ('pa', 'float32', 'PA on sky (degrees)'),
    ('obs_run', 'str', 'Run group title'),
    ('date_start', 'str', 'Date at the start of the run'),
    ('utc_start', 'str', 'UTC at the start of the run'),
    ('utc_end', 'str', 'UTC at the end of the run'),
    ('mjd_start', 'float64', 'MJD (UTC) at the start of the run'),
    ('mjd_end', 'float64', 'MJD (UTC) at the end of the run'),
    ('total', 'float32', 'Total exposure time (seconds)'),
    ('cadence', 'float32', 'Sampling time between exposures (seconds)'),
    ('exposure', 'float32', 'Actual exposure time (seconds)'),
    ('nframe', 'float32', 'Total number of frames in file'),
    ('nok', 'float32', 'Number of frames from first to last with OK times'),
    ('filters', 'str', 'Colour filters used'),
    ('run_type', 'str', 'Provisional type of the run, indicative only'),
    ('read_mode', 'str', 'ULTRASPEC readout mode'),
    ('win1', 'str', 'Format of window 1'),
    ('win2', 'str', 'Format of window 2'),
    ('win3', 'str', 'Format of window 3'),
    ('win4', 'str', 'Format of window 4'),
    ('binning', 'str', 'X by Y binning'),
    ('clr', 'str', 'Clear enabled or not'),
    ('read_speed', 'str', 'Readout speed'),
    ('output', 'str', 'amplifier output in use (0 = normal, 1 = avalanche)'),
    ('hv_gain', 'str', 'High voltage gain setting'),
    ('fpslide', 'float32', 'Focal plane slide'),
    ('observers', 'str', 'Observers'),
    ('pi', 'str', 'PI of data'),
    ('pid', 'str', 'Proposal ID'),
    ('tel', 'str', 'Telescope'),
    ('size', 'float32', 'Data file size, MB'),
    ('nlink', 'str', 'Link to html log of night containing run (spreadsheet only)'),
    ('comment', 'str', 'Hand comments from night'),
    ('alt_start', 'float32', 'Target altitude at start of run, degrees'),
    ('alt_middle', 'float32', 'Target altitude in middle of run, degrees'),
    ('alt_end', 'float32', 'Target altitude at end of run, degrees'),
    ('az_start', 'float32', 'Target azimuth at start of run, degrees'),
    ('az_middle', 'float32', 'Target azimuth in middle of run, degrees'),
    ('az_end', 'float32', 'Target azimuth at end of run, degrees'),
    ('sun_dist', 'float32', 'Distance from Sun, degrees, middle of run'),
    ('moon_dist', 'float32', 'Distance from Moon, degrees, middle of run'),
    ('sun_alt_start', 'float32', 'Altitude of Sun at start of run'),
    ('sun_alt_end', 'float32', 'Altitude of Sun at end of run'),
    ('moon_alt_start', 'float32', 'Altitude of Moon at start of run'),
    ('moon_alt_end', 'float32', 'Altitude of Moon at end of run'),
    ('moon_phase', 'float32', 'Lunar phase (fraction illuminated')
)


def noval(value):
    return None if value == '' else value


# Header of the main file

INDEX_HEADER = """<html>
<head>

<title>{instrument} data logs</title>
</head>

<body>
<h1>{instrument} data logs</h1>

<p> These on-line logs summarise all {instrument} runs. They were
automatically generated from the data files and hand-written logs.

<p>For a searchable file summarising the same information, plus
calculated data on the Sun and Moon for each run, see this <a
href="{linstrument}-log.xlsx">spreadsheet</a>. This is also
available as an <a href="sqldb.html">sqlite3 database</a> to
allow SQL queries.

<p>
<h2>Vikcam team observing nights</h2>

<p>
<table border=1 cellspacing="1" cellpadding="4">
<tr><th>Run</th><th>Telescope</th><th>Dates at the start of the nights</th></tr>
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

<link rel="stylesheet" type="text/css" href="../ultra.css" />
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

ULTRACAM_TABLE_HEADER = """
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
<td><button id="hidden28" onclick="hide(29)"></button></td>
<td><button id="hidden29" onclick="hide(29)"></button></td>
<td align="left"><button id="hidden30" onclick="hide(30)"></button></td>
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
<th class="cen">xl,xr,ys,nx,ny<br>(1)</th>
<th class="cen">xl,xr,ys,nx,ny<br>(2)</th>
<th class="cen">xl,xr,ys,nx,ny<br>(3)</th>
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

ULTRASPEC_TABLE_HEADER = """
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
<td><button id="hidden28" onclick="hide(29)"></button></td>
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
<th class="left">Auto<br>ID</th>
<th class="left">RA (J2000)</th>
<th class="left">Dec&nbsp;(J2000)</th>
<th class="left">Tel RA</th>
<th class="left">Tel Dec</th>
<th class="left">PA</th>
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
<th class="cen">xs,ys,nx,ny<br>(1)</th>
<th class="cen">xs,ys,nx,ny<br>(2)</th>
<th class="cen">xs,ys,nx,ny<br>(3)</th>
<th class="cen">xs,ys,nx,ny<br>(4)</th>
<th class="cen">XxY<br>bin</th>
<th class="cen">Clr</th>
<th class="cen">Read<br>speed</th>
<th class="cen">Output</th>
<th class="cen">HV Gain</th>
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

<p> 'PA' is the PA on the sky (ULTRASPEC). 'Clr' indicates whether
clears were enabled; 'Read mode' is the readout mode which can be one
of several options: 'FFCLR' for full frames with clear; 'FFNCLR' full
frames with no clear; '1-PAIR', '2-PAIR', '3-PAIR', for standard
windowed modes, etc. The 'xl,xr,..' column gives the 5 parameters
defining the window pair ULTRACAM, or 4 for ULTRASPEC. The 'cadence'
is the time between exposures, the 'exposure' the actual time spent
integrating, although note that these numbers are not always entirely
accurate. 'Nok' is the number of OK frames judged as having OK
times. If it grossly disagrees with Nframe, there might be timing
problems present, such as the GPS dropping out.
</p>

<address>Tom Marsh, Warwick</address>
</body>
</html>
"""



