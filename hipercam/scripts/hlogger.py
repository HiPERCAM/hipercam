import sys
import traceback
import os
import shutil
import re
import warnings
import argparse
import sqlite3
import json
import datetime

import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta, TimeISO
from astropy.coordinates import get_sun, get_moon, EarthLocation, SkyCoord, AltAz
import astropy.units as u

import hipercam as hcam
from hipercam.utils import format_hlogger_table, target_lookup, dec2sexg, str2radec, \
    LOG_CSS, LOG_MONTHS

__all__ = [
    "hlogger",
]

############################################################################
#
# hlogger -- generates html logs and spreadsheets for ultracam and ultraspec
#
############################################################################

###########
#
# Some constant stuff
#

def hlogger(args=None):
    description = \
    """``hlogger``

    Generates html logs for |hiper|, ULTRACAM and ULTRASPEC runs, a
    searchable spreadsheet and an sqlite3 database for programmatic
    SQL enquiries. The latter are used by the pipeline command |logsearch|.

    hlogger expects to work in a directory containing the runs for
    each night. I call this 'raw_data' for each camera and this will
    be checked as a guard against problems. It extracts information
    from each run file.

    hlogger tries to read target data from three files: TARGETS,
    FAILED_TARGETS and SKIP_TARGETS. TARGETS contains details of
    targets, particularly if not available on SIMBAD; FAILED_TARGETS
    are those of interest for tracking down; SKIP_TARGETS are ones of
    no interest.

    If run at Warwick, it writes data to the web pages. Otherwise it
    writes them to a sub-directory of raw_data called "logs". hlogger
    is a specialist routine this one; if you have access to the
    on-line logs at Warwick, it should be unnecessary for you to run
    it. It requires the installation of the python module xlsxwriter
    in order to write an excel spreadsheet; this is not a required
    module for the whole pipeline and may not exist. The script will
    fail early on if it is not installed.

    It also writes information to a sub-directory "meta" of each night
    directory, and can pick up information stored there by the related
    script |redplt| which should be run first.

    To save time, hlogger does not by default re-do everything, so it
    does not by default recreate the timing and position data files if
    they already exist. Even so, a default run takes a while and if
    you just want to see the logs for a few new nights, the '-q'
    option could help as that skips re-creating already existing logs
    and the spreadsheet and database (although it means the night-to-night
    links might not be quite right). Similarly the '-n' option is
    useful to simply focus on the timing and positional data of a
    single night.

    """

    warnings.filterwarnings("ignore")

    # Options
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        dest="full",
        action="store_true",
        help="carry out full re-computation of times, positional data, html logs, spreadsheet and SQL database",
    )
    parser.add_argument(
        "-n",
        dest="night", default=None,
        help="use this with a YYYY-MM-DD date to update times and positions for a specific night (lots of diagnostic ouput; no html created)",
    )
    parser.add_argument(
        "-p",
        dest="positions",
        action="store_true",
        help="re-compute positional data, html logs, spreadsheet and SQL database, but not the times",
    )
    parser.add_argument(
        "-q",
        dest="quick",
        action="store_true",
        help="quick update, useful for seeing the log of a new night. Skips olds logs or creation of the spreadsheet and SQL database",
    )
    parser.add_argument(
        "-r",
        dest="retry",
        action="store_true",
        help="re-try any missing positions loaded from a file as UNDEF",
    )
    parser.add_argument(
        "--nowrite",
        dest="nowrite",
        action="store_true",
        help="write nothing to disk (for speed)",
    )

    # Get command line options
    args = parser.parse_args()

    if (args.full or args.positions) and args.night is not None:
        print('-n switch is not compatible with either -f or -p; please check help')
        return

    # start with defining a few regular expressions to match run directories, nights,
    # and run files

    # runs of the form  '2021A', '2022B', '2020-21', '2021-P106', 'Others' supported
    rre = re.compile("^(Others|\d\d\d\d-\d\d|\d\d\d\d-P\d\d\d|\d\d\d\d[AB])$")

    # Night must be of form 2021-05-12 [YYYY-MM-DD]
    nre = re.compile("^\d\d\d\d-\d\d-\d\d$")

    # Files of form 'run123.xml' or 'run1234.fits'
    fre = re.compile("^(run\d\d\d\.xml|run\d\d\d\d\.fits)$")

    cwd = os.getcwd()
    if os.path.basename(cwd) != "raw_data":
        print("** hlogger must be run in a directory called 'raw_data'")
        print("hlogger aborted")
        return

    # POSITIONS are any positional data from phase II entries,
    # of use for coordinate lookup when simbad fails.
    POSITIONS = {}
    if cwd.find("ultracam") > -1:
        instrument = "ULTRACAM"
        COLNAMES = ULTRACAM_COLNAMES
        if 'ULTRACAM_PHASEII' in os.environ:
            with open(os.environ['ULTRACAM_PHASEII']) as fp:
                POSITIONS = json.load(fp)
            print(
                f'Loaded Phase II data from {os.environ["ULTRACAM_PHASEII"]}'
            )
        else:
            print('No phase II environment variable set')

    elif cwd.find("ultraspec") > -1:
        instrument = "ULTRASPEC"
        COLNAMES = ULTRASPEC_COLNAMES
        if 'ULTRASPEC_PHASEII' in os.environ:
            with open(os.environ['ULTRASPEC_PHASEII']) as fp:
                POSITIONS = json.load(fp)
            print(
                f'Loaded Phase II data from {os.environ["ULTRASPEC_PHASEII"]}'
            )
        else:
            print('No phase II environment variable set')

    elif cwd.find("hipercam") > -1:
        instrument = "HiPERCAM"
        COLNAMES = HIPERCAM_COLNAMES
        if 'HIPERCAM_PHASEII' in os.environ:
            with open(os.environ['HIPERCAM_PHASEII']) as fp:
                POSITIONS = json.load(fp)
            print(
                f'Loaded Phase II data from {os.environ["HIPERCAM_PHASEII"]}'
            )
        else:
            print('No phase II environment variable set')
    else:
        print("** hlogger: cannot find hipercam, ultracam or ultraspec in path")
        print("hlogger aborted")
        return
    linstrument = instrument.lower()
    print(f'Identified instrument as "{instrument}"')

    # create date-time string to use when writing files
    iso_time = datetime.datetime.utcnow().isoformat(timespec='seconds')

    # buffers of messages on target lookups to print at the end
    smessages, fmessages = [], []

    if args.night:

        # identify observatory
        tlink = os.path.join(args.night, "telescope")
        with open(os.path.join(args.night, "telescope")) as tel:
            telescope = tel.read().strip()

        if telescope == 'WHT' or telescope == 'GTC':
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

        tpath = os.readlink(tlink)
        rname = os.path.basename(os.path.dirname(tpath))
        # read and store the hand written log
        handlog = os.path.join(args.night, f"{args.night}.dat")
        hlog = Log(handlog)

        # read target info from standard locations
        targets = Targets('TARGETS')
        skip_targets = load_skip_fail()

        # Just re-do timing and positions for a particular night.
        runs = [os.path.splitext(run)[0] for run in os.listdir(args.night) if fre.match(run)]
        runs.sort()
        if len(runs) == 0:
            print(f'Found no runs in night = {args.night}')
            return
        else:
            print(f'Found {len(runs)} runs in night = {args.night}\n')

        # create directory for any meta info such as the times
        meta = os.path.join(args.night, 'meta')
        os.makedirs(meta, exist_ok=True)

        # make the times
        times = os.path.join(meta, 'times')
        tdata = make_times(args.night, runs, observatory, times, True, instrument, True)
        print(f'Created & wrote timing data for {args.night} to {times}\n')

        # make the positions
        posdata = os.path.join(meta, 'posdata')

        crname = CORRECTORS.get(rname,rname)
        p2positions = POSITIONS.get(crname,{})

        make_positions(
            args.night, runs, observatory, instrument, hlog, targets,
            skip_targets, tdata, posdata, False, None, True, rname,
            smessages, fmessages, p2positions, True
        )
        print(f'Created & wrote positional data for {args.night} to {posdata}')
        print(f'Finished creating time & position data for {args.night}')
        print('Note that the html log for this night has not been created or updated')

        if len(smessages):
            print('\nYou may want to add the following to TARGETS to short-circuit SIMBAD lookups:\n')
            print('\n'.join(smessages))

        if len(fmessages):
            print('\nYou may want to check these or add them to either SKIP_TARGETS or FAILED_TARGETS to short-circuit SIMBAD lookups:\n')
            print('\n'.join(fmessages))

        # finish specific night at this point
        return

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
            "there were no run directories of the forms "
            "'YYYY-MM', 'YYYY[A|B]', 'YYYY-PNNN' or 'Others'"
            " with a file called 'telescope' in them"
        )
        print("hlogger aborted")
        return

    # If 'Others' exists, ensure it comes last.
    if 'Others' in rnames:
        rnames.remove('Others')
        rnames += ['Others']

    # Get started. First make sure all files are created with the
    # right permissions
    os.umask(0o022)

    # shortcut
    okwrite = not args.nowrite

    # Ensure the root directory exists.
    if okwrite:
        os.makedirs(root, exist_ok=True)

        print(f'Will write to directory = "{root}".')

    # initialise storage list
    barr = []

    # Load target positional information. Hand written, automatically
    # looked up, target names to skip and failed ones potentially
    # recoverable with some work
    targets = Targets('TARGETS')
    skip_targets = load_skip_fail()

    # Index file. To avoid too much down time, write to a
    # a temporary that will be switched in at the end,
    if okwrite:
        index_tmp = os.path.join(root, 'index.html.tmp')
        index = os.path.join(root, 'index.html')
    else:
        index_tmp = os.devnull

    with open(index_tmp, "w") as ihtml:

        # start the top level index html file
        ihtml.write(
            INDEX_HEADER.format(
                instrument=instrument,
                linstrument=linstrument
            )
        )

        for rname in rnames:
            # loop through the runs
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
            if telescope == 'WHT' or telescope == 'GTC':
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

            # scan through the nights to work out which to re-do (time saving)
            redo = {}
            one_before_is_new = False
            for nn, night in enumerate(nnames):
                if args.full or args.positions or not args.quick :
                    redo[night] = True
                else:
                    already_there = os.path.exists(os.path.join(root, night, 'index.html'))
                    if one_before_is_new:
                        # Has to be re-done under any circumstances
                        # because of the previous/next links
                        redo[night] = True
                        one_before_is_new = not already_there
                    else:
                        redo[night] = one_before_is_new = not already_there
                        if one_before_is_new and nn > 0:
                            # Need to re-do one before because of previous/next links
                            redo[nnames[nn-1]] = True

            for nn, night in enumerate(nnames):

                # load all the run names without extensions
                runs = [os.path.splitext(run)[0] for run in os.listdir(night) if fre.match(run)]
                runs.sort()
                if len(runs) == 0:
                    continue

                # Write an entry in the main index for each night
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

                if not redo[night]:
                    # can save a lot time by not re-making the log
                    # file more often than not
                    print(f"  night {night} log already exists and will not be re-created")
                    continue

                print(f"  night {night}")

                # create directory for any meta info such as the times
                meta = os.path.join(night, 'meta')
                if okwrite:
                    os.makedirs(meta, exist_ok=True)

                # Add links to previous and next night of the run
                links = '<p><a href="../index.html">Run index</a>'
                if nn > 0:
                    links += f', <a href="../{nnames[nn-1]}/">Previous night</a>'
                else:
                    links += f', Previous night'

                if nn < len(nnames) - 1:
                    links += f', <a href="../{nnames[nn+1]}/">Next night</a></p>\n'
                else:
                    links += f', Next night</p>\n'

                # Create the directory for the night
                date = f"{night}, {telescope}"
                ndir = os.path.join(root, night)
                if okwrite:
                    os.makedirs(ndir, exist_ok=True)
                    fname = os.path.join(ndir, 'index.html')
                    fname_tmp = os.path.join(ndir, 'index.html.tmp')
                else:
                    fname_tmp = os.devnull

                with open(fname_tmp, "w") as nhtml:

                    # write header of night file
                    nhtml.write(
                        NIGHT_HEADER.format(
                            date=date,
                            instrument=instrument,
                            links=links
                        )
                    )

                    # Instrument specific table header
                    if linstrument == 'hipercam':
                        nhtml.write(HIPERCAM_TABLE_HEADER)
                    elif linstrument == 'ultracam':
                        nhtml.write(ULTRACAM_TABLE_HEADER)
                    elif linstrument == 'ultraspec':
                        nhtml.write(ULTRASPEC_TABLE_HEADER)

                    # read and store the hand written log
                    handlog = os.path.join(night, f"{night}.dat")
                    hlog = Log(handlog)

                    ####################################
                    #
                    # Get or create timing info

                    times = os.path.join(meta, 'times')
                    if not args.full and os.path.exists(times):
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
                        tdata = make_times(night, runs, observatory, times, False, instrument, okwrite)

                    ##################################
                    # Get or create positional info

                    posdata = os.path.join(meta, 'posdata')
                    load_old = not args.full and not args.positions

                    crname = CORRECTORS.get(rname,rname)
                    p2positions = POSITIONS.get(crname,{})

                    pdata = make_positions(
                        night, runs, observatory, instrument, hlog, targets, skip_targets,
                        tdata, posdata, load_old, args.retry, False, rname, smessages,
                        fmessages, p2positions, okwrite
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

                        if len(tdata[run]) == 1:
                            # means its a power on/off
                            continue

                        # open the run file as an Rhead or Rtime
                        runname = os.path.join(night, run)
                        try:
                            if instrument == 'HiPERCAM':
                                rthead = hcam.hcam.Rtime(runname)
                            else:
                                rthead = hcam.ucam.Rhead(runname)
                        except Exception as err:

                            # So many ways to fail while this takes so
                            # long to run, that we basically want to
                            # soldier on if we possibly can but print
                            # stuff out. But there are a few ULTRACAM
                            # files with invalid framesizes so we
                            # don't bother saying anything with these

                            if instrument == 'HiPERCAM' or str(err).find('Framesize')== -1:
                                exc_type, exc_value, exc_traceback = sys.exc_info()
                                traceback.print_tb(
                                    exc_traceback, limit=1, file=sys.stdout
                                )
                                traceback.print_exc(file=sys.stdout)
                                print(f"Problem on run = {runname}")
                            else:
                                print(f"Problem on run = {runname} [framesize]")

                            # dummy info line just to allow us to proceed
                            nhtml.write("<tr>\n")
                            # run number
                            nhtml.write(f'<td class="lalert">{run}</td>')
                            nhtml.write("</tr>\n")
                            if instrument == 'ULTRACAM':
                                brow = [rname, night, run[3:]] + 50*[None]
                            elif instrument == 'ULTRASPEC':
                                brow = [rname, night, run[3:]] + 57*[None]
                            elif instrument == 'HiPERCAM':
                                brow = [rname, night, run[3:]] + 62*[None]
                            continue

                        hd = rthead.header

                        # start the row
                        nhtml.write("<tr>\n")

                        # run number
                        png = os.path.join(meta,f'{run}.png')
                        if os.path.exists(png):
                            npng = os.path.join(ndir, f'{run}.png')
                            nhtml.write(
                                f'<td class="left"><a href="{run}.png">' +
                                f'{run}</a></td>'
                            )
                            shutil.copyfile(png, npng)
                        else:
                            nhtml.write(f'<td class="left">{run}</td>')

                        # start list to append to array for spreadsheet
                        brow = [rname, night, run,]

                        # object name
                        if hlog.format == 1:
                            target = hlog.target[run]
                        elif instrument == 'HiPERCAM':
                            target = hd.get("OBJECT",'')
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

                        nhtml.write(
                            f'<td class="left">{target}</td>' +
                            f'<td class="left">{autoid}</td>'
                        )
                        nhtml.write(
                            f'<td class="left">{rastr}' +
                            f'</td><td class="left">{decstr}</td>'
                        )

                        if instrument == 'HiPERCAM':
                            # instr PA
                            tel_ra = tradec(hd.get("RA", ""),False)
                            tel_dec = tradec(hd.get("Dec", ""),True)
                            tel_pa = hd.get("INSTRPA", "")
                            nhtml.write(
                                f'<td class="cen">{tel_ra}</td><td class="cen">' +
                                f'{tel_dec}</td><td class="cen">{tel_pa}</td>'
                            )
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
                            brow += [
                                target, autoid, rastr, decstr, ra, dec, tel_ra,
                                tel_dec, tel_ra_deg, tel_dec_deg, noval(tel_pa)
                            ]

                        elif instrument == 'ULTRASPEC':
                            # instr PA
                            tel_ra = hd.get("RA", "")
                            tel_dec = hd.get("Dec", "")
                            tel_pa = hd.get("PA", "")
                            nhtml.write(
                                f'<td class="cen">{tel_ra}</td><td class="cen">'
                                f'{tel_dec}</td><td class="cen">{tel_pa}</td>'
                            )
                            if tel_ra != '' and tel_dec != '':
                                try:
                                    tel_ra_deg, tel_dec_deg, syst = \
                                        str2radec(tel_ra + ' ' + tel_dec)
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
                            brow += [
                                target, autoid, rastr, decstr,
                                ra, dec, tel_ra, tel_dec,
                                tel_ra_deg, tel_dec_deg, noval(tel_pa)
                            ]
                        else:
                            brow += [target, autoid, rastr, decstr, ra, dec]

                        # Timing info
                        ut_start, mjd_start, ut_end, mjd_end, \
                            cadence, expose, nok, ntotal = tdata[run]

                        try:
                            # Start time, date
                            mjd_start = float(mjd_start)
                            tstart = Time(mjd_start, format='mjd').isot
                            date_start = tstart[:tstart.find("T")]
                            nhtml.write(
                                f'<td class="cen">{ut_start}</td>'
                            )
                            brow += [date_start, ut_start]
                        except:
                            nhtml.write('<td></td>')
                            brow += ['','']
                            mjd_start = None

                        try:
                            # End time, total time, cadence (assumes we have a valid start time)
                            mjd_end = float(mjd_end)
                            ttime = round(86400.*(mjd_end - mjd_start))
                            nhtml.write(
                                f'<td class="cen">{ut_end}</td> <td class="right">'
                                f'{ttime}</td> <td class="right">{cadence}</td>'
                            )
                            brow += [ut_end, mjd_start, mjd_end, ttime, noval(cadence)]
                        except:
                            nhtml.write(3*'<td></td>')
                            brow += ['',mjd_start,None,None,None]

                        # Actual exposure time per frame
                        nhtml.write(f'<td class="right">{expose}</td>')
                        brow.append(noval(expose))

                        # number of frames, number ok
                        nhtml.write(
                            f'<td class="right">{ntotal}</td>' +
                            f'<td class="right">{nok}</td>' +
                            f'<td class="right">{pdat[11]}</td>'
                        )
                        brow += [noval(ntotal), noval(nok)]

                        # filters used
                        if hlog.format == 1:
                            filters = hlog.filters.get(run, '')
                        else:
                            filters = hd.get("filters", '')

                        # shrink the "Super" names
                        filters = re.sub("Super ([ugriz])'", r"\1s'", filters)
                        nhtml.write(f'<td class="cen">{filters}</td>')
                        brow.append(filters)

                        # run type
                        if instrument == 'HiPERCAM':
                            itype = hd.get("IMAGETYP", '')
                        else:
                            itype = hd.get("DTYPE", '')
                        if itype == 'acquisition' or itype == 'data caution':
                            itype = 'acquire'
                        nhtml.write(f'<td class="left">{itype}</td>')
                        brow.append(itype)

                        # readout mode
                        if instrument == 'HiPERCAM':
                            mode = TRANSLATE_MODE[rthead.mode]
                            nhtml.write(f'<td class="cen">{mode}</td>')
                            brow.append(mode)

                            texps, toffs, nskips, tdead = rthead.tinfo()
                            skips = ",".join([str(nskip) for nskip in nskips])
                            nhtml.write(f'<td class="cen">{skips}</td>')
                            brow.append(skips)

                            # window formats
                            quad1 = rthead.wforms[0]
                            nhtml.write(f'<td class="cen">{quad1}</td>')
                            quad1 = quad1.replace('&nbsp;',' ')
                            brow.append(quad1)

                            quad2 = rthead.wforms[1] if len(rthead.wforms) > 1 else ""
                            nhtml.write(f'<td class="cen">{quad2}</td>')
                            quad2 = quad2.replace('&nbsp;',' ')
                            brow.append(quad2)

                        else:
                            mode = rthead.mode
                            nhtml.write(f'<td class="cen">{mode}</td>')
                            brow.append(mode)

                            # nblue
                            if instrument == 'ULTRACAM':
                                nblue = hd['NBLUE']
                                nhtml.write(f'<td class="cen">{nblue}</td>')
                                brow.append(nblue)

                            # window formats
                            nwmax = 3 if instrument == 'ULTRACAM' else 4
                            for n in range(nwmax):
                                if n < len(rthead.wforms):
                                    win = rthead.wforms[n]
                                    nhtml.write(f'<td class="cen">{win}</td>')
                                    brow.append(win)
                                else:
                                    nhtml.write(f'<td></td>')
                                    brow.append('')

                        # binning
                        binning = f'{rthead.xbin}x{rthead.ybin}'
                        nhtml.write(f'<td class="cen">{binning}</td>')
                        brow.append(binning)

                        if instrument == 'HiPERCAM':
                            # clear
                            clear = "Y" if rthead.clear else "N"
                            nhtml.write(f'<td class="cen">{clear}</td>')
                            brow.append(clear)

                            # dummy output in use
                            dummy = "Y" if rthead.dummy else "N"
                            nhtml.write(f'<td class="cen">{dummy}</td>')
                            brow.append(dummy)

                            # LED on
                            led = hd.get("ESO DET EXPLED", "----")
                            led = "Y" if led == 1 else "N"

                            nhtml.write(f'<td class="cen">{led}</td>')
                            brow.append(led)

                            # over-scan
                            oscan = "Y" if rthead.oscan else "N"
                            nhtml.write(f'<td class="cen">{oscan}</td>')
                            brow.append(oscan)

                            # pre-scan
                            pscan = "Y" if rthead.pscan else "N"
                            nhtml.write(f'<td class="cen">{pscan}</td>')
                            brow.append(pscan)

                            # Nodding
                            nod = hd.get("ESO DET SEQ1 TRIGGER", 0)
                            nod = "Y" if nod == 1 else "N"
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
                                fclock = "N"
                            elif fclock == 1:
                                fclock = "Y"
                            nhtml.write(f'<td class="cen">{fclock}</td>')
                            brow.append(fclock)

                            # Tbytes problem
                            tbytes = "OK" if rthead.ntbytes == 36 else "NOK"
                            nhtml.write(f'<td class="cen">{tbytes}</td>')
                            brow.append(tbytes)

                            # Focal plane slide
                            fpslide = hd.get("FPslide", "----")
                            nhtml.write(f'<td class="cen">{fpslide}</td>')
                            brow.append(fpslide)

                            # CCD temps
                            t1,t2,t3,t4,t5 = hd.get("CCD1TEMP", 0.0), \
                                hd.get("CCD2TEMP", 0.0), hd.get("CCD3TEMP", 0.0), \
                                hd.get("CCD4TEMP", 0.0), hd.get("CCD5TEMP", 0.0)

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

                        else:
                            # ultracam / ultraspec

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

                        # Some stuff common to all

                        # Telescope name
                        brow.append(telescope)

                        try:
                            # file size
                            if instrument == 'HiPERCAM':
                                mbytes = round(os.stat(runname + '.fits').st_size / (1024*1024),1)
                            else:
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

                        if len(brow) != len(COLNAMES):
                            print(
                                f'{runname}: data items vs colnames' +
                                f' = {len(brow)} vs {len(COLNAMES)}'
                            )

                        barr.append(brow)

                    # finish off the night file
                    nhtml.write(NIGHT_FOOTER.format(links=links))

                if okwrite:
                    # rename night file
                    shutil.move(fname_tmp, fname)

            ihtml.write('</td></tr>\n')

        if 'Others' not in rnames:
            ihtml.write('</table>\n')

        # finish the main index
        ihtml.write(INDEX_FOOTER)

    if okwrite:
        # rename index file (means that any old index
        # file still exists while new one is being written)
        shutil.move(index_tmp, index)

        # write out the css file
        css = os.path.join(root, f"hiper.css")
        with open(css, "w") as fout:
            fout.write(LOG_CSS)

        # write out the js file
        css = os.path.join(root, f"hiper.js")
        with open(css, "w") as fout:
            fout.write(LOGS_JS)

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
The sqlite3 database file <a href="{linstrument}.db">{linstrument}.db</a> is
designed to allow SQL queries to be applied using Python's sqlite3 module.
The database contains a single table called "{linstrument}". The columns of
the database are defined below. Note although some are naturally integers,
they are converted to floats to allow null values due to details of pandas.
The pipeline command "logsearch" searches the database files and is a good
place to look to see how it is done.

<p>
Here is a simple example of carrying out a search for ULTRACAM runs. The position
selected matches AR Sco and in this case only runs longer than 200 seconds are
returned (37 runs in the ULTRACAM database as of Feb 2021).
<pre>
import sqlite3
import pandas as pd

# Connect to the database
cnx = sqlite3.connect('ultracam.db')

# Query string
query = """
SELECT * FROM ultracam
WHERE ra_deg > 245 AND ra_deg < 246
AND dec_deg > -23 AND dec_deg < -22 and total > 200"""

# return the result as a pandas DataFrame
res = pd.read_sql_query(query, cnx)
print(res)
cnx.close()
</pre>

<p>
Here is a more complex example in which we carry out the same simple
search as before to derive a sub-table labelled t2, but then search
through the original table (t1) for all runs taken within 2 days of
the t2 ones which match them in terms of format. This returns matching
calibration frames, checks that the filters are right for flats and
the red speed is OK for biases. It also writes the results to disk.

<pre>
import sqlite3
import pandas as pd
from hipercam.utils import format_hlogger_table

# Connect to the database
cnx = sqlite3.connect('ultracam.db')

query = """
SELECT t1.*
FROM ultracam AS t1, (
   SELECT *
   FROM ultracam
   WHERE ra_deg > 245 AND ra_deg < 246 AND dec_deg > -23 AND dec_deg < -22
   AND total > 200
) as t2
WHERE ABS(t1.mjd_start-t2.mjd_start) < 2
AND (t1.run_type != 'flat' OR t1.filters == t2.filters)
AND (t1.run_type != 'bias' OR t1.read_speed == t2.read_speed)
AND (
     (t1.win1 == t2.win1 AND t1.binning == t2.binning
     AND ((t1.ra_hms == t2.ra_hms) OR t1.ra_deg is NULL))
     OR (t1.win1 == '1,513,1,512,1024' AND t1.binning == '1x1')
         AND ((t1.ra_hms == t2.ra_hms) OR t1.ra_deg is NULL)
    )
GROUP BY t1.night, t1.run_no
"""

res = pd.read_sql_query(query, cnx)
cnx.close()
print(res)
format_hlogger_table('results.xlsx', res, 'ultracam')
</pre>

<p>
This search comes up with 283 runs as of 22 Feb; there tend to be a fair
few matching calibrations for each bit of real data, but this does cut
down significantly on what to look through.

<h2>Column names, definitions and data types</h2>

<p>
The database is called {linstrument}.db and contains a single table called
'{linstrument}' with the following columns:

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

        print('\nFinished generation of the web pages.')

        if args.full or args.positions or not args.quick:
            print('Now generating a spreadsheet and an SQL database')

            # create pd.DataFrame containing all info
            ptable = pd.DataFrame(data=barr,columns=cnames)

            # enforce data types
            ptable = ptable.astype(dtypes)

            spreadsheet = os.path.join(root, f"{linstrument}.xlsx")
            format_hlogger_table(spreadsheet, ptable, linstrument)
            print(f'Written spreadsheet to {linstrument}.xlsx')

            # write out sqlite database
            sqldb = os.path.join(root, f'{linstrument}.db')
            cnx = sqlite3.connect(sqldb)
            ptable.to_sql(name=f'{linstrument}', con=cnx, if_exists='replace')
            cnx.commit()
            cnx.close()
            print(f'Written sqlite database to {linstrument}.db')
            print(f'Table dimensions (rows,columns) = {ptable.shape}')
        else:
            print("""
Skipping generation of the spreadsheet and SQL database as not
all nights have been processed. Use the "-l" switch to get a
slower update that includes the spreadsheet and database, and
see other switches for even slower and more in-depth options.""")

        print(f'\nAll done. Look in {root} for the outputs.')

    if len(smessages):
        print('\nYou may want to add the following to TARGETS to short-circuit SIMBAD lookups:\n')
        print('\n'.join(smessages))

        # also save to disk
        ofile = f'lookup_success_{iso_time}'
        with open(ofile,'w') as fsim:
            fsim.write('\n'.join(smessages))
        print(f'\nList of successful lookups written to {ofile}')

    if len(fmessages):
        print('\nYou may want to check these or add them to either SKIP_TARGETS or FAILED_TARGETS to short-circuit SIMBAD lookups:\n')
        print('\n'.join(fmessages))

        # also save to disk
        ofile = f'lookup_failure_{iso_time}'
        with open(ofile,'w') as fsim:
            fsim.write('\n'.join(fmessages))
        print(f'\nList of failed lookups written to {ofile}')

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
            rec = re.compile('file\s+object\s+filter', re.I)
            rre = re.compile('(run\d\d\d\d?)(.*)')
            old = re.compile('\s*(\S+)\s+(\S+)\s+(.*)$')
            oldii = re.compile('\s*(\S+)\s*$')
            with open(fname) as f:
                for line in f:
                    m = rec.search(line)
                    if m:
                        self.format = 1
                        if len(self.comment):
                            raise Exception('Error in night log = ' + fname + ', line = ' + line)

                    mr = rre.match(line)
                    if mr:
                        run = mr.group(1)
                        if self.format == 2:
                            self.comment[run] = mr.group(2).strip()
                        else:
                            m = old.search(mr.group(2))
                            if m:
                                self.target[run]  = m.group(1)
                                self.filters[run] = m.group(2)
                                self.comment[run] = m.group(3)
                            else:
                                m = oldii.search(mr.group(2))
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


def load_skip_fail():
    """Looks for a file called SKIP_TARGETS which simply contains a column
    of target names to skip (e.g. "Tungsten flat"), and another called
    FAILED_TARGETS which contains a list of targets that have failed
    target lookup, which will also be skipped for now.  FAILED_TARGETS
    contains other information following the target name, hence any
    spaces in the target name are replaced by '~' so that it can be
    clearly identified.

    Simply returns with a list of target names that will be skipped
    for speed. Returns with an empty list if neither file is found.

    """

    if os.path.exists('SKIP_TARGETS'):
        with open('SKIP_TARGETS') as fp:
            skip_targets = fp.readlines()
            skip_targets = [name.strip() for name in skip_targets if not name.startswith('#')]
        print(f'Loaded {len(skip_targets)} target names to skip from SKIP_TARGETS')
    else:
        print('No file called SKIP_TARGETS')
        skip_targets = []

    if os.path.exists('FAILED_TARGETS'):
        nfail = 0
        with open('FAILED_TARGETS') as fp:
            for line in fp:
                if not line.startswith('#'):
                    arr = line.split()
                    target = arr[0].replace('~',' ')
                    skip_targets.append(target.replace('~',' '))
                    nfail += 1

        print(f'Loaded {nfail} target names from FAILED_TARGETS')
    else:
        print('No file called FAILED_TARGETS')

    return skip_targets


class TimeHMSCustom(TimeISO):
    """
    Just the HMS part of a Time as "<HH>:<MM>:<SS.sss...>".
    """

    name = "hms_custom"  # Unique format name
    subfmts = (("date_hms", "%H%M%S", "{hour:02d}:{min:02d}:{sec:02d}"),)


def make_times(night, runs, observatory, times, full, instrument, okwrite):
    """Generates timing data for a set of runs from a particular night.
    Results are stored in a file called "times", and returned as
    dictionary keyed on run numbers.

    """

    # use this to check times are vaguely right. time of runs
    # must lie between 06.00 local time on date corresponding to
    # start of night date and 1.5 days later. Has picked up a
    # few erroneously dated nights on the TNT.
    mjd_ref = Time(night).mjd - observatory.lon.degree/360 + 0.25

    tdata = {}
    with open(times if okwrite else os.devnull,'w') as tout:
        for run in runs:
            if full:
                print(f'Analysing times for run {run}')
            dfile = os.path.join(night, run)
            try:
                ntotal = 0
                if instrument == 'HiPERCAM':
                    rtime = hcam.hcam.Rtime(dfile)
                else:
                    rtime = hcam.ucam.Rtime(dfile)

                # Find first good time, has to roughly match the start
                # date of the night because some times can just be
                # junk
                not_alerted = True
                for n, tdat in enumerate(rtime):
                    if instrument == 'HiPERCAM':
                        time, tinfo, tflag = tdat
                        expose = 1000000
                        for tmid,texp,tiflag in tinfo:
                            expose = min(round(texp,3),expose)
                    else:
                        time, tinfo = tdat[:2]
                        tflag = time.good
                        expose = round(time.expose,3)

                    if instrument == 'HiPERCAM' or tflag:
                        mjd_start = time.mjd
                        tdelta = mjd_start-mjd_ref
                        if tdelta > 0 and tdelta < 1.5:
                            ts = Time(mjd_start, format="mjd", precision=2)
                            ut_start = ts.hms_custom
                            n_start = n+1
                            if expose >= 0 and expose < 2000:
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
                if instrument == 'HiPERCAM':
                    nback = 4
                elif rtime.header['MODE'] == 'DRIFT':
                    # ultracam or hipercam
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

                if instrument == 'HiPERCAM':
                    ntotal = rtime.ntotal()
                else:
                    nbytes = os.stat(dfile + '.dat').st_size
                    ntotal = nbytes // rtime.framesize

                if instrument != 'HiPERCAM' and ntotal > 20000:
                    # this is a risk-reducing strategy in case the end
                    # of a long ultracam or ultraspec run is
                    # corrupt. Better to look at more than the
                    # necessary number of frames if it prevents us
                    # from having to wind through the whole lot.
                    nback = max(nback, 500)

                # next statement basically resets the frame
                # we are on
                nreset = max(1, ntotal - nback)
                rtime.set(nreset)

                flast = False
                for n, tdat in enumerate(rtime):
                    if instrument == 'HiPERCAM':
                        time, tinfo, tflag = tdat
                        nexpose = 1000000
                        for tmid,texp,tiflag in tinfo:
                            nexpose = min(round(texp,3),expose)
                    else:
                        time, tinfo = tdat[:2]
                        tflag = time.good
                        nexpose = round(time.expose,3)

                    if instrument == 'HiPERCAM' or tflag:
                        mjd = time.mjd
                        if mjd >= mjd_start and mjd < mjd_start + 0.4:
                            mjd_end = mjd
                            ts = Time(mjd_end, format="mjd", precision=2)
                            ut_end = ts.hms_custom
                            n_end = nreset + n
                            if nexpose < 2000:
                                expose = max(expose, nexpose)
                                flast = True

                if not flast:
                    # no good time found near end. There must be
                    # one or we wouldn't get to this point, so
                    # grind it out the hard way by going through
                    # the whole run, which can be slow.
                    rtime.set()
                    for n, tdat in enumerate(rtime):
                        if instrument == 'HiPERCAM':
                            time, tinfo, tflag = tdat
                            nexpose = 1000000
                            for tmid,texp,tiflag in tinfo:
                                nexpose = min(round(texp,3),expose)
                        else:
                            time, tinfo = tdat[:2]
                            tflag = time.good
                            nexpose = round(time.expose,3)

                        if tflag:
                            mjd = time.mjd
                            if mjd >= mjd_start and mjd < mjd_start + 0.4:
                                mjd_end = mjd
                                ts = Time(mjd_end, format="mjd", precision=2)
                                ut_end = ts.hms_custom
                                n_end = n + 1
                                if nexpose < 2000:
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

    if okwrite:
        print('Written timing data to',times)

    return tdata

def make_positions(
        night, runs, observatory, instrument, hlog, targets,
        skip_targets, tdata, posdata, load_old, retry,
        full, rname, smessages, fmessages, p2positions, okwrite
):
    """Determine positional info, write to podata, return as dictionary
    keyed on the runs. Uses pre-determined timing data from
    make_times.

    smessages and fmessages are lists that should be initialised to []
    that are used to accumulate successful and failed target lookup
    messages.

    Arguments::

      night : night name
      runs : the runs to process
      observatory : telescope
      instrument : instrument name
      hlog : hand-written log
      targets : targets with pre-existing position data to avoid SIMBAD
      skip_targets : list of target names not to bother with
      tdata : timing data
      posdata : name of file containing positional data
      load_old : whether to start by trying to read pre-stored data
      retry : whether to attempt to re-do any targets with undefined positions.
              irrelevant if load_old == False. If not, and retry = False, then
              having loaded data, the function returns immediately, else it goes
              through the targets retrying any that have no positional data.
      full : lots of info printed
      rname : run name
      smessages : buffer of successful target lookup messages
      fmessages : buffer of failed target lookup messages
      p2positions : positions keyed by target name (RA, Dec pairs, sexagesimal)
                    used as a last resort

    Returns with dictionary of positional data. Each entry is keyed by run name
    and contains the following 19 items:

    ra dec autoid alt1 alt2 alt3 az1 az2 az3 seczmin seczmax seczdelta
    sund moond salt1 salt2 malt1 malt2 sm.

    Corresponding data files start with the run number so have 20 items per line.

    """

    pdata = {}

    if load_old and os.path.exists(posdata):
        # pre-existing file found
        with open(posdata) as pin:
            for line in pin:
                arr = line.split()
                if len(arr) != 20:
                    raise ValueError(
                        f'Line = "{line.strip()}" from {posdata} had {len(arr)}!=20 items'
                    )
                arr[3] = arr[3].replace('~',' ')
                pdata[arr[0]] = [
                    '' if val == 'UNDEF' else val for val in arr[1:]
                ]
        print('Read position data from',posdata)

        if not retry:
            return pdata

    with open(posdata if okwrite else os.devnull,'w') as pout:
        for run in runs:

            if len(tdata[run]) == 1:
                # means its a power on/off
                continue

            if run in pdata and pdata[run][0] != '':
                # Already have positional data which we will
                # not re-do, so just write out to disk
                arr = ['UNDEF' if val == '' else val for val in pdata[run]]
                arr[2] = arr[2].replace(' ','~')
                pout.write(
                    f"{run} {arr[0]} {arr[1]} {arr[2]} {arr[3]} {arr[4]} " +
                    f"{arr[5]} {arr[6]} {arr[7]} {arr[8]} {arr[9]} {arr[10]} " +
                    f"{arr[11]} {arr[12]} {arr[13]} {arr[14]} {arr[15]} " +
                    f"{arr[16]} {arr[17]} {arr[18]}\n"
                )
                continue

            recomp = True

            # Now going to try to work stuff out

            if full:
                print(f'Analysing positions for run {run}')

            # open the run file as an Rhead
            runname = os.path.join(night, run)
            try:
                if instrument == 'HiPERCAM':
                    rhead = hcam.hcam.Rhead(runname)
                else:
                    rhead = hcam.ucam.Rhead(runname)
            except:
                if full:
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_tb(exc_traceback, limit=1)
                    traceback.print_exc()
                    print(f"Failed to open {runname} as an Rhead")
                continue

            # object name
            if hlog.format == 1:
                target = hlog.target[run]
            elif instrument == 'HiPERCAM':
                target = rhead.header.get("OBJECT",'')
            else:
                target = rhead.header.get("TARGET",'')
            target = target.strip().replace('~',' ')

            # RA, Dec lookup
            if target == '' or target in skip_targets:
                # don't even try
                autoid, ra, dec = 'UNDEF', 'UNDEF', 'UNDEF'
                recomp = False
            else:
                try:
                    # See if we already have the info stored
                    autoid, ra, dec = targets(target)
                except:
                    # apparently we don't ...
                    try:
                        # attempt simbad lookup here
                        autoid, ra, dec = target_lookup(target)
                        targets.add_target(target, ra, dec, autoid)
                        print(f'  Added {target} to targets')
                        pos = SkyCoord(f'{ra} {dec}',unit=(u.hourangle, u.deg))

                        # save successful SIMBAD-based lookup
                        smessages.append(
                            f"{autoid.replace(' ','~'):32s} " +
                            f"{pos.to_string('hmsdms',sep=':',precision=2)} " +
                            f"{target.replace(' ','~')}"
                        )

                    except:
                        if target in p2positions:
                            # data loaded at the phase II stage -- last resort
                            ra, dec = p2positions[target]
                            print(f'  Found {target} in phaseII data at RA={ra}, Dec={dec}')
                            pos = SkyCoord(f'{ra} {dec}',unit=(u.hourangle, u.deg))
                            targets.add_target(target, pos.ra.hour, pos.dec.value, target)
                            autoid, ra, dec = targets(target)

                            # save successful lookups
                            smessages.append(
                                f"{target.replace(' ','~'):32s} " +
                                f"{pos.to_string('hmsdms',sep=':',precision=2)} " +
                                f"{target.replace(' ','~')}"
                            )

                        else:
                            # nothing worked
                            print(
                                f'  No position found for {runname}, target = "{target}"'
                            )
                            autoid, ra, dec = 'UNDEF', 'UNDEF', 'UNDEF'
                            skip_targets.append(target)

                            # save in suitable format for adding to FAILED_TARGETS if wanted.
                            fmessages.append(
                                f"{target.replace(' ','~'):32s} {rname} {night} {run}"
                            )
                            recomp = False

            if not recomp and run in pdata:
                # can save a stack of time by not recomputing any Sun / Moon stuff
                arr = ['UNDEF' if val == '' else val for val in pdata[run]]
                arr[2] = arr[2].replace(' ','~')
                pout.write(
                    f"{run} {arr[0]} {arr[1]} {arr[2]} {arr[3]} {arr[4]} " +
                    f"{arr[5]} {arr[6]} {arr[7]} {arr[8]} {arr[9]} {arr[10]} " +
                    f"{arr[11]} {arr[12]} {arr[13]} {arr[14]} {arr[15]} " +
                    f"{arr[16]} {arr[17]} {arr[18]}\n"
                )
                continue

            # start accumulating stuff to write out
            arr = [ra, dec, autoid]

            if ra == 'UNDEF' and dec == 'UNDEF' and instrument == 'ULTRASPEC':
                # for altitude / Sun / Moon stuff, telescope position
                # is good enough, so this is one final go at getting a
                # usable position.
                hd = rhead.header

                ra = hd.get("RA", "UNDEF")
                dec = hd.get("Dec", "UNDEF")
                if ra != 'UNDEF' and dec != 'UNDEF':
                    try:
                        ra, dec, syst = str2radec(ra + ' ' + dec)
                    except:
                        pass

            # time-dependent info
            ut_start, mjd_start, ut_end, mjd_end, cadence, \
                expose, nok, ntotal = tdata[run]

            try:

                mjd_start = float(mjd_start)
                mjd_end = float(mjd_end)
                tstart = Time(mjd_start, format='mjd')
                tmid = Time((mjd_start+mjd_end)/2, format='mjd')
                tend = Time(mjd_end, format='mjd')

                # Scale Sun-Moon angle at mid time (0 = New Moon, 1 =
                # Full)
                sun_mid = get_sun(tmid)
                moon_mid = get_moon(tmid)
                sun_moon = sun_mid.separation(moon_mid).degree / 180

                if ra != 'UNDEF' and dec != 'UNDEF':

                    # Calculate the Alt, Az at start, middle, end
                    frames = AltAz(obstime=[tstart,tmid,tend], location=observatory)
                    pos = SkyCoord(f'{ra} {dec}',unit=(u.hourangle, u.deg))
                    points = pos.transform_to(frames)
                    alts = [round(alt,1) for alt in points.alt.degree]
                    azs = [round(az,1) for az in points.az.degree]
                    arr += alts + azs

                    # Calculate range of airmasses
                    seczs = np.array([float(secz) for secz in points.secz])
                    secz_min, secz_max = seczs.min(), seczs.max()

                    # Need to check for meridian crossing, and if it happens
                    # we need to close in on it
                    sinas = [np.sin(az) for az in points.az]
                    if sinas[0] > 0 and sinas[2] < 0:
                        s1, s2 = sinas[0], sinas[2]
                        t1, t2 = tstart, tend
                        if sinas[1] > 0:
                            s1 = sinas[1]
                            t1 = tmid
                        else:
                            s2 = sinas[1]
                            t2 = tmid
                        while s1 - s2 > 0.0005:
                            tguess = t1 + s1/(s1-s2)*(t2-t1)
                            frame = AltAz(obstime=tguess, location=observatory)
                            point = pos.transform_to(frame)
                            sina = np.sin(point.az)
                            if sina > 0:
                                s1 = sina
                                t1 = tguess
                            else:
                                s2 = sina
                                t2 = tguess
                        secz_min = float(point.secz)

                    dsecz = round(secz_max-secz_min,2)
                    arr += [round(secz_min,2), round(secz_max,2), dsecz]

                    # Now calculate the angular distance from the Sun
                    # and Moon at the mid-time
                    sun_mid_trans = sun_mid.transform_to(frames[1])
                    moon_mid_trans = moon_mid.transform_to(frames[1])
                    point_mid = points[1]
                    sun_dist = point_mid.separation(sun_mid_trans).degree
                    moon_dist = point_mid.separation(moon_mid_trans).degree
                    arr += [round(sun_dist,1),round(moon_dist,1)]

                else:
                    arr = arr[:3] + 11*['UNDEF']

                # Now some data on the altitude of the Sun & Moon
                frame = AltAz(obstime=tstart, location=observatory)
                sun_start = get_sun(tstart).transform_to(frame)
                moon_start = get_moon(tstart).transform_to(frame)

                # end
                frame = AltAz(obstime=tend, location=observatory)
                sun_end = get_sun(tend).transform_to(frame)
                moon_end = get_moon(tend).transform_to(frame)

                arr += [
                    round(sun_start.alt.degree,1), round(sun_end.alt.degree,1),
                    round(moon_start.alt.degree,1), round(moon_end.alt.degree,1),
                    round(sun_moon,3),
                ]

            except:
                if full:
                    print(f"Problem on run = {run}")
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_tb(
                        exc_traceback, limit=1, file=sys.stdout
                    )
                    traceback.print_exc(file=sys.stdout)

                # write out info
                arr = arr[:3] + 16*['UNDEF']

            arr[2] = arr[2].replace(' ','~')
            pout.write(
                f"{run} {arr[0]} {arr[1]} {arr[2]} {arr[3]} {arr[4]} " +
                f"{arr[5]} {arr[6]} {arr[7]} {arr[8]} {arr[9]} {arr[10]} " +
                f"{arr[11]} {arr[12]} {arr[13]} {arr[14]} {arr[15]} " +
                f"{arr[16]} {arr[17]} {arr[18]}\n"
            )

            arr[2] = arr[2].replace('~',' ')
            pdata[run] = [
                '' if val == 'UNDEF' else val for val in arr
            ]

    if okwrite:
        print('Written positional data to',posdata)

    return pdata

HIPERCAM_COLNAMES = (
    ('obs_run', 'str', 'Run group title'),
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
    ('date_start', 'str', 'Date at the start of the run'),
    ('utc_start', 'str', 'UTC at the start of the run'),
    ('utc_end', 'str', 'UTC at the end of the run'),
    ('mjd_start', 'float64', 'MJD (UTC) at the start of the run'),
    ('mjd_end', 'float64', 'MJD (UTC) at the end of the run'),
    ('total', 'float32', 'Total exposure time (seconds)'),
    ('cadence', 'float32', 'Sampling time between exposures (seconds)'),
    ('exposure', 'float32', 'Actual exposure time (seconds)'),
    ('nframe', 'float32', 'Total number of frames in file'),
    ('nok', 'float32', 'Total number of frames in file'),
    ('filters', 'str', 'Colour filters used'),
    ('run_type', 'str', 'Provisional type of the run, indicative only'),
    ('read_mode', 'str', 'ULTRACAM readout mode'),
    ('nskips', 'str', 'nskips for each CCD'),
    ('quad1', 'str', 'Format of quad 1'),
    ('quad2', 'str', 'Format of quad 2'),
    ('binning', 'str', 'X by Y binning'),
    ('clr', 'str', 'Clear enabled or not'),
    ('dummy', 'str', 'Dummy output being used or not'),
    ('led', 'str', 'LED on or not'),
    ('oscan', 'str', 'Over-scan on or not'),
    ('pscan', 'str', 'Pre-scan on or not'),
    ('nodding', 'str', 'Telescope nodding or not'),
    ('read_speed', 'str', 'Readout speed'),
    ('fclocks', 'str', 'Fast clocks or not'),
    ('tbytes', 'str', 'Timing bytes problem or not'),
    ('fpslide', 'float32', 'Focal plane slide'),
    ('ccdtemps', 'str', 'CCD temps'),
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
    ('secz_min', 'float32', 'Minimum airmass'),
    ('secz_max', 'float32', 'Maximum airmass'),
    ('secz_del', 'float32', 'Change in airmass'),
    ('sun_dist', 'float32', 'Distance from Sun, degrees, middle of run'),
    ('moon_dist', 'float32', 'Distance from Moon, degrees, middle of run'),
    ('sun_alt_start', 'float32', 'Altitude of Sun at start of run'),
    ('sun_alt_end', 'float32', 'Altitude of Sun at end of run'),
    ('moon_alt_start', 'float32', 'Altitude of Moon at start of run'),
    ('moon_alt_end', 'float32', 'Altitude of Moon at end of run'),
    ('moon_phase', 'float32', 'Angle between Sun and Moon / 180')
)

ULTRACAM_COLNAMES = (
    ('obs_run', 'str', 'Run group title'),
    ('night' , 'str', 'Date at the start of the night'),
    ('run_no', 'str', 'Run number'),
    ('target', 'str', 'Target name'),
    ('auto_id', 'str', 'Lookup name, from disk files or SIMBAD'),
    ('ra_hms', 'str', 'Target RA (J2000), HMS'),
    ('dec_dms', 'str', 'Target Dec (J2000), DMS'),
    ('ra_deg', 'float64', 'Target RA (J2000), degrees'),
    ('dec_deg', 'float64', 'Target Dec (J2000), degrees'),
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
    ('secz_min', 'float32', 'Minimum airmass'),
    ('secz_max', 'float32', 'Maximum airmass'),
    ('secz_del', 'float32', 'Change in airmass'),
    ('sun_dist', 'float32', 'Distance from Sun, degrees, middle of run'),
    ('moon_dist', 'float32', 'Distance from Moon, degrees, middle of run'),
    ('sun_alt_start', 'float32', 'Altitude of Sun at start of run'),
    ('sun_alt_end', 'float32', 'Altitude of Sun at end of run'),
    ('moon_alt_start', 'float32', 'Altitude of Moon at start of run'),
    ('moon_alt_end', 'float32', 'Altitude of Moon at end of run'),
    ('moon_phase', 'float32', 'Angle between Sun and Moon / 180')
)

ULTRASPEC_COLNAMES = (
    ('obs_run', 'str', 'Run group title'),
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
    ('secz_min', 'float32', 'Minimum airmass'),
    ('secz_max', 'float32', 'Maximum airmass'),
    ('secz_del', 'float32', 'Change in airmass'),
    ('sun_dist', 'float32', 'Distance from Sun, degrees, middle of run'),
    ('moon_dist', 'float32', 'Distance from Moon, degrees, middle of run'),
    ('sun_alt_start', 'float32', 'Altitude of Sun at start of run'),
    ('sun_alt_end', 'float32', 'Altitude of Sun at end of run'),
    ('moon_alt_start', 'float32', 'Altitude of Moon at start of run'),
    ('moon_alt_end', 'float32', 'Altitude of Moon at end of run'),
    ('moon_phase', 'float32', 'Angle between Sun and Moon / 180')
)


def noval(value):
    return None if value == '' else value

# Header of the main file

INDEX_HEADER = """<!DOCTYPE html>
<html>
<head>
<title>{instrument} data logs</title>
</head>

<body>
<h1>{instrument} data logs</h1>

<p> These on-line logs summarise all {instrument} runs. They were
automatically generated from the data files and hand-written logs.

<p>For a searchable file summarising the same information, plus
calculated data on the Sun and Moon along with a few summary statistics
generated by "hmeta" for each run, see this <a
href="{linstrument}.xlsx">spreadsheet</a>. This is also
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

NIGHT_HEADER = """<!DOCTYPE html>
<html>
    <head>
       <link rel="stylesheet" type="text/css" href="../hiper.css" />
       <title>{instrument} log {date}</title>
       <script src="../hiper.js"></script>
    </head>

  <body>

    <h1>{instrument} log {date}</h1>

    <p>
      The table below lists information on runs from the night starting
      on {date}.  See the end of the table details on the meanings of the
      various columns. 
    </p>

{links}

<p>
<button id="More" onClick="showDetails(true)">More detail</button>
<button id="Less" onClick="showDetails(false)">Less detail</button>
</p>

<p>
<div class="tableFixHead">
<table>
"""

HIPERCAM_TABLE_HEADER = """
<col span="2">
<col span="6" class="hide">
<col>
<col class="hide">
<col span="4">
<col class="hide">
<col>
<col class="hide">
<col span="2">
<col span="3" class="hide">
<col span="2">
<col span="5" class="hide">
<col>
<col span="9" class="hide">
<col>

<thead>
<tr>
<th class="left">Run<br>no.</th>
<th class="left">Target name</th>
<th class="left">Auto ID</th>
<th class="left">RA (J2000)</th>
<th class="left">Dec&nbsp;(J2000)</th>
<th class="left">Tel RA</th>
<th class="left">Tel Dec</th>
<th class="left">PA</th>
<th class="cen">Start<br>UTC</th>
<th class="cen">End<br>UTC</th>
<th class="cen">Total<br>(sec)</th>
<th class="right">Cad.<br>(sec)</th>
<th class="right">Exp.<br>(sec)</th>
<th class="right">Nfrm</th>
<th class="right">Nok</th>
<th class="right">Delta<br>secz</th>
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
<th class="cen">Read<br>speed</th>
<th class="cen">Fast<br>clcks</th>
<th class="cen">Tbytes</th>
<th class="cen">FPslide</th>
<th class="cen">CCD temps</th>
<th class="cen">Observers</th>
<th class="left">PI</th>
<th class="left">PID</th>
<th class="cen">Size<br>(MB)</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
</thead>

<tbody>
"""

ULTRACAM_TABLE_HEADER = """
<col span="2">
<col span="3" class="hide">
<col>
<col class="hide">
<col span="4">
<col class="hide">
<col span="11">
<col span="6" class="hide">
<col>

<thead>
<tr>
<th class="left">Run<br>no.</th>
<th class="left">Target name</th>
<th class="left">Auto ID</th>
<th class="left">RA (J2000)</th>
<th class="left">Dec&nbsp;(J2000)</th>
<th class="cen">Start<br>UTC</th>
<th class="cen">End<br>UTC</th>
<th class="right">Total<br>(sec)</th>
<th class="right">Cad.<br>(sec)</th>
<th class="right">Exp.<br>(sec)</th>
<th class="right">Nfrm</th>
<th class="right">Nok</th>
<th class="right">Delta<br>secz</th>
<th class="cen">Filters</th>
<th class="left">Type</th>
<th class="cen">Read<br>mode</th>
<th class="left">Nb</th>
<th class="cen">xl,xr,ys,nx,ny<br>window pair 1</th>
<th class="cen">xl,xr,ys,nx,ny<br>window pair 2</th>
<th class="cen">xl,xr,ys,nx,ny<br>window pair 3</th>
<th class="cen">XxY<br>bin</th>
<th class="cen">Clr</th>
<th class="cen">Read<br>speed</th>
<th class="cen">FPslide</th>
<th class="cen">Observers</th>
<th class="left">PI</th>
<th class="left">PID</th>
<th class="cen">Size<br>(MB)</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
</thead>

<tbody>
"""

ULTRASPEC_TABLE_HEADER = """
<col span="2">
<col span="6" class="hide">
<col>
<col class="hide">
<col span="4">
<col class="hide">
<col span="11">
<col span="8" class="hide">
<col>

<thead>
<tr>
<th class="left">Run<br>no.</th>
<th class="left">Target name</th>
<th class="left">Auto ID</th>
<th class="left">RA (J2000)</th>
<th class="left">Dec&nbsp;(J2000)</th>
<th class="left">Tel RA</th>
<th class="left">Tel Dec</th>
<th class="left">PA</th>
<th class="cen">Start<br>UTC</th>
<th class="cen">End<br>UTC</th>
<th class="right">Total<br>(sec)</th>
<th class="right">Cad.<br>(sec)</th>
<th class="right">Exp.<br>(sec)</th>
<th class="right">Nfrm</th>
<th class="right">Nok</th>
<th class="right">Delta<br>secz</th>
<th class="cen">Filters</th>
<th class="left">Type</th>
<th class="cen">Read<br>mode</th>
<th class="cen">xs,ys,nx,ny<br>window 1</th>
<th class="cen">xs,ys,nx,ny<br>window 2</th>
<th class="cen">xs,ys,nx,ny<br>window 3</th>
<th class="cen">xs,ys,nx,ny<br>window 4</th>
<th class="cen">XxY<br>bin</th>
<th class="cen">Clr</th>
<th class="cen">Read<br>speed</th>
<th class="cen">Output</th>
<th class="cen">HV Gain</th>
<th class="cen">FPslide</th>
<th class="cen">Observers</th>
<th class="left">PI</th>
<th class="left">PID</th>
<th class="cen">Size<br>(MB)</th>
<th class="left">Run<br>no.</th>
<th class="left">Comment</th>
</tr>
</thead>

<tbody>
"""

NIGHT_FOOTER = """
</tbody>
</table>
</div>

{links}

<p> The night date is the date at the start of the night; 'Total' is
the run length; 'Cad.' is the cadence or sampling time; 'Exp.'  is the
exposure time per point; 'PA' is the PA on the sky (ULTRASPEC); 'Clr'
indicates whether clears were enabled; 'Read mode' is the readout mode
which can be one of several options: 'FFCLR' for full frames with
clear; 'FFNCLR' full frames with no clear; '1-PAIR', '2-PAIR',
'3-PAIR', for standard windowed modes, etc. The 'xl,xr,..' columns
gives the parameters defining the windows (5-per-pair for ULTRACAM;
4-per-window for ULTRASPEC); 'Nfrm' is the number of frame; 'Nok' is
the number of OK frames judged as having OK times. If it grossly
disagrees with 'Nframe', there might be timing problems present, such
as the GPS dropping out; 'Nb' is the nblue paramter of ULTRACAM.  </p>

<address>Tom Marsh, Warwick</address>
</body>
</html>
"""

# Javascript controlling hiding / uncovering of parts of the night logs
# for clarity / detail. Various elements have ids which are picked out by
# this little script. The brief form is selected by default at the start
LOGS_JS = """

function showDetails(showDetail) {

  // shows more or less detail according to the value of "showDetail"
  let els = document.getElementsByClassName("hide");
  for(let i = 0; i < els.length; i++) {
     els[i].style.visibility = showDetail ? "visible" : "collapse";
  }

  // options to hide / restore stuff
  let x = document.getElementById("More");
  x.style.display = showDetail ? "none" : "inline";

  x = document.getElementById("Less");
  x.style.display = showDetail ? "inline" : "none";

};

function initialise() {
  // turn off stuff at start
  showDetails(false);
};

document.addEventListener('DOMContentLoaded', initialise);

"""

def tradec(radec, isdec):
    """
    Correct for improper formatted RA, Dec
    """

    try:
        v1,v2,v3 = radec.split(':')
        if v1.strip().startswith('-'):
            # careful because of '-00' possibility
            v1 = -abs(int(v1))
        else:
            v1 = int(v1)
        v2,v3 = int(v2),float(v3)
        if isdec:
            return f'{v1:+03d}:{v2:02d}:{v3:04.1f}'
        else:
            return f'{v1:02d}:{v2:02d}:{v3:05.2f}'
    except:
        return ''

# For hipercam
TRANSLATE_MODE = {
    "FullFrame": "FULL",
    "OneWindow": "1-QUAD",
    "TwoWindows": "2-QUAD",
    "DriftWindow": "DRIFT",
}


# Run name correctors for picking up positions
CORRECTORS = {
    '2021-P107' : 'P107',
    '2021-P108' : 'P108',
    '2022-P109' : 'P109'
}
