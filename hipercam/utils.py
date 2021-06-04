"""
Classes and functions of general use
"""

import os
import sys
import math
import re
import requests
import hipercam as hcam

import numpy as np
import pandas as pd
from .core import *

__all__ = (
    "Vec2D", "add_extension", "sub_extension", "script_args", "rgb",
    "format_hlogger_table", "what_flags", "temp_dir"
)


class Vec2D:
    """Simple 2D vector class."""

    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y

    def length(self):
        """Returns the Euclidean length"""
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def unit(self):
        """Returns a unit vector version of the vector"""
        dist = self.length()
        if dist > 0:
            x = self.x / dist
            y = self.y / dist
            return Vec2D(x, y)
        else:
            raise ValueError("cannot normalise a zero vector")

    def __iadd__(self, vec):
        """+=: in place addition"""
        self.x += vec.x
        self.y += vec.y
        return self

    def __add__(self, vec):
        """+: addition"""
        return Vec2D(self.x + vec.x, self.y + vec.y)

    def __isub__(self, vec):
        """-=: in place subtraction"""
        self.x -= vec.x
        self.y -= vec.y
        return self

    def __sub__(self, vec):
        """-: subtraction"""
        return Vec2D(self.x - vec.x, self.y - vec.y)

    def __imul__(self, const):
        """*=: in place multiplication by a constant"""
        self.x *= const
        self.y *= const
        return self

    def __mul__(self, const):
        """*: post multiplication by a constant"""
        return Vec2D(const * self.x, const * self.y)

    def __rmul__(self, const):
        """*: pre multiplication by a constant"""
        return Vec2D(const * self.x, const * self.y)

    def __itruediv__(self, const):
        """Divides by a constant"""
        self.x /= const
        self.y /= const

    def __repr__(self):
        return "Vec2D(x={:f}, y={:f})".format(self.x, self.y)


def dot(v1, v2):
    """Returns the scalar or 'dot' product of two vectors"""
    return v1.x * v2.x + v1.y * v2.y


def add_extension(fname, ext):
    """Add extension ext to a file name if it is not already there, and returns
    the revised name

    """
    if len(ext) and not fname.endswith(ext):
        return "{}{}".format(fname, ext)
    else:
        return fname


def sub_extension(fname, ext):
    """Subtracts extension ext from a file name if it is present, and returns
    the revised name

    """
    if fname.endswith(ext):
        return fname[: -len(ext)]
    else:
        return fname


def rgb(cname):
    """Returns the RGB tuple associated with colour name 'cname', following
    the colour definitions used by 'reduce'.

    Returns a ValueError if cname is not recognised.
    """
    if cname not in CNAMS:
        raise ValueError("colour = {:s} not recognised".format(cname))
    else:
        return CIS[CNAMS[cname]]


def script_args(args):
    """
    This is a small helper method that is used at the start of most
    of the HiPERCAM entry point scripts such as 'reduce' and 'grab'.

    If 'args' is None on entry, it is replaced by 'sys.argv'. It is then
    assumed that the first argument is the command name or else None. The
    command name argument is removed from args, and if not None, converted to just the
    file part of the path. It is then returned along with the remaining list
    of arguments in a two element tuple, i.e. (commands, args).
    """
    if args is None:
        args = sys.argv.copy()

    command = args.pop(0)
    if command is not None:
        # just keep filename part, not any path
        command = os.path.split(command)[1]

    return (command, args)


def print_stats(ccd, cnam, x, y, hsbox, warn=True):
    """
    Prints a few statistics of pixels around an x,y position
    in a CCD. It returns a tuple containing a label and a reference
    to the Window containing the x, y position, or None, None
    if the x,y position is outside any window. This is used by 'hplot'
    and 'setdefect'.

    Arguments::

       ccd    : :class:`hipercam.CCD`
          the CCD of interest

       cnam   : string
          the name of the CCD

       x, y   : float, float
          the X,Y position. Lower-left unbinned pixel is centred on 1,1

       hsbox  : int
          half-width in binned pixels of the box

       warn   : bool
          prints out a message if no enclosing Window found. If you want to
          supply your own message, you might want to set this to False.

    Returns:: (wnam, wind)

       wnam  : string
          Window label, None if not found

       wind  : :class:`hipercam.Window`
          Window enclosing x,y, none if not found.
    """

    wnam = ccd.inside(x, y, 0)
    if wnam is not None:
        wind = ccd[wnam]
        ix = int(round(wind.x_pixel(x)))
        iy = int(round(wind.y_pixel(y)))
        ix1 = max(0, ix - hsbox)
        ix2 = min(wind.nx, ix + hsbox + 1)
        iy1 = max(0, iy - hsbox)
        iy2 = min(wind.ny, iy + hsbox + 1)

        print(
            "\nClicked on x,y = {:.2f},{:.2f} in CCD {:s}, window {:s}".format(
                x, y, cnam, wnam
            )
        )

        print(
            " Stats box in window pixels, X,Y = [{:d}:{:d},{:d}:{:d}] ({:d}x{:d}), central pixel = [{:d},{:d}], value = {:.2f}".format(
                ix1, ix2, iy1, iy2, ix2 - ix1, iy2 - iy1, ix, iy, wind.data[iy, ix]
            )
        )

        box = wind.data[iy1:iy2, ix1:ix2]
        print(
            " Mean = {:.4g}, RMS = {:.4g}, min = {:.4g}, max = {:.4g}, median = {:.4g}".format(
                box.mean(), box.std(), box.min(), box.max(), np.median(box)
            )
        )

    else:
        wind = None
        if warn:
            print(
                "\n *** selected position ({:.1f},{:.1f}) not in any window".format(
                    x, y
                )
            )

    return (wnam, wind)


def old_format_hlogger_table(fname, table):
    """

    Formats the column widths etc of the excel file produced by
    |hlogger| ('hipercam-log.xlsx') and writes it to disk. It is made
    available as a utility function to allow the same styling to be
    applied to derived versions of this file.

    Arguments::

      fname : str
         The file name to write to

      table : pandas.DataFrame
         DataFrame containing the data to be written that is assumed to be
         derived from the |hlogger| file. In particular it is expected to
         have columns called 'Nlink', 'Run no.' and 'Night'. The latter two
         must lie somewhere in the column range 1-26 [A-Z].
    """

    writer = pd.ExcelWriter(fname, engine='xlsxwriter')
    table.to_excel(writer, sheet_name='Hipercam logs', index=False)

    # format
    worksheet = writer.sheets["Hipercam logs"]

    # loop through all columns, work out maximum length needed to span
    # columns
    cnames = table.columns

    # integer index of 'Nlink' column
    idx_nlink = cnames.get_loc('Nlink')

    # Labels of 'Night' and 'Run no.' columns
    lab_night = chr(ord('A')+cnames.get_loc('Night'))
    lab_runno = chr(ord('A')+cnames.get_loc('Run no.'))

    def clen(mlen):
        return min(150,int(math.ceil(0.88*mlen+1.5)))

    for idx, col in enumerate(table):
        series = table[col]
        max_len = max(
            series.astype(str).map(len).max(),
            max(map(len, str(series.name).split('\n')))
        )
        if idx == idx_nlink:
            width_nlink = clen(max_len)

        # set column width
        worksheet.set_column(idx, idx, clen(max_len))

    # Fancy stuff with links column
    cell_format = writer.book.add_format({'font_color': 'blue'})

    worksheet.set_column(idx_nlink, idx_nlink, width_nlink, cell_format)

    # set links using a formula
    for nr in range(1,len(table)+1):
        worksheet.write_formula(
            nr, idx_nlink,
            f'=HYPERLINK("http://deneb.astro.warwick.ac.uk/phsaap/hipercam/logs/" & {lab_night}{nr+1}, {lab_runno}{nr+1})'
        )

    # make first row high enough, freeze it, set a decent zoom
    worksheet.set_row(0,28)
    worksheet.freeze_panes(1, 0)
    worksheet.set_zoom(150)
    writer.save()

def format_hlogger_table(fname, table, instrument):
    """Formats the column widths etc of the excel files produced by
    |hlogger| ('ultracam.xlsx', 'ultraspec.xlsx' or 'hipercam.xlsx') and writes
    it to disk. It is made available as a utility function to allow
    the same styling to be applied to derived versions of this file.

    Arguments::

      fname : str
         The file name to write to

      table : pandas.DataFrame
         DataFrame containing the data to be written that is assumed to be
         derived from the |hlogger| file. In particular it is expected to
         have columns called 'nlink', 'run_no' and 'night'. The latter two
         must lie somewhere in the column range 1-26 [A-Z].

      instrument : str
         Either ultracam or ultraspec

    """

    writer = pd.ExcelWriter(fname, engine='xlsxwriter')
    table.to_excel(writer, sheet_name=f'{instrument} logs', index=False)

    # format
    worksheet = writer.sheets[f'{instrument} logs']

    # loop through all columns, work out maximum length needed to span
    # columns
    cnames = table.columns

    # integer index of 'Nlink' column
    idx_nlink = cnames.get_loc('nlink')

    # Labels of 'Night' and 'Run no.' columns
    lab_night = chr(ord('A')+cnames.get_loc('night'))
    lab_runno = chr(ord('A')+cnames.get_loc('run_no'))

    def clen(mlen):
        if mlen is None or np.isnan(mlen):
            return 5
        else:
            return min(150,int(math.ceil(0.88*mlen+1.5)))

    for idx, col in enumerate(table):
        series = table[col]
        max_len = max(
            series.astype(str).map(len).max(),
            max(map(len, str(series.name).split('\n')))
        )
        if idx == idx_nlink:
            width_nlink = clen(max_len)

        # set column width
        worksheet.set_column(idx, idx, clen(max_len))

    # Fancy stuff with links column
    cell_format = writer.book.add_format({'font_color': 'blue'})

    worksheet.set_column(idx_nlink, idx_nlink, width_nlink, cell_format)

    # set links using a formula
    for nr in range(1,len(table)+1):
        worksheet.write_formula(
            nr, idx_nlink,
            f'=HYPERLINK("https://cygnus.astro.warwick.ac.uk/phsaap/{instrument}/logs/" & {lab_night}{nr+1}, {lab_runno}{nr+1})'
        )

    # make first row high enough, freeze it, set a decent zoom
    worksheet.set_row(0,28)
    worksheet.freeze_panes(1, 0)
    worksheet.set_zoom(150)
    writer.save()

def what_flags(bitmask):
    """
    Given a bitmask value, this routine prints which flags have been set, ignoring
    the two catch alls, ANY_FLAG and ALL_GOOD.
    """
    flags_raised = {}
    for fname, flag in FLAGS:
        if flag != ANY_FLAG and flag != ALL_GOOD: 
            if bitmask & flag:
                flags_raised.append(fname)
                flags_raised = ', '.join(flags_raised)

    if len(flags_raised):
        flags_raised = ", ".join(flags_raised)
        print(f"{bitmask} corresponds to flags = {flags_raised}")
    else:
        print(f"{bitmask} means no flags have been set")

def target_lookup(target, ra_tel=None, dec_tel=None, dist=5.5, url='http://simbad.u-strasbg.fr'):
    """
    Tries to determine coordinates of a target in
    one of three ways:

      1) First it attempts to find coordinates in simbad by running
         a query ID with target (unless target is None or ''). 

      2) If (1) fails, try to find coordinates within the name in the
         form JHHMMSS.S[+/-]DDMMSS [can be more precise, but not less]

      3) If 1 and 2 fail or target == '' or None, then attempt a coordinate
         search in simbad within a limit of dist arcmin (default equals half
         diagonal of ULTRASPEC. It needs a unique match in this case. The
         search position comes from ra_tel and dec_tel

    If successful, it comes back with (name, ra, dec), the matched name and position.
    ra, dec are in decimal form (hours, degrees). If it fails, it returns ('UNDEF','UNDEF','UNDEF').
    This is for use by the logging scripts.
    """

    err_msgs = []

    if url.endswith('/'):
        url += 'simbad/sim-script'
    else:
        url += '/simbad/sim-script'
    RESPC = re.compile('\s+')

    if target is not None and target != '':
        # simbad query script
        query = f"""set limit 2
format object form1 "Target: %IDLIST(1) | %COO(A D;ICRS)"
query id {target}"""

        payload = {'submit' : 'submit script', 'script' : query}
        rep = requests.post(url, data=payload)

        data = False
        found = False
        fail = False
        nres = 0
        for line in rep.text.split('\n'):
            if line.startswith('::data::'):
                data = True

            if line.startswith('::error::'):
                err_msgs.append(
                    f'  Error occurred querying simbad for {target}'
                )
                fail = True
                break

            if data:

                if line.startswith('Target:'):
                    nres += 1
                    if nres > 1:
                        err_msgs.append(
                            f'More than one target returned when querying simbad for {target}'
                        )
                        fail = True
                        break

                    name, coords = line[7:].split(' | ')
                    found = True

        if found and not fail:
            # OK we have found one, but we are still not done --
            # some SIMBAD lookups are no good
            name = name.strip()
            name = re.sub(RESPC, ' ', name)
            try:
                ra, dec, syst = str2radec(coords)
                return (name,ra,dec)
            except:
                err_msgs.append(
                    f'Matched "{target}" with "{name}", but failed to translate position = {coords}'
                )

        # simbad ID lookup failed, so now try to read position from
        # coordinates
        REPOS = re.compile(r'J(\d\d)(\d\d)(\d\d\.\d(?:\d*)?)([+-])(\d\d)(\d\d)(\d\d(?:\.\d*)?)$')

        m = REPOS.search(target)
        if m:
            rah,ram,ras,decsgn,decd,decm,decs = m.group(1,2,3,4,5,6,7)
            rah,ram,ras,decd,decm,decs = int(rah),int(ram),float(ras),int(decd),int(decm),float(decs)
            if rah > 23 or ram > 59 or ras >= 60. or decd > 89 or decm > 59 or decs >= 60.:
                err_msgs.append(
                    f'{target} matched the positional regular expression but the numbers were out of range'
                )
            else:
                name = re.sub(RESPC, ' ', target)
                ra = rah+ram/60+ras/3600
                dec = decd+decm/60+decs/3600
                dec = dec if decsgn == '+' else -dec
                return (name,ra,dec)

        err_msgs.append(
            f'Could not extract position from {target}'
        )

    if ra_tel is not None and dec_tel is not None and dist is not None:
        # Second
        query = f"""set limit 2
format object form1 "Target: %IDLIST(1) | %COO(A D;ICRS)"
query coo {ra_tel} {dec_tel} radius={dist}m"""

        payload = {'submit' : 'submit script', 'script' : query}
        rep = requests.post(url, data=payload)

        data = False
        found = False
        fail = False
        nres = 0
        for line in rep.text.split('\n'):
            if line.startswith('::data::'):
                data = True

            if line.startswith('::error::'):
                err_msgs.append(
                    f'Error occured during simbad coordinate query for {target}'
                )
                fail = True
                break

            if data:
                if line.startswith('Target:'):
                    nres += 1
                    if nres > 1:
                        err_msgs.append(
                            f'More than one target returned when querying simbad for telescope position = {ra_tel} {dec_tel} ({target})'
                        )
                        fail = True
                        break

                    name, coords = line[7:].split(' | ')
                    found = True

        if found and not fail:
            # OK we have found one, but we are still not done --
            # some SIMBAD lookups are no good
            name = name.strip()
            name = re.sub(RESPC, ' ', name)
            try:
                ra, dec, syst = str2radec(coords)
                ra_str= dec2sexg(ra,False,2)
                dec_str = dec2sexg(dec,True,1)
                print(f'  Matched telescope position = {coords} {dec} with "{name}", position = {ra_str} {dec_str}')
                return (name,ra,dec)
            except:
                err_msgs.append(
                    f'Matched telescope position = {ra} {dec} ({target}) with "{name}", but failed to translate position = {coords}'
                )

    raise hcam.HipercamError('\n'.join(err_msgs))


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

        if system is None: system = 'ICRS'

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

def dec2sexg(value, sign, dp, sep=':'):
    """
    Converts deciaml to sexagesimal
    """
    v = abs(value)
    v1 = int(v)
    v2 = int(60*(v-v1))
    v3 = 3600*(v-v1-v2/60)
    if sign:
        if value >= 0.:
            return f'+{v1:02d}{sep}{v2:02d}{sep}{v3:0{3+dp}.{dp}f}'
        else:
            return f'-{v1:02d}{sep}{v2:02d}{sep}{v3:0{3+dp}.{dp}f}'
    else:
        return f'{v1:02d}{sep}{v2:02d}{sep}{v3:0{3+dp}.{dp}f}'

def temp_dir():
    """Generates the name of the hipercam temp directory which will be
    used to store temporary files. Creates the directory if it does
    not exist. The temp directory is called "tmp" and is created as a
    subdirectory of the directory containing the command default files
    (".hipercam" by default)

    """

    # Generate temp directory name
    home = os.environ["HOME"]
    hdir = os.environ.get("HIPERCAM_ENV", os.path.join(home,".hipercam"))
    tmpdir = os.path.join(hdir, "tmp")

    # create it (ok if not already there)
    os.makedirs(tmpdir, exist_ok=True)

    return tmpdir

# Used to translate month numbers
LOG_MONTHS = {
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
LOG_CSS = """

body {
    background-color: black;
    font: 11pt sans-serif;
    color: #FFFFFF;
    margin: 10px;
    border: 1px;
    overflow: auto;
    height: 100%;
    max-height: 100%;
}

button {
   font-family: Arial, Helvetica, sans-serif;
   font-size: 11pt;
   padding: 4px 8px;
   text-align: center;
   display: inline-block;
   border-radius: 6px;
   border: 3px solid darkblue;
}

button:hover {
   background-color: #ddd;
}

p {
   color: #ffffe0
}

h1 {
   color: white
}

input.text {
   margin-right: 20px;
   margin-left: 10px
}

span {
   margin-right: 20px;
   margin-left: 20px
}

.tableFixHead {
    overflow-y: auto;
    height: 700px;
}

.tableFixHead thead th {
    position: sticky;
    top: 0;
    background: black;
}

.tableFixHead td {
    box-shadow: inset 2px 2px #888, 2px 2px #888;
}

.tableFixHead th {
    box-shadow: inset 2px 2px #888, 2px 2px #888;
}

table {
    border-collapse: collapse;
    font: 11pt sans-serif;
    padding: 4px;
    margin: 10px;
    border: 2px;
    border-color: white;
}

tr:hover {
    background-color: #303030;
}

td {
    vertical-align: top;
    text-align: center;
    white-space: nowrap;
    padding: 4px;
}

td.left {
    text-align: left;
}

td.lalert {
    color: #FFAAAA;
    text-align: left;
}

td.right {
    text-align: right;
}

td.cen {
    text-align: center;
}

td.undef {
    background-color: #100000;
}

td.format {
    text-align: center;
}

td.long {
    text-align: left;
    font: 9pt sans-serif;
    white-space: normal;
}

td.bleft {
    color: #ffffa0;
    text-align: left;
    font: 9pt sans-serif;
    font-weight: bold;
}

th {
    vertical-align:top;
    text-align: center;
    padding: 4px;
    color: #ffffaa;
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

# end of CSS

