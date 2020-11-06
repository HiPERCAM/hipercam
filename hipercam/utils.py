"""
Classes and functions of general use
"""

import os
import sys
import math
import numpy as np
import pandas as pd
from .core import *

__all__ = ("Vec2D", "add_extension", "sub_extension", "script_args", "rgb", "format_hlogger_table")


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
        return CIS[CNAMS[cnam]]


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


def format_hlogger_table(fname, table):
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
        return min(150,int(math.ceil(0.88*max_len+1.5)))

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

