import sys
import traceback
import os
import shutil
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import astropy.units as u

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline
from hipercam.utils import format_ulogger_table, target_lookup

__all__ = [
    "logsearch",
]

#######################################################################
#
# logsearch -- carries out a search for particular objects in data logs
#
#######################################################################

def logsearch(args=None):
    description = \
    """``logsearch target tmin [dmax regex (nocase) instrument hcamlog ucamlog uspeclog]``

    Searches for ``target`` in the |hiper| or |ucam| logs. Carries out
    a coordinate lookup given a name. Then searches for runs with
    matching RAs and Decs, either as identified from target names or
    from telescope pointing values (by ending J123456.7-123456 or
    similar). Can also search by regular expression matching. It tries
    not to miss things, i.e. be on the generous side, so it will
    return a match if either of two sets of coordinates match, even if
    one does not.

    Arguments::

       target : str
          Target name. On the command line, must be enclosed in quotes if it
          contains space. /this will be used first to carry out a lookup in
          SIMBAD to find the RA and Dec. Failing this it tries to identify
          coordinates from a final strength of the form JHHMMSS.S[+/-]DDMMSS

       tmin : float
          Minimum exposure duration seconds.

       dmax : float [hidden]
          Maximum distance from lookup position, arcminutes

       regex : str [hidden]
          Regular expression to use to try to match target names in addition to
          the coordinate matching. Default to "none"

       nocase : bool [hidden, if regex is not "none"]
          True for case-insensitive matching, else case-sensitive used with regex

       instrument : str [hidden]
          'hipercam', 'ultracam', or 'ultrapec'

       hcamlog : log
          Path to hipercam spreadsheet

       ucamlog : log
          Path to ultracam spreadsheet

       uspeclog : log
          Path to ultraspec spreadsheet

    """

    command, args = utils.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("target", Cline.LOCAL, Cline.PROMPT)
        cl.register("tmin", Cline.LOCAL, Cline.PROMPT)
        cl.register("dmax", Cline.LOCAL, Cline.HIDE)
        cl.register("regex", Cline.LOCAL, Cline.HIDE)
        cl.register("nocase", Cline.LOCAL, Cline.HIDE)
        cl.register("instrument", Cline.LOCAL, Cline.HIDE)
        cl.register("hcamlog", Cline.LOCAL, Cline.HIDE)
        cl.register("ucamlog", Cline.LOCAL, Cline.HIDE)
        cl.register("uspeclog", Cline.LOCAL, Cline.HIDE)

        # get inputs
        target = cl.get_value(
            "target", "target name for simbad lookup",
            "AR Sco"
        )

        tmin = cl.get_value(
            "tmin", "minimum exposure duration for a run to be included [seconds, < 0 for no check at all]", -1
        )

        dmax = cl.get_value(
            "dmax", "maximum distance from target [arcminutes]",
            12., 0.
        )

        cl.set_default('regex','none')
        regex = cl.get_value(
            "regex", "regular expression to match target name [none to ignore]",
            "none", ignore="none"
        )

        if regex is not None:
            nocase = cl.get_value(
                "nocase", "case insensitive match?", True
            )

        instrument = cl.get_value(
            "instrument",
            "instrument of interest [hipercam, ultracam, ultraspec]",
            "ultraspec", lvals=['hipercam', 'ultracam', 'ultraspec']
        )

        if instrument == 'hipercam':
            log = cl.get_value(
                "hcamlog", "path to hipercam spreadsheet",
                cline.Fname("hipercam.xlsx", ".xlsx")
            )

        elif instrument == 'ultracam':
            log = cl.get_value(
                "ucamlog", "path to ultracam spreadsheet",
                cline.Fname("ultracam.xlsx", ".xlsx")
            )

        elif instrument == 'ultraspec':
            log = cl.get_value(
                "uspeclog", "path to ultraspec spreadsheet",
                cline.Fname("ultraspec.xlsx", ".xlsx")
            )

    name, ra, dec = target_lookup(target)
    if name == 'UNDEF':
        print(f'Coordinate lookup for "{target}" failed')
        return

    print(f'Coordinate lookup for "{target}" returned name = "{name}", RA [hrs] = {ra}, Dec [deg] = {dec}')

    ra = np.radians(15*ra)
    dec = np.radians(dec)
    cat = np.cos(ra)
    sat = np.sin(ra)
    cdt = np.cos(dec)
    sdt = np.sin(dec)
    cmin = np.cos(np.radians(dmax/60))

    print(f'\nLoading log = {log}')
    tab = pd.read_excel(log)
    print(f'Loaded {len(tab)} runs from {log}')

    # Extract positions
    ras = tab['ra_deg']
    decs = tab['dec_deg']
    ras = np.array(np.radians(ras))
    decs = np.array(np.radians(decs))

    ca = np.cos(ras)
    sa = np.sin(ras)
    cd = np.cos(decs)
    sd = np.sin(decs)
    cdiff = ca*cd*cat*cdt + cd*sa*cdt*sat + sd*sdt
    # set logical array all False
    ok = np.zeros_like(cdiff, dtype=np.bool)
    good = ~np.isnan(cdiff)
    ok[good] = cdiff[good] > cmin

    if instrument == 'ultraspec':
        ras_tels = np.array(np.radians(tab['ra_tel_deg']))
        decs_tels = np.array(np.radians(tab['dec_tel_deg']))
        ca = np.cos(ras_tels)
        sa = np.sin(ras_tels)
        cd = np.cos(decs_tels)
        sd = np.sin(decs_tels)
        cdiff = ca*cd*cat*cdt + cd*sa*cdt*sat + sd*sdt
        good = ~np.isnan(cdiff)
        ok[good] |= cdiff[good] > cmin

    if regex is not None:
        # regular expression matching
        rec = re.compile(regex)
        targets = tab['target']
        for n, targ in enumerate(targets):
            if isinstance(targ, str):
                m = rec.match(targ, re.I) \
                    if nocase else \
                       rec.match(targ)
                if m: ok[n] = True

    tab = tab[ok]

    output = f"{target.replace(' ','_')}.xlsx"
    format_ulogger_table(output, tab, instrument)
    print(f'Wrote {len(tab)} runs to {output}')
