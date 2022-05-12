import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from trm.pgplot import *

from trm import cline
from trm.cline import Cline

import hipercam as hcam

##################################################
#
# redanal -- some basic stats on a reduce log file
#
##################################################


def redanal(args=None):
    """``redanal aper log``

    This provides some basic stats on a log file such as the fraction of duff
    points for each aperture of each CCD, how much targets have moved and the
    like. The aim is to help with setting parameters in the reduce file in
    difficult cases and to diagnose problems. It will probably be added to with
    time.

    Parameters:

      aper : string
         the aperture file used for the reduction.

      log : string
         the log file.

    """

    command, args = cline.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("aper", Cline.LOCAL, Cline.PROMPT)
        cl.register("log", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        aper_file = cl.get_value(
            "aper", "aperture file used for reduction", cline.Fname("run001", hcam.APER)
        )
        aper = hcam.MccdAper.read(aper_file)

        log_file = cl.get_value(
            "log",
            "ASCII reduction log file to analyse",
            cline.Fname("run001", hcam.LOG),
        )
        log = hcam.hlog.Hlog.rascii(log_file)

    for cnam in sorted(log):
        print()
        mjds = log[cnam]["MJD"]
        mjdoks = log[cnam]["MJDok"]
        diffs = mjds[1:] - mjds[:-1]
        tgap_max = diffs.max()
        tgap_mean = diffs.mean()
        print(
            "CCD {:s} mean / maximum time gap = {:.2f} / {:.2f} seconds".format(
                cnam, 86400.0 * tgap_mean, 86400.0 * tgap_max
            )
        )
        apnams = log.apnames[cnam]
        for apnam in apnams:
            if not aper[cnam][apnam].linked:
                x = log[cnam]["x_{:s}".format(apnam)]
                xe = log[cnam]["xe_{:s}".format(apnam)]
                y = log[cnam]["y_{:s}".format(apnam)]
                ye = log[cnam]["ye_{:s}".format(apnam)]
                ok = (xe > 0) & (ye > 0)
                xdiffs = x[ok][1:] - x[ok][:-1]
                ydiffs = y[ok][1:] - y[ok][:-1]
                jitt = np.sqrt(xdiffs ** 2 + ydiffs ** 2)
                print(
                    "CCD {:s}, ap {:s} has {:d}/{:d}/{:d} OK/NOK/total points. Min/mean/median/max jitter = {:.2f}/{:.2f}/{:.2f}/{:.2f} pixels {:s}".format(
                        cnam,
                        apnam,
                        len(x[ok]),
                        len(x) - len(x[ok]),
                        len(x),
                        jitt.min(),
                        jitt.mean(),
                        np.median(jitt),
                        jitt.max(),
                        "[reference]" if aper[cnam][apnam].ref else "[non-reference]",
                    )
                )
            else:
                print(
                    "CCD {:s}, aperture {:s} is linked to aperture {:s}".format(
                        cnam, apnam, aper[cnam][apnam].link
                    )
                )
