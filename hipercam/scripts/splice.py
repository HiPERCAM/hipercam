import sys
import os
import time

import numpy as np

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ = [
    "splice",
]

##################################################################
#
# mstats -- lists stats of a series of images from a raw data file
#
##################################################################


def splice(args=None):
    """``splice input1 input2 ccd [win] output``

    Merges two multi-CCD images so that parts of input1 are replaced
    by equivalent parts of input2. This can be useful e.g. to create
    a combined flat field out of different frames.

    Parameters:

       input1 : str
           first input frame

       input2 : str
           second input frame

       ccd : str
           the CCD(s) to transfer, '0' for all of them. '1 3', i.e.
           separate with spaces for more than one CCD.

       win : str [hidden]
           the window(s) to transfer, '0' for all of them (default). The
           same windows must exist in the selected CCDs of both inputs 

       output : str
           output file.

    """

    command, args = utils.script_args(args)

    # get inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("input1", Cline.GLOBAL, Cline.PROMPT)
        cl.register("input2", Cline.GLOBAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.PROMPT)
        cl.register("win", Cline.LOCAL, Cline.HIDE)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        infile1 = cl.get_value(
            "input1", "first input frame", cline.Fname("frame1", hcam.HCAM)
        )
        mccd1 = hcam.MCCD.read(infile1)

        infile2 = cl.get_value(
            "input2", "second input frame", cline.Fname("frame2", hcam.HCAM)
        )
        mccd2 = hcam.MCCD.read(infile2)

        if len(mccd1) > 1:
            ccd = cl.get_value("ccd", "CCD(s) to transfer ('0' for all)", "1")
            if ccd == "0":
                ccds = list(mccd1.keys())
            else:
                ccds = ccd.split()
        else:
            ccds = list(mccd1.keys())

        # Find window names in common across all selected CCDs
        # (intersection of sets)
        wins = None
        for cnam in ccds:
            if wins is None:
                wins = set(mccd1[cnam].keys())
            else:
                wins &= set(mccd1[cnam].keys())

        if len(wins) == 0:
            raise hcam.HipercamError(
                "There were no windows in common across the selected CCDs"
            )
        elif len(wins) > 1:
            cl.set_default("win", "0")
            win = cl.get_value("win", "windows to transfer ('0' for all)", "0")
            if win != "0":
                wins = wins.split()

        outfile = cl.get_value(
            "output", "output file", cline.Fname("output", hcam.HCAM, cline.Fname.NEW)
        )

    # Now the (trivial) actual work.

    for cnam in ccds:
        for wnam in wins:
            mccd1[cnam][wnam] = mccd2[cnam][wnam]

    mccd1.write(outfile, True)
