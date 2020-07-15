import sys
import os

import numpy as np

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = ["cadd", "csub", "cdiv", "cmul"]

#############################################
#
# carith -- arithematic with multi-CCD images
#
#############################################


def carith(args=None):
    """Carries out operations of the form output = input [op] constant where [op]
    is '+', '-', '*' or '/' referred to by 'cadd', 'csub', 'cmul' and
    'cdiv'. This is meant to be used as an entry point script. If called
    within a program, then the first argument should be 'cadd', 'csub', 'cmul'
    or 'cdiv' to define the default parameter file name.

    Parameters:

       input   : string
          input hcm file name

       const   : float
          constant number to apply

       ccd     : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. '2 4' or just
          one '3'

       win     : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. 'E2 G1' or just
          one 'H1'. If you specify windows in this manner, it is assumed that
          all the CCDs chosen in the previous input have the named windows;
          'all' just applies the operation to all windows regardless.

       output  : string
          output hcm file name. Can be same as input in which case
          the file will be over-written.

    """

    command, args = utils.script_args(args)

    # get inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("input", Cline.LOCAL, Cline.PROMPT)
        cl.register("const", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.HIDE)
        cl.register("win", Cline.LOCAL, Cline.HIDE)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        prompts = {
            "cadd": "add",
            "csub": "subtract",
            "cmul": "multiply by",
            "cdiv": "divide by",
        }

        infile = cl.get_value("input", "input file", cline.Fname("hcam", hcam.HCAM))
        mccd = hcam.MCCD.read(infile)

        constant = cl.get_value("const", "constant to " + prompts[command], 0.0)

        if len(mccd) > 1:
            cl.set_default("ccd", "all")
            ccd = cl.get_value("ccd", "CCD(s) to process", "all")
            if ccd == "all":
                ccds = list(mccd.keys())
            else:
                ccds = ccd.split()
        else:
            ccd = "all"
            ccds = list(mccd.keys())

        tccd = mccd[ccds[0]]
        if len(tccd) > 1:
            cl.set_default("win", "all")
            win = cl.get_value("win", "window(s) to process", "all")
            if win == "all":
                wins = "all"
            else:
                wins = win.split()
        else:
            win = "all"
            wins = "all"

        outfile = cl.get_value(
            "output", "output file", cline.Fname("hcam", hcam.HCAM, cline.Fname.NEW)
        )

    # carry out operation
    if command == "cadd":
        # addition
        for cnam in ccds:
            tccd = mccd[cnam]
            if wins == "all":
                tccd += constant
            else:
                for wnam in wins:
                    tccd[wnam] += constant

    elif command == "csub":
        # subtraction
        for cnam in ccds:
            tccd = mccd[cnam]
            if wins == "all":
                tccd -= constant
            else:
                for wnam in wins:
                    tccd[wnam] -= constant

    elif command == "cmul":
        # multiplication
        for cnam in ccds:
            tccd = mccd[cnam]
            if wins == "all":
                tccd *= constant
            else:
                for wnam in wins:
                    tccd[wnam] *= constant

    elif command == "cdiv":
        # division
        for cnam in ccds:
            tccd = mccd[cnam]
            if wins == "all":
                tccd /= constant
            else:
                for wnam in wins:
                    tccd[wnam] /= constant

    # Add a history line
    mccd.head.add_history(
        "{:s} {:s} {:f} {:s} {:s} {:s}".format(
            command,
            utils.sub_extension(infile, hcam.HCAM),
            constant,
            ccd,
            win,
            utils.sub_extension(outfile, hcam.HCAM),
        )
    )

    # save result
    mccd.write(outfile, True)


def cadd(args=None):
    """``cadd input const [ccd win] output``

    Adds a constant to a HiPERCAM hcm frame. Can be applied only to particular
    CCDs and windows if wanted.

    Parameters:

       input  : string
          input hcm file name

       const  : float
          constant to add

       ccd    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. '2 4' or just
          one '3'

       win    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. 'E2 G1' or just
          one 'H1'. If you specify windows in this manner, it is assumed that
          all the CCDs chosen in the previous input have the named windows;
          'all' just applies the operation to all windows regardless.

       output : string
          output hcm file name. Can be same as input in which case
          the file will be over-written.

    """

    # send to carith
    carith(args)


def csub(args=None):
    """``csub input const [ccd win] output``

    Subtracts a constant from a HiPERCAM hcm frame. Can be applied only to
    particular CCDs and windows if wanted.

    Parameters:

       input  : string
          input hcm file name

       const  : float
          constant to subtract

       ccd    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. '2 4' or just
          one '3'

       win    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. 'E2 G1' or just
          one 'H1'. If you specify windows in this manner, it is assumed that
          all the CCDs chosen in the previous input have the named windows;
          'all' just applies the operation to all windows regardless.

       output : string
          output hcm file name. Can be same as input in which case
          the file will be over-written.

    """

    # send to carith
    carith(args)


def cdiv(args=None):
    """``cdiv input const [ccd win] output``

    Divides a HiPERCAM hcm frame by a constant. Can be applied only to
    particular CCDs and windows if wanted.

    Parameters:

       input  : string
          input hcm file name

       const  : float
          constant to divide by

       ccd    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. '2 4' or just
          one '3'

       win    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. 'E2 G1' or just
          one 'H1'. If you specify windows in this manner, it is assumed that
          all the CCDs chosen in the previous input have the named windows;
          'all' just applies the operation to all windows regardless.

       output : string
          output hcm file name. Can be same as input in which case
          the file will be over-written.

    """

    # send to carith
    carith(args)


def cmul(args=None):
    """``cmul input const [ccd win] output``

    Multiplies a HiPERCAM hcm frame by a constant. Can be applied only to
    particular CCDs and windows if wanted.

    Parameters:

       input  : string
          input hcm file name

       const  : float
          constant to multiply by

       ccd    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. '2 4' or just
          one '3'

       win    : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. 'E2 G1' or just
          one 'H1'. If you specify windows in this manner, it is assumed that
          all the CCDs chosen in the previous input have the named windows;
          'all' just applies the operation to all windows regardless.

       output : string
          output hcm file name. Can be same as input in which case
          the file will be over-written.

    """

    # send to carith
    carith(args)
