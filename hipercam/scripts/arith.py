import sys
import os

import numpy as np

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = ["add", "div", "mul", "sub"]

#############################################
#
# arith -- arithematic with multi-CCD images
#
#############################################


def arith(args=None):
    """
    Carries out operations of the form output = input1 [op] input2 where [op]
    is '+', '-', '*' or '/' referred to by 'add', 'sub', 'mul' and
    'div'.

    Parameters:

       input1 : string
          first input hcm file

       input2 : string
          second input hcm file

       ccd : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. '2 4' or just
          one '3'

       win : string [hidden, defaults to 'all']
          the CCD or CCDs to apply the operation to. 'all' for the whole lot
          which it returns to by default.  Can be several e.g. 'E2 G1' or just
          one 'H1'. If you specify windows in this manner, it is assumed that
          all the CCDs chosen in the previous input have the named windows;
          'all' just applies the operation to all windows regardless.

       crop : bool [hidden, defaults to False]
         set True to try to crop input2 to have the same format as input1 (it
         must enclose it and have compatible binning for this to work).

       output : string
          output hcm file name. Can be same as either input1 or input2
          in which case the input file will be over-written.

    """

    command, args = utils.script_args(args)

    # get inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("input1", Cline.LOCAL, Cline.PROMPT)
        cl.register("input2", Cline.LOCAL, Cline.PROMPT)
        cl.register("ccd", Cline.LOCAL, Cline.HIDE)
        cl.register("win", Cline.LOCAL, Cline.HIDE)
        cl.register("crop", Cline.LOCAL, Cline.HIDE)
        cl.register("output", Cline.LOCAL, Cline.PROMPT)

        prompts = {
            "add": "add",
            "sub": "subtract",
            "mul": "multiply by",
            "div": "divide by",
        }

        infile1 = cl.get_value(
            "input1", "first input file", cline.Fname("hcam", hcam.HCAM)
        )
        mccd1 = hcam.MCCD.read(infile1)

        infile2 = cl.get_value(
            "input2", "second input file", cline.Fname("hcam", hcam.HCAM)
        )
        mccd2 = hcam.MCCD.read(infile2)

        if len(mccd1) > 1:
            cl.set_default("ccd", "all")
            ccd = cl.get_value("ccd", "CCD(s) to process", "all")
            if ccd == "all":
                ccds = list(mccd1.keys())
            else:
                ccds = ccd.split()
        else:
            ccd = "all"
            ccds = list(mccd1.keys())

        tccd = mccd1[ccds[0]]
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

        cl.set_default("crop", False)
        crop = cl.get_value(
            "crop", "try to crop input2 to the same format as input1", False
        )
        if crop:
            mccd2 = mccd2.crop(mccd1)
            print(
                "cropped {:s} to match {:s} before operation".format(infile2, infile1)
            )

        outfile = cl.get_value(
            "output", "output file", cline.Fname("hcam", hcam.HCAM, cline.Fname.NEW)
        )

    # carry out operation
    if command == "add":
        # addition
        for cnam in ccds:
            ccd1 = mccd1[cnam]
            ccd2 = mccd2[cnam]
            if wins == "all":
                ccd1 += ccd2
            else:
                for wnam in wins:
                    ccd1[wnam] += ccd2[wnam]

    elif command == "sub":
        # subtraction
        for cnam in ccds:
            ccd1 = mccd1[cnam]
            ccd2 = mccd2[cnam]
            if wins == "all":
                ccd1 -= ccd2
            else:
                for wnam in wins:
                    ccd1[wnam] -= ccd2[wnam]

    elif command == "mul":
        # multiplication
        for cnam in ccds:
            ccd1 = mccd1[cnam]
            ccd2 = mccd2[cnam]
            if wins == "all":
                ccd1 *= ccd2
            else:
                for wnam in wins:
                    ccd1[wnam] *= ccd2[wnam]

    elif command == "div":
        # multiplication
        for cnam in ccds:
            ccd1 = mccd1[cnam]
            ccd2 = mccd2[cnam]
            if wins == "all":
                ccd1 /= ccd2
            else:
                for wnam in wins:
                    ccd1[wnam] /= ccd2[wnam]

    # Add a history line
    mccd1.head.add_history(
        "{:s} {:s} {:s} {:s} {:s} {:s}".format(
            command,
            utils.sub_extension(infile1, hcam.HCAM),
            utils.sub_extension(infile2, hcam.HCAM),
            ccd,
            win,
            utils.sub_extension(outfile, hcam.HCAM),
        )
    )

    # save result
    mccd1.write(outfile, True)


def add(args=None):
    """``add input1 input2 [ccd win] output``

    Adds two hcm frames and outputs the result. Can be applied only to
    particular CCDs and windows if wanted.

    Parameters:

       input1  : string
          first input hcm file

       input2  : string
          second input hcm file to add to the first.

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
          output hcm file name. Can be same as either input1 or input2
          in which case the input file will be over-written.
    """
    arith(args)


def sub(args=None):
    """``sub input1 input2 [ccd win] output``

    Subtracts two hcm frames and outputs the result. Can be applied only to
    particular CCDs and windows if wanted.

    Parameters:

       input1  : string
          first input hcm file

       input2  : string
          second input hcm file to subtract from the first.

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
          output hcm file name. Can be same as either input1 or input2
          in which case the input file will be over-written.

    """
    arith(args)


def div(args=None):
    """``div input1 input2 [ccd win] output``

    Divides two hcm frames and outputs the result. Can be applied only to
    particular CCDs and windows if wanted.

    Parameters:

       input1  : string
          first input hcm file

       input2  : string
          second input hcm file to divide into the first.

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
          output hcm file name. Can be same as either input1 or input2
          in which case the input file will be over-written.

    """
    arith(args)


def mul(args=None):
    """``mul input1 input2 [ccd win] output``

    Multiplies two hcm frames and outputs the result. Can be applied only to
    particular CCDs and windows if wanted.

    Parameters:

       input1  : string
          first input hcm file

       input2  : string
          second input hcm file to multiply into the first.

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
          output hcm file name. Can be same as either input1 or input2
          in which case the input file will be over-written.

    """
    arith(args)
