import sys
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from trm import cline
from trm.cline import Cline

import hipercam as hcam

__all__ = [
    "hinfo",
]

#########################################################
#
# hinfo -- prints essential information about an hcm file
#
#########################################################


def hinfo(args=None):
    """``hinfo input``

    Prints out an hcm file along with its name and the number of CCDs. You
    will get information CCD-by-CCD first followed by the top level header.
    You will also get at least partial printouts of the data arrays. If
    you want still more detail, I recommend the FITS viewer tool 'fv', or
    'ds9' for examining images (also :hplot:).

    Parameters:

      input  : string
         name of the MCCD file

    """

    command, args = cline.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("input", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        frame = cl.get_value("input", "frame to plot", cline.Fname("hcam", hcam.HCAM))
        mccd = hcam.MCCD.read(frame)

    print("Name of file = {:s}".format(frame))
    print("Number of CCDs = {:d}\n".format(len(mccd)))
    print(mccd)
