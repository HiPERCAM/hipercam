import sys
import os

import numpy as np

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = ['stats',]

#################################################################
#
# stats -- lists basic stats of each window of a multi-CCD image.
#
#################################################################

def stats(args=None):
    """Lists basic stats of a multi-CCD image, i.e. the minimum, maximum,
    mean, median and standard deviation of each window of each CCD. The output
    format can be altered to suit preference.

    Arguments::

      input  : (string)
         name of the MCCD file

      format : (string) [hidden, defaults to 9.3f]
         C-style format code as used in Python format statements for output of
         the numerical values. e.g. '300.00' is '6.2f' (6 characters toal, 2 after
         the decimal point), '1.22e24' is '.2e' (as many characters as needed, 2
         after the decimal point)

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('input', Cline.LOCAL, Cline.PROMPT)
        cl.register('format', Cline.LOCAL, Cline.HIDE)

        # get inputs
        frame = cl.get_value('input', 'frame to lists stats of',
                             cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.read(frame)

        cl.set_default('format','9.3f')
        form = cl.get_value('format', 'output format for numbers', '9.3f')

    for cnam, ccd in mccd.items():
        for wnam, wind in ccd.items():
            print(
                'CCD {0:s}, window {1:s}: min = {3:{2:s}}, max = {4:{2:s}}, mean = {5:{2:s}}, median = {6:{2:s}}, std = {7:{2:s}}'.format(
                    cnam, wnam, form, wind.min(), wind.max(), wind.mean(), wind.median(), wind.std()
                    )
                )


