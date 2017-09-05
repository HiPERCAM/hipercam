import sys
import os

import numpy as np

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ = ['carith',]

#############################################
#
# carith -- arithematic with multi-CCD images
#
#############################################

def carith(args=None):
    """
    Carries out operations of the form output = input [op] constant where [op]
    is '+', '-', '*' or '/' referred to by 'cadd', 'csub', 'cmul' and
    'cdiv'. This is meant to be used as an entry point script. If called
    within a program, then the first argument should be 'cadd', 'csub', 'cmul'
    or 'cdiv' to define the default parameter file name.

    Arguments::

       input  : (string)
          input hcm file name

       const  : (float)
          constant number to apply

       output  : (string)
          output hcm file name. Can be same as input in which case
          the file will be over-written.

    """

    if args is None:
        args = sys.argv
    command = os.path.split(args.pop(0))[1]

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('input', Cline.LOCAL, Cline.PROMPT)
        cl.register('const', Cline.LOCAL, Cline.PROMPT)
        cl.register('output', Cline.LOCAL, Cline.PROMPT)

        prompts = {'cadd' : 'add', 'csub' : 'subtract',
                   'cmul' : 'multiply by', 'cdiv' : 'divide by'}

        infile = cl.get_value('input', 'input file',
                              cline.Fname('hcam', hcam.HCAM))
        mccd = hcam.MCCD.rfits(infile)

        constant = cl.get_value('const', 'constant to ' + prompts[command], 0.)

        outfile = cl.get_value('output', 'output file',
                               cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW))

    # carry out operation
    if command == 'cadd':
        mccd += constant
    elif command == 'csub':
        mccd -= constant
    elif command == 'cmul':
        mccd *= constant
    elif command == 'cdiv':
        mccd /= constant

    # Add a history line
    mccd.head.add_history(
        '{:s} {:s} {:f} {:s}'.format(command,infile,constant,outfile)
    )

    # save result
    mccd.wfits(outfile, True)

