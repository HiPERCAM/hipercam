import sys
import os

import numpy as np

import hipercam as hcam
import hipercam.cline as cline
from hipercam.cline import Cline

__all__ = ['arith',]

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

    Arguments::

       input1  : (string)
          first input hcm file

       input2  : (string)
          second input hcm file

       output  : (string)
          output hcm file name. Can be same as either input1 or input2
          in which case the input file will be over-written.

    """

    if args is None:
        args = sys.argv
    command = os.path.split(args.pop(0))[1]

    # get inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('input1', Cline.LOCAL, Cline.PROMPT)
        cl.register('input2', Cline.LOCAL, Cline.PROMPT)
        cl.register('output', Cline.LOCAL, Cline.PROMPT)

        prompts = {'add' : 'add', 'sub' : 'subtract',
                   'mul' : 'multiply by', 'div' : 'divide by'}

        infile1 = cl.get_value('input1', 'first input file',
                               cline.Fname('hcam', hcam.HCAM))
        mccd1 = hcam.MCCD.rfits(infile1)

        infile2 = cl.get_value('input2', 'second input file',
                               cline.Fname('hcam', hcam.HCAM))
        mccd2 = hcam.MCCD.rfits(infile2)

        outfile = cl.get_value('output', 'output file',
                               cline.Fname('hcam', hcam.HCAM, cline.Fname.NEW))

    # carry out operation
    if command == 'add':
        mccd1 += mccd2
    elif command == 'sub':
        mccd1 -= mccd2
    elif command == 'mul':
        mccd1 *= mccd2
    elif command == 'div':
        mccd1 /= mccd2

    # Add a history line
    mccd1.head.add_history(
        '{:s} {:s} {:s} {:s}'.format(command,infile1,infile2,outfile)
    )

    # save result
    mccd1.wfits(outfile, True)

