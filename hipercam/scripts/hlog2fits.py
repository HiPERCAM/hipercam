import sys
import os
import struct
import numpy as np
import fitsio

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

###############################################
#
# hlog2fits -- convert reduce log files to FITS
#
###############################################

def hlog2fits(args=None):
    """``hlog2fits log``

    Converts a |hiper| ASCII log into a FITS file. As well as a modest
    reduction in file size (~40%, the ASCII logs are written relatively
    efficiently), the resulting file is much faster to read than the ASCII log
    so this may be useful for very large log files [test of 78,000 frame file:
    12.9 seconds to read the ASCII file, 1.9 to read the FITS version]. At the
    moment no significant header information is transferred beyond the CCD
    names. Each CCD appears as a single binary table, starting at the second
    HDU (or HDU 1 if you number them 0,1,2 ..). This can be read using
    :meth:`hipercam.hlog.Hlog.from_fits`.

    Parameters:

      log : string
         name of the log file (should end .log). The output FITS file
         will have the same root name but end .fits. The routine will abort
         if there is a pre-existing file of the same name.

    """

    command, args = utils.script_args(args)

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('log', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        log = cl.get_value(
            'log', 'name of log file from reduce to convert to FITS',
            cline.Fname('red', hcam.LOG)
        )

    oname = os.path.basename(log)
    oname = oname[:oname.find('.')] + '.fits'
    if os.path.exists(oname):
        raise hcam.HipercamError(
            ('A file called {:s} already exists and'
             ' will not be over-written; aborting').format(oname)
            )

    # Read in the ASCII log
    hlg = hcam.hlog.Hlog.from_ascii(log)

    # Open fits file for output, write out results CCD by CCD as binary tables
    with fitsio.FITS(oname,'rw') as fout:
        for cnam in sorted(hlg):
            print('Writing CCD =',cnam)
            fout.write(
                hlg[cnam], extname='CCD {:s}'.format(cnam),
                header={'CCDNAME' : cnam}
                )

    print('Converted {:s} to {:s}'.format(log,oname))

