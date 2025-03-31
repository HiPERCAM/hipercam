import sys
import os
import struct
import numpy as np
from astropy.io import fits

from trm import cline
from trm.cline import Cline

import hipercam as hcam

__all__ = [
    "hlog2fits",
]

###############################################
#
# hlog2fits -- convert reduce log files to FITS
#
###############################################


def hlog2fits(args=None):
    """``hlog2fits log [origin dir]``

    Converts a |hiper| ASCII log into a FITS file. As well as a modest
    reduction in file size (~40%, the ASCII logs are written relatively
    efficiently), the resulting file is faster to read than the ASCII log
    so this may be useful for very large log files [test of 78,000 frame file:
    12.9 seconds to read the ASCII file, 1.9 to read the FITS version]. The FITS
    log is also much easier to understand than the ASCII files, but they don't
    have all the header information, so are not a replacement. At the
    moment no significant header information is transferred beyond the CCD
    names. Each CCD appears as a single binary table, starting at the second
    HDU (or HDU 1 if you number them 0,1,2 ..). This can be read using
    :meth:`hipercam.hlog.Hlog.from_fits`.

    Parameters:

      log : str
         name of the log file (should end .log). The output FITS file
         will have the same root name but end .fits. The routine will abort
         if there is a pre-existing file of the same name.

      origin : str [hidden]
         'h' or 'u' depending upon whether the log file was created with
         the hipercam or old ultracam pipeline. Defaults to 'h'.

      dir : str [hidden]
         directory for the output; defaults to the present working directory

    NB Because of the danger of over-writing raw data (also ends
    .fits), this routine will not over-write pre-existing files. You
    should delete clashing files if you really want to proceed.

    """

    command, args = cline.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("log", Cline.LOCAL, Cline.PROMPT)
        cl.register("origin", Cline.LOCAL, Cline.HIDE)
        cl.register("dir", Cline.LOCAL, Cline.HIDE)

        # get inputs
        log = cl.get_value(
            "log",
            'name of log file from "reduce" to convert to FITS',
            cline.Fname("red", hcam.LOG),
        )

        cl.set_default('origin','h')
        origin = cl.get_value(
            "origin", "h(ipercam) or u(ltracam) pipeline?", "h",
            lvals=["h", "u"]
        )

        cl.set_default('dir','.')
        dir = cl.get_value(
            "dir", "directory for output", ".",
        )

    root = os.path.splitext(os.path.basename(log))[0]
    oname = os.path.join(dir, root + ".fits")
    if os.path.exists(oname):
        raise hcam.HipercamError(
            f"A file called {oname} already exists and"
            " will not be over-written; aborting"
        )

    # Read in the ASCII log
    if origin == "h":
        hlg = hcam.hlog.Hlog.rascii(log)
    elif origin == "u":
        hlg = hcam.hlog.Hlog.rulog(log)

    print(f"Loaded ASCII log = {log}")

    # Generate HDU list

    # First the primary HDU (no data)
    phdr = fits.Header()
    phdr["LOGFILE"] = (os.path.basename(log), "Original log file")
    phdu = fits.PrimaryHDU(header=phdr)
    hdul = [
        phdu,
    ]

    # Now a BinTable for each CCD
    for cnam in sorted(hlg):
        hdr = fits.Header()
        hdr["CCDNAME"] = (cnam, "CCD name")
        hdul.append(
            fits.BinTableHDU(
                hlg[cnam], header=hdr,
                name=f"CCD {cnam}"
            )
        )

    hdul = fits.HDUList(hdul)

    # finally write to disk
    print(f"Writing to disk in file = {oname}")
    hdul.writeto(oname)

    print(f"Converted {log} to {oname}")
