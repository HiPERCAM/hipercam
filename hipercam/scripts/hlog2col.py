import sys
import os
import struct
import numpy as np
from astropy.time import Time

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = [
    "hlog2col",
]

#######################################################
#
# hlog2col -- converts reduce log files to columns data
#
#######################################################


def hlog2col(args=None):
    """``hlog2col log ap1 ap2``

    Converts a |hiper| ASCII log into one or more column format files
    of one aperture versus another for each CCD. This tacks on the date
    of the start of the night in order to distinguish run log names of the
    the same number.

    Parameters:

      log : str
         name of the log file (should end .log). The output FITS file
         will have the same root name but end .fits. The routine will
         abort if there is a pre-existing file of the same name.

      origin : str
         'h' or 'u' depending upon whether the log file was created with
         the hipercam or old ultracam pipeline.

      ap1 : str
         name of first aperture

      ap2 : str
         name of second aperture

      type : str
         output type 'l' for linear; 'm' for magnitudes. Linear has the
         advantage of handling negative values.

    The output files will get names like 2021-01-12_run014_2_3_4.dat meaning
    CCD 2, aperture 3 divided by 4 of run014 from the night starting 2021-01-12.
    """

    command, args = utils.script_args(args)

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("log", Cline.LOCAL, Cline.PROMPT)
        cl.register("origin", Cline.LOCAL, Cline.PROMPT)
        cl.register("ap1", Cline.LOCAL, Cline.PROMPT)
        cl.register("ap2", Cline.LOCAL, Cline.PROMPT)
        cl.register("type", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        log = cl.get_value(
            "log",
            'name of log file from "reduce" to convert to ASCII column data',
            cline.Fname("run010", hcam.LOG),
        )

        origin = cl.get_value(
            "origin", "h(ipercam) or u(ltracam) pipeline?", "h", lvals=["h", "u"]
        )

        ap1 = cl.get_value(
            "ap1", "first aperture", "1"
        )

        ap2 = cl.get_value(
            "ap2", "second aperture", "2"
        )

        otype = cl.get_value(
            "type", "l(inear) or m(agnitudes)?", "l", lvals=["l", "m"]
        )


    # Read in the ASCII log
    if origin == "h":
        hlg = hcam.hlog.Hlog.read(log)
    elif origin == "u":
        hlg = hcam.hlog.Hlog.fulog(log)

    print(f"Loaded ASCII log = {log}")

    bname = os.path.basename(log[:-4])
    for cnam in hlg.cnames:
        apnames = hlg.apnames[cnam]
        if ap1 in apnames and ap2 in apnames:
            lc = hlg.tseries(cnam, ap1) / hlg.tseries(cnam, ap2)
            time = Time(round(lc.t[0]-0.5), format='mjd', scale='utc')
            start = time.iso[:10]
            oname = f'{start}_{bname}_{cnam}_{ap1}_{ap2}.asc'

            if otype == 'm':
                lc.to_mag()
                fmt = '%.10f %.3e %.4f %.4f %d'
                header = \
                    f"""Data derived from {log}, taken on night starting {start}.

CCD = {cnam}, mag (aperture {ap1}) - mag(aperture {ap2})

Five columns:

MJD (centre of exposure), exposure time (sec), mags, uncertainty (mags), flag
                    """

            else:
                fmt = '%.10f %.3e %.4e %.3e %d'
                header = \
                    f"""Data derived from {log}, taken on night starting {start}.

CCD = {cnam}, aperture {ap1} / aperture {ap2}

Five columns:

MJD (centre of exposure), exposure time (sec), linear ratio, uncertainty, flag
                    """


            data = np.column_stack([lc.t,86400.*lc.te,lc.y,lc.ye,lc.bmask])
            np.savetxt(oname,data,fmt,header=header)
