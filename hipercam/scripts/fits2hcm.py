import sys
import os

import numpy as np
from astropy.io import fits
from astropy.time import Time

from trm import cline
from trm.cline import Cline

import hipercam as hcam
from hipercam.core import *

############################################
#
# fits2hcm -- convert non-native FITS to hcm
#
############################################


def fits2hcm(args=None):
    """``fits2hcm flist origin``

    This routine converts "foreign" data into a format suitable for the
    pipeline. |hiper|'s hcm files are in fact FITS-format so this is mostly
    a case of re-organising the files.

    Parameters:

      flist : str
         name of list of input FITS-format files.  'flist' should end
         '.lis'. The output file names will have the same rootname as
         the input files (but without any leading directories) but end
         in '.hcm'

      origin : str
         origin of the data. Currently recognised:

           HICKS :
             University of Sheffield 16" Hicks Telescope
             SBIG ST10-XME CCD

           INTWFC :
             Wide field Camera on the INT. Just operates on a single CCD's-worth
             of data.

           LCOGT :
             Las Cumbres Observatory, based on some FTN data taken in 2019.

           LTIO :
             Liverpool Telescope IO camera.

           LTRISE :
             Liverpool Telescope RISE camera. I don't know how general my
             code is in this case. I assume XBIN=YBIN=2.

           PT5M :
             University of Sheffield/Durham pt5m telescope on La Palma
             QSI 532 CCD

           ROSA :
             University of Sheffield 10" Hicks Telescope
             ATIK ONE 6.0 CCD

           WHTCAM :
             Camera used by Richard Ashley in Feb 2021.

      overwrite : bool
         overwrite files on output

    """

    command, args = cline.script_args(args)

    FORMATS = ["HICKS", "INTWFC", "LTRISE", "LTIO", "PT5M", "ROSA", "LCOGT", "WHTCAM"]

    # get input section
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("flist", Cline.LOCAL, Cline.PROMPT)
        cl.register("origin", Cline.LOCAL, Cline.PROMPT)
        cl.register("overwrite", Cline.LOCAL, Cline.PROMPT)

        # get inputs
        flist = cl.get_value(
            "flist",
            "list of input FITS files to convert to hcm",
            cline.Fname("files", hcam.LIST),
        )

        origin = cl.get_value("origin", "origin of data", "LTRISE", lvals=FORMATS)

        overwrite = cl.get_value("overwrite", "overwrite data on output", True)

    with open(flist) as fin:
        counter = 0
        for line in fin:
            fname = line.strip()
            counter += 1
            with fits.open(fname) as hdul:
                bname = os.path.basename(fname)
                if bname.find(".") > -1:
                    oname = bname[: bname.find(".")] + hcam.HCAM
                else:
                    oname = bname + hcam.HCAM

                if origin == "LTRISE":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "Number of CCDs")
                    ophdu.header["TIMSTAMP"] = ihead["DATE-OBS"]

                    # Copy data into first HDU
                    ofhdu = fits.ImageHDU(hdul[0].data)

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")
                    ofhdu.header["NXTOT"] = (1026, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (1026, "Total unbinned X dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (
                        ihead["CCDWXOFF"] + 1,
                        "X-ordinate of lower-left pixel",
                    )
                    ofhdu.header["LLY"] = (
                        ihead["CCDWYOFF"] + 1,
                        "Y-ordinate of lower-left pixel",
                    )
                    ofhdu.header["XBIN"] = (2, "X-binning factor")
                    ofhdu.header["YBIN"] = (2, "Y-binning factor")
                    exptime = ihead["EXPTIME"]
                    # correct MJD to mid-exposure value, assuming value
                    # from header is at start
                    mjd = ihead["MJD"] + exptime / 2 / 86400

                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (
                        mjd,
                        "MJD at centre of exposure; fits2hcm",
                    )
                    ophdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ophdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                elif origin == "INTWFC":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "CCD number; fits2hcm")
                    exptime = ihead["EXPTIME"]
                    mjd = ihead["MJD-OBS"] + exptime / 2 / 86400
                    time = Time(mjd + exptime / 2 / 86400, format="mjd")
                    ophdu.header["TIMSTAMP"] = (time.isot, "Time stamp; fits2hcm")

                    # Copy data into first HDU
                    if hdul[0].data is not None:
                        ofhdu = fits.ImageHDU(hdul[0].data)
                    elif hdul[1].data is not None:
                        ofhdu = fits.ImageHDU(hdul[1].data)
                    else:
                        raise HipercamError(
                            'Failed to find any data in first two HDUs'
                        )

                    try:
                        wsec1 = ihead['WINSEC1']
                        wsec2 = ihead['WINSEC2']
                        wsec3 = ihead['WINSEC3']
                        wsec4 = ihead['WINSEC4']
                        nenabled = 0
                        for wsec in (wsec1,wsec2,wsec3,wsec4):
                            if wsec.find('enabled') > -1:
                                nenabled += 1
                                win = wsec

                        if nenabled != 1:
                            raise HipercamError(
                                f'File = {fname}: incorrect number window segments enabled ({nenabled})'
                            )
                        xrng,yrng,rest = wsec.split(',')
                        xl,xr = xrng[1:].split(':')
                        yb,yt = yrng.split(':')
                    except:
                        xl,yb = 1,1

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")
                    ofhdu.header["NXTOT"] = (2154, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (4200, "Total unbinned X dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (int(xl), "X-ordinate of lower-left pixel")
                    ofhdu.header["LLY"] = (int(yb), "Y-ordinate of lower-left pixel")
                    ofhdu.header["XBIN"] = (ihead["CCDXBIN"], "X-binning factor")
                    ofhdu.header["YBIN"] = (ihead["CCDYBIN"], "Y-binning factor")
                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (
                        mjd,
                        "MJD at centre of exposure; fits2hcm",
                    )
                    ofhdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ofhdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                elif origin == "PT5M":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "CCD number; fits2hcm")
                    ophdu.header["NFRAME"] = (counter, "NFRAME number; fits2hcm")
                    date_obs = ihead["DATE-OBS"]
                    t = Time(date_obs, format="isot", scale="utc")
                    exptime = ihead["EXPTIME"]
                    mjd = t.mjd + exptime / 2 / 86400
                    time = Time(mjd, format="mjd")
                    ophdu.header["TIMSTAMP"] = (time.isot, "Time stamp; fits2hcm")

                    # Copy data into first HDU
                    ofhdu = fits.ImageHDU(hdul[0].data)

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")
                    ofhdu.header["NXTOT"] = (2184, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (1472, "Total unbinned Y dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (1, "X-ordinate of lower-left pixel")
                    ofhdu.header["LLY"] = (1, "Y-ordinate of lower-left pixel")
                    ofhdu.header["XBIN"] = (ihead["XBINNING"], "X-binning factor")
                    ofhdu.header["YBIN"] = (ihead["YBINNING"], "Y-binning factor")
                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (
                        mjd,
                        "MJD at centre of exposure; fits2hcm",
                    )
                    ofhdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ofhdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                elif origin == "HICKS":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "CCD number; fits2hcm")
                    ophdu.header["NFRAME"] = (counter, "NFRAME number; fits2hcm")
                    date_obs = ihead["DATE-OBS"]
                    t = Time(date_obs, format="isot", scale="utc")
                    exptime = ihead["EXPTIME"]
                    mjd = t.mjd + exptime / 2 / 86400
                    time = Time(mjd, format="mjd")
                    ophdu.header["TIMSTAMP"] = (time.isot, "Time stamp; fits2hcm")

                    # Copy data into first HDU
                    ofhdu = fits.ImageHDU(hdul[0].data)

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")
                    ofhdu.header["NXTOT"] = (2184, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (1472, "Total unbinned Y dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (1, "X-ordinate of lower-left pixel")
                    ofhdu.header["LLY"] = (1, "Y-ordinate of lower-left pixel")
                    ofhdu.header["XBIN"] = (ihead["XBINNING"], "X-binning factor")
                    ofhdu.header["YBIN"] = (ihead["YBINNING"], "Y-binning factor")
                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (
                        mjd,
                        "MJD at centre of exposure; fits2hcm",
                    )
                    ofhdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ofhdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                elif origin == "ROSA":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "CCD number; fits2hcm")
                    ophdu.header["NFRAME"] = (counter, "NFRAME number; fits2hcm")
                    date_obs = ihead["DATE-OBS"]
                    t = Time(date_obs, format="isot", scale="utc")
                    exptime = ihead["EXPTIME"]
                    mjd = t.mjd + exptime / 2 / 86400
                    time = Time(mjd, format="mjd")
                    ophdu.header["TIMSTAMP"] = (time.isot, "Time stamp; fits2hcm")

                    # Copy data into first HDU
                    ofhdu = fits.ImageHDU(hdul[0].data)

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")
                    ofhdu.header["NXTOT"] = (2748, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (2198, "Total unbinned Y dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (1, "X-ordinate of lower-left pixel")
                    ofhdu.header["LLY"] = (1, "Y-ordinate of lower-left pixel")
                    ofhdu.header["XBIN"] = (ihead["XBINNING"], "X-binning factor")
                    ofhdu.header["YBIN"] = (ihead["YBINNING"], "Y-binning factor")
                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (
                        mjd + exptime / 2 / 86400,
                        "MJD at centre of exposure; fits2hcm",
                    )
                    ofhdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ofhdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                elif origin == "LTIO":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "Number of CCDs")
                    ophdu.header["TIMSTAMP"] = ihead["DATE-OBS"]

                    # Copy data into first HDU
                    ofhdu = fits.ImageHDU(hdul[0].data)

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")

                    # some temporaries
                    xbin = ihead["CCDXBIN"]
                    ybin = ihead["CCDYBIN"]
                    nxtot = xbin * ihead["CCDXIMSI"]
                    nytot = ybin * ihead["CCDYIMSI"]

                    ofhdu.header["NXTOT"] = (nxtot, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (nytot, "Total unbinned Y dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (
                        ihead["CCDWXOFF"] + 1,
                        "X-ordinate of lower-left pixel",
                    )
                    ofhdu.header["LLY"] = (
                        ihead["CCDWYOFF"] + 1,
                        "Y-ordinate of lower-left pixel",
                    )
                    ofhdu.header["XBIN"] = (xbin, "X-binning factor")
                    ofhdu.header["YBIN"] = (ybin, "Y-binning factor")
                    exptime = ihead["EXPTIME"]
                    mjd = ihead["MJD"] + exptime / 2 / 86400
                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ofhdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ofhdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                elif origin == "LCOGT":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[1].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "CCD number; fits2hcm")
                    exptime = ihead["EXPTIME"]
                    mjd = ihead["MJD-OBS"] + exptime / 2 / 86400
                    time = Time(mjd + exptime / 2 / 86400, format="mjd")
                    ophdu.header["TIMSTAMP"] = (time.isot, "Time stamp; fits2hcm")

                    # Copy data into first HDU
                    ofhdu = fits.ImageHDU(hdul[1].data)

                    NXTOT = ihead['NAXIS1']
                    NYTOT = ihead['NAXIS2']

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")
                    ofhdu.header["NXTOT"] = (NXTOT, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (NYTOT, "Total unbinned Y dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (1, "X-ordinate of lower-left pixel")
                    ofhdu.header["LLY"] = (1, "Y-ordinate of lower-left pixel")
                    ofhdu.header["XBIN"] = (1, "X-binning factor")
                    ofhdu.header["YBIN"] = (1, "Y-binning factor")
                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (
                        mjd,
                        "MJD at centre of exposure; fits2hcm",
                    )
                    ofhdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ofhdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                elif origin == "WHTCAM":

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header["NUMCCD"] = (1, "CCD number; fits2hcm")
                    exptime = ihead["EXPTIME"]
                    mjd = ihead["MJD-OBS"] + exptime / 2 / 86400
                    time = Time(mjd + exptime / 2 / 86400, format="mjd")
                    ophdu.header["TIMSTAMP"] = (time.isot, "Time stamp; fits2hcm")

                    # Copy data into first HDU
                    if hdul[0].data is not None:
                        ofhdu = fits.ImageHDU(hdul[0].data)
                    elif hdul[1].data is not None:
                        ofhdu = fits.ImageHDU(hdul[1].data)
                    else:
                        raise HipercamError(
                            'Failed to find any data in first two HDUs'
                        )

                    # no window info

                    # Get header into right format
                    ofhdu.header["CCD"] = ("1", "CCD label")
                    ofhdu.header["NXTOT"] = (3200, "Total unbinned X dimension")
                    ofhdu.header["NYTOT"] = (3200, "Total unbinned Y dimension")
                    ofhdu.header["NUMWIN"] = (1, "Total number of windows")
                    ofhdu.header["WINDOW"] = ("1", "Window label")
                    ofhdu.header["LLX"] = (1, "X-ordinate of lower-left pixel")
                    ofhdu.header["LLY"] = (1, "Y-ordinate of lower-left pixel")
                    ofhdu.header["XBIN"] = (ihead["XBINNING"], "X-binning factor")
                    ofhdu.header["YBIN"] = (ihead["YBINNING"], "Y-binning factor")
                    ofhdu.header["MJDUTC"] = (mjd, "MJD at centre of exposure")
                    ophdu.header["MJDUTC"] = (
                        mjd,
                        "MJD at centre of exposure; fits2hcm",
                    )
                    ofhdu.header["MJDINT"] = (
                        int(mjd),
                        "Integer part of MJD at centre of exposure",
                    )
                    ofhdu.header["MJDFRAC"] = (
                        mjd - int(mjd),
                        "Fractional part of MJD at centre of exposure",
                    )
                    ofhdu.header["EXPTIME"] = (exptime, "Exposure time, seconds")
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                else:
                    raise HipercamError(
                        (
                            "Origin = {:s} unrecognised. Codes " "recognised: {:s}"
                        ).format(origin, ", ".join(FORMATS))
                    )

                print(fname, "-->", oname)
