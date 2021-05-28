import sys
import os
import tempfile

import numpy as np
from astropy.time import Time
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

import hipercam as hcam
from hipercam import cline, utils, spooler, defect, fringe
from hipercam.cline import Cline

__all__ = [
    "hpackage",
]

##############################################################
#
# hpackage -- bundles up standard data from a run on an object
#
##############################################################


def hpackage(args=None):
    """``hpackage run``

    ``hpackage`` looks for standard data products run#.hcm, run#.ape,
    run#.red, run#.log, (where # is an integer identifier), spits out
    joined up images of the CCDs in run#.hcm, converts the log to fits,
    and writes out a "region" file representing the apertures readable by
    ds9. It then tars the whole lot up. For simplicity, it deliberately runs
    with a minimal set or arguments: the run name. It won't overwrite any
    file.

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("run", Cline.LOCAL, Cline.PROMPT)
        run = cl.get_value("run", "run name", cline.Fname("run005", hcam.HCAM))
        mccd = hcam.MCCD.read(run)

    # check for existence of various files
    if not os.path.exists(run + HCAM.APE):
        raise hcam.HipercamError(f'Could not find file = {run + hcam.APE}')

    # check for existence of various files
    if not os.path.exists(run + HCAM.LOG):
        raise hcam.HipercamError(f'Could not find file = {run + hcam.LOG}')

    with tempfile.TemporaryDirectory() as tmpdir:
        # create temporary directory where everything
        # gets written
        print('Will write files to {tmpdir}')

        # First make individual CCDs from .hcm file

        for nc, cnam in enumerate(ccds):
            ccd = mccd[cnam]

            # We need to generate a single data array for all data. First
            # check that it is even possible

            for n, wind in enumerate(ccd.values()):
                if n == 0:
                    llx = llxmin = wind.llx
                    lly = llymin = wind.lly
                    urxmax = wind.urx
                    urymax = wind.ury
                    xbin = wind.xbin
                    ybin = wind.ybin
                else:
                    # Track overall dimensions
                    llxmin = min(llxmin, wind.llx)
                    llymin = min(llymin, wind.lly)
                    urxmax = max(urxmax, wind.urx)
                    urymax = max(urymax, wind.ury)
                    if xbin != wind.xbin or ybin != wind.ybin:
                        raise hcam.HipercamError('Found windows with clashing binning factors')
                    if (wind.llx - llx) % xbin != 0 or (wind.lly - lly) % ybin != 0:
                        raise hcam.HipercamError('Found windows which are out of sync with each other')

            # create huge array of nothing
            ny = (urymax-llymin+1) // ybin
            nx = (urxmax-llxmin+1) // xbin
            data = np.zeros((ny,nx))

            # fill it with data
            for n, wind in enumerate(ccd.values()):
                xstart = (wind.llx - llxmin) // xbin
                ystart = (wind.lly - llymin) // ybin
                data[ystart:ystart+wind.ny,xstart:xstart+wind.nx] = wind.data
            data = data.astype(np.float32)

            # Header
            phead = mccd.head.copy()

            # Add some extra stuff
            phead["CCDLABEL"] = (cnam, "CCD label")
            phead["NFRAME"] = (nf, "Frame number")
            phead["LLX"] = (llxmin, "X of lower-left unbinned pixel (starts at 1)")
            phead["LLY"] = (llymin, "Y of lower-left unbinned pixel (starts at 1)")
            phead["NXTOT"] = (ccd.nxtot, "Total unbinned X dimension")
            phead["NYTOT"] = (ccd.nytot, "Total unbinned Y dimension")
            phead.add_comment('Written by HiPERCAM script "hpackage"')

            # Make header
            header=fits.Header(phead.cards)

            # attempt to generate a WCS
            instrume = header.get('INSTRUME','UNKNOWN')
            if instrume == 'HIPERCAM' and 'RADEG' in header and \
               'DECDEG' in header and 'INSTRPA' in header:

                ZEROPOINT = 69.3 # rotator zeropoint
                SCALE = 0.081 # pixel scale, "/unbinned pixel

                # ra, dec in degrees at rotator centre
                ra = header['RADEG']
                dec = header['DECDEG']

                # position angle, degrees
                pa = header['INSTRPA'] - ZEROPOINT

                # position within binned array of rotator centre
                x0, y0 = (1020.-llxmin+1)/xbin, (524.-llymin+1)/ybin

                # Build the WCS
                w = wcs.WCS(naxis=2)
                w.wcs.crpix = [x0, y0]
                w.wcs.crval = [ra,dec]
                cpa = np.cos(np.radians(pa))
                spa = np.sin(np.radians(pa))
                cd = np.array([[xbin*cpa,-ybin*spa],[xbin*spa,ybin*cpa]])
                cd = np.array([[xbin*cpa,ybin*spa],[-xbin*spa,ybin*cpa]])
                cd *= SCALE/3600
                w.wcs.cd = cd
                w.wcs.ctype = ['RA---TAN','DEC--TAN']
                w.wcs.cunit = ["deg", "deg"]
                header.update(w.to_header())
                print(f'   CCD {cnam}: added WCS')
            else:
                print(f'   CCD {cnam}: missing positional data; no WCS added')

            # make the first & only HDU
            hdul = fits.HDUList()
            hdul.append(fits.PrimaryHDU(header=header))
            compressed_hdu = fits.CompImageHDU(
                data=data, compression_type='RICE_1'
            )
            hdul.append(compressed_hdu)

            root = os.path.basename(os.path.splitext(run])[0])
            oname = os.path.join(tmpdir, f'{root}_ccd{cnam}.fits')

            hdul.writeto(oname, overwrite=overwrite)
            print(f'   CCD {cnam}: written to {oname}')

