import sys
import os

import numpy as np
from astropy.io import fits

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

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

      flist : string
         name of list of input FITS-format files.
         'flist' should end '.lis'. The output file names will
         have the same rootname (without any leading directories)
         but ends '.hcm' 

      origin : string
         origin of the data. Currently recognised:

           LTRISE :
             Liverpool telescope RISE camera. I don't know how general my
             code is in this case. I assume XBIN=YBIN=2. If you find it does
             not work, let me know.

      overwrite : bool
         overwrite files on output
    """

    command, args = utils.script_args(args)

    # get input section
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('flist', Cline.LOCAL, Cline.PROMPT)
        cl.register('origin', Cline.LOCAL, Cline.PROMPT)
        cl.register('overwrite', Cline.LOCAL, Cline.PROMPT)

        # get inputs
        flist = cl.get_value(
            'flist', 'list of input FITS files to convert to hcm',
            cline.Fname('files', hcam.LIST)
        )

        origin = cl.get_value(
            'origin', 'origin of data', 'LTRISE',
            lvals=['LTRISE',]
        )

        overwrite = cl.get_value(
            'overwrite', 'overwrite data on output', True
        )

    with open(flist) as fin:
        for line in fin:
            fname = line.strip()
            with fits.open(fname) as hdul:
                bname = os.path.basename(fname)
                if bname.find('.') > -1:
                    oname = bname[:bname.find('.')] + hcam.HCAM
                else:
                    oname = bname + hcam.HCAM

                if origin == 'LTRISE':

                    # Copy main header into primary data-less HDU
                    ihead = hdul[0].header
                    ophdu = fits.PrimaryHDU(header=ihead)
                    ophdu.header['NUMCCD'] = (1, 'Number of CCDs')
                    ophdu.header['TIMSTAMP'] = ihead['DATE-OBS']

                    # Copy data into first HDU
                    ofhdu = fits.ImageHDU(hdul[0].data)

                    # Get header into right format
                    ofhdu.header['CCD'] = ('1','CCD label')
                    ofhdu.header['NXTOT'] = (1026, 'Total unbinned X dimension')
                    ofhdu.header['NYTOT'] = (1026, 'Total unbinned X dimension')
                    ofhdu.header['NUMWIN'] = (1, 'Total number of windows')
                    ofhdu.header['WINDOW'] = ('1', 'Window label')
                    ofhdu.header['LLX'] = (
                        ihead['CCDWXOFF']+1,
                        'X-ordinate of lower-left pixel'
                    )
                    ofhdu.header['LLY'] = (
                        ihead['CCDWYOFF']+1,
                        'Y-ordinate of lower-left pixel'
                    )
                    ofhdu.header['XBIN'] = (2, 'X-binning factor')
                    ofhdu.header['YBIN'] = (2, 'Y-binning factor')
                    mjd = ihead['MJD']
                    exptime = ihead['EXPTIME']
                    ofhdu.header['MJDUTC'] = (
                        mjd+exptime/2/86400,'MJD at centre of exposure'
                    )
                    ophdu.header['MJDUTC'] = (
                        mjd+exptime/2/86400,'MJD at centre of exposure'
                    )
                    ofhdu.header['EXPTIME'] = (
                        exptime, 'Exposure time, seconds'
                    )
                    ohdul = fits.HDUList([ophdu, ofhdu])
                    ohdul.writeto(oname, overwrite=overwrite)

                else:
                    raise HipercamError(
                        'Origin = {:s} unrecognised'.format(origin)
                    )

                print(fname,'-->',oname)
