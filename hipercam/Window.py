"""
Defines a class to represent one window of a CCD
"""

# Standard pre-amble from astropy
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

class Window(object):
    """
    Class representing one window of a CCD. This represents an arbitrary
    rectangular region of pixels out of some fixed total. The lower-left
    pixel of the CCD is assumed to be have coordinates (x,y) = (1,1). Window
    dimensions are in binned pixels.
    """

    def __init__(self, llx, lly, nx, ny, xbin, ybin):
        """
        Constructor. Arguments::

        llx : (int)
            X position of lower-left pixel of window

        lly : (int)
            Y position of lower-left pixel of window

        nx : (int)
            X dimension of window, binned pixels (> 0)

        ny : (int)
            Y dimension of window, binned pixels (> 0)

        xbin : (int)
            Binning factor in X (> 0)

        ybin : (int)
            Binning factor in Y (> 0)
        """

        self.llx   = llx
        self.lly   = lly
        self.nx    = nx
        self.ny    = ny
        self.xbin  = xbin
        self.ybin  = ybin

    @property
    def urx(self):
        """
        Returns (unbinned) X pixel at upper-right of Window
        """
        return self.llx-1+self.nx*self.xbin

    @property
    def ury(self):
        """
        Returns (unbinned) Y pixel at upper-right of Window
        """
        return self.lly-1+self.ny*self.ybin

    def wfhead(self, nccd, nwin, fhead, nxccd, nyccd, NX, NY):
        """
        Writes window data into a FITS header. Assumes HiPERCAM type
        data in the sense that the amplifier is computed according to
        the window number. Non-HiPERCAM data will still work but the
        amplifier assignment will be wrong.

        Arguments::

          nccd : (int)
             CCD number, starting from 0

          nwin : (int)
             Window number, starting from 0

          fhead : (astropy.io.fits.Header)
             FITS header object. Will be modified on exit.

          nxccd : (int)
             X dimension of CCD

          nyccd : (int)
             Y dimension of CCD

          NX : (int)
             X dimension for iraf mosaic format

          NY : (int)
             Y dimension for iraf mosaic format

        Returns::

          fhead : (astropy.io.fits.Header)
             The modified FITS header object.
        """

        # The next line will likely need modifying
        namp = nwin % 4 + 1

        fhead['NCCD']    = (nccd+1,'CCD number')
        fhead['NWIN']    = (nwin+1,'Window number')
        fhead['LLX']     = (self.llx,  'X-coordinate of lower-left pixel')
        fhead['LLY']     = (self.lly,  'Y-coordinate of lower-left pixel')
        fhead['XBIN']    = (self.xbin, 'Binning factor in X')
        fhead['YBIN']    = (self.ybin, 'Binning factor in Y')

        fhead['CCDNAME'] = (str(nccd+1),'CCD identifier')
        fhead['AMPNAME'] = (str(namp),'Amplifier identifier')
        fhead['CCDSIZE'] = ('[1:' + str(nxccd) + ',1:' + str(nyccd) + ']', 'CCD size')

        # Flags to define the amplifier quadrant
        left   = nwin % 4 == 0 or nwin % 4 == 2
        bottom = nwin % 4 == 0 or nwin % 4 == 1

        # location within the NX by NY array of CCDs
        nx = nccd % NX
        ny = nccd // NY

        # Now we start on the iraf mosaic stuff, see
        # http://iraf.noao.edu/projects/ccdmosaic/imagedef/imagedef.html
        # for more detail.
        # Basically, 4 coordinate systems: CCD, Amplifier, Data, Detector
        #
        # CCD  -- basic unbinned 1,1 in lower-left corner
        # Amp  -- same, but measured with respect to amplifier
        # Data -- 

        # head
        ccdsec = ampsec = datasec = detsec = '['

        # x range
        ccdsec += str(self.llx) + ':' + str(self.llx+self.nx*self.xbin-1)
        if left:
            ampsec += str(self.llx) + ':' + str(self.llx+self.nx*self.xbin-1)
        else:
            ampsec += str(nxccd-self.llx+1) + ':' + str(nxccd-self.llx-self.nx*self.xbin+2)
        datasec += '1:' + str(self.nx)
        detsec  += str(self.llx+nx*nxccd) + ':' + str(self.llx+nx*nxccd+self.nx*self.xbin-1)

        # middle
        ccdsec += ','
        ampsec += ','
        datasec += ','
        detsec += ','

        # y range
        ccdsec += str(self.lly) + ':' + str(self.lly+self.nx*self.ybin-1)
        if bottom:
            ampsec += str(self.lly) + ':' + str(self.lly+self.ny*self.ybin-1)
        else:
            ampsec += str(nyccd-self.lly+1) + ':' + str(nyccd-self.lly-self.ny*self.ybin+2)
        datasec += '1:' + str(self.ny)
        detsec  += str(self.lly+ny*nyccd) + ':' + str(self.lly+ny*nyccd+self.ny*self.ybin-1)

        # tail
        ccdsec += ']'
        ampsec += ']'
        datasec += ']'
        detsec += ']'


        fhead['CCDSEC']   = ccdsec
        fhead['AMPSEC']   = ccdsec
        fhead['DATASEC']  = datasec
        fhead['DETSEC']   = detsec
        fhead['CCDSUM']   = str(self.xbin) + ' ' + str(self.ybin)

#        fhead['CTYPE1']   = 'Pixel'
#        fhead['CTYPE2']   = 'Pixel'
#        fhead['CRPIX1']   = 1
#        fhead['CRPIX2']   = 1
#        fhead['CRVAL1']   = xref
#        fhead['CRVAL2']   = yref
#        fhead['CDELT1']   = self.xbin
#        fhead['CDELT2']   = self.ybin

        if left:
            fhead['ATM1_1']   = 1.
            fhead['ATV1']     = 0.
        else:
            fhead['ATM1_1']   = -1.
            fhead['ATV1']     = nxccd+1

        if bottom:
            fhead['ATM2_2']   = 1.
            fhead['ATV2']     = 0.
        else:
            fhead['ATM2_2']   = -1.
            fhead['ATV2']     = nyccd+1

        fhead['LTM1_1']   = 1./self.xbin
        fhead['LTM2_2']   = 1./self.ybin
        fhead['LTV1']     = (1-self.llx)/self.xbin
        fhead['LTV2']     = (1-self.lly)/self.ybin
        fhead['DTM1_1']   = 1.
        fhead['DTM2_2']   = 1.
        fhead['DTV1']     = nx*nxccd
        fhead['DTV2']     = ny*nyccd

        return fhead

def overlap(win1, win2):
    """
    Returns True if windows win1 and win2 overlap, False if they don't. Overlap
    means that they have 1 or more pixels in common.
    """
    return win1.llx <=  win2.urx and win1.urx >= win2.llx and \
        win1.lly <=  win2.ury and win1.ury >= win2.lly

if __name__ == '__main__':
    win1 = Window(1, 2, 10, 10, 3, 4)
    win2 = Window(10, 20, 10, 10, 3, 5)
    print('overlap(win1,win2) =',overlap(win1,win2))

