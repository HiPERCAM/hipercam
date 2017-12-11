import unittest
import copy

import numpy as np
from astropy.io import fits

from hipercam import Winhead, Window, CCD, MCCD, HipercamError

class TestMCCD(unittest.TestCase):
    """Provides simple tests of MCCD methods and attributes.

    """

    def setUp(self):

        # build an MCCD
        self.level1 = 10.
        win1 = Winhead(31,41,5,5,1,2)
        data1 = np.empty((win1.ny,win1.nx))
        data1[:] = self.level1
        wind1 = Window(win1,data1)

        self.level3 = 20.
        win2 = Winhead(250,561,5,5,1,2)
        data2 = np.empty((win1.ny,win1.nx))
        data2[:] = self.level3
        wind2 = Window(win2,data2)

        winds = ((1,wind1),(3,wind2))

        head = fits.Header()
        head['OBJECT'] = ('Fake 1','Name')
        head['FILTER'] = ('r','Filter')
        ccd1 = CCD(winds, 2048, 1024, head)

        ccd2 = ccd1.copy()
        ccd2.head['FILTER'] = ('g','Filter')

        ccds = ((3,ccd1),(5,ccd2))
        head = fits.Header()
        head['TEST'] = 'This is HiPERCAM'

        self.mccd = MCCD(ccds,head)

    def test_mccd_copy(self):

        mccd = self.mccd.copy()
        mccd[3][1].data[0,0] = 0.
        self.assertEqual(mccd[3][1].data[0,0]+self.level1,self.mccd[3][1].data[0,0],
                         'copy of MCCD failed')

        mccd = copy.copy(self.mccd)
        mccd[3][1].data[0,0] = 0.
        self.assertEqual(mccd[3][1].data[0,0]+self.level1,self.mccd[3][1].data[0,0],
                         'copy.copy of MCCD failed')

        mccd = copy.deepcopy(self.mccd)
        mccd[3][1].data[0,0] = 0.
        self.assertEqual(mccd[3][1].data[0,0]+self.level1,self.mccd[3][1].data[0,0],
                         'copy.deepcopy of MCCD failed')


if __name__ == '__main__':
    unittest.main()
