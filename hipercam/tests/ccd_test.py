import unittest
import copy

import numpy as np
from astropy.io import fits


from hipercam import Window, Windat, CCD, HipercamError

class TestCCD(unittest.TestCase):
    """Provides simple tests of CCD methods and attributes.

    """

    def setUp(self):

        # build a CCD
        win1 = Window(31,41,5,5,1,2)
        data1 = np.empty((win1.ny,win1.nx))
        data1[:] = 10.
        wind1 = Windat(win1,data1)

        win2 = Window(250,561,5,5,1,2)
        data2 = np.empty((win1.ny,win1.nx))
        data2[:] = 20.
        wind2 = Windat(win2,data2)

        winds = ((1,wind1),(3,wind2))

        head = fits.Header()
        head['OBJECT'] = ('Fake 1','Name')
        self.ccd = CCD(winds, 2048, 1024, head)

    def test_ccd_copy(self):
        ccd = self.ccd.copy()
        ccd[1].data += 1.
        self.assertNotEqual(ccd[1].data[0,0],self.ccd[1].data[0,0],
                            'copy of CCD failed')

        ccd = copy.copy(self.ccd)
        ccd[1].data += 1.
        self.assertNotEqual(ccd[1].data[0,0],self.ccd[1].data[0,0],
                            'copy.copy of CCD failed')

        ccd = copy.deepcopy(self.ccd)
        ccd[1].data += 1.
        self.assertNotEqual(ccd[1].data[0,0],self.ccd[1].data[0,0],
                            'copy.deepcopy of CCD failed')

    def test_ccd_iadd(self):
        ccd = self.ccd.copy()
        ccd += 5.
        self.assertTrue(ccd[1].data[0,0] == 15. and ccd[3].data[0,0] == 25.,
                        'in place addition of constant to CCD failed')

    def test_ccd_isub(self):
        ccd = self.ccd.copy()
        ccd -= 5.
        self.assertTrue(ccd[1].data[0,0] == 5. and ccd[3].data[0,0] == 15.,
                        'in place subtraction of constant from CCD failed')

if __name__ == '__main__':
    unittest.main()
