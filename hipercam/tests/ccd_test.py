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
        self.level1 = 10.
        win1 = Window(31,41,5,5,1,2)
        data1 = np.empty((win1.ny,win1.nx))
        data1[:] = self.level1
        wind1 = Windat(win1,data1)

        self.level3 = 20.
        win2 = Window(250,561,5,5,1,2)
        data2 = np.empty((win1.ny,win1.nx))
        data2[:] = self.level3
        wind2 = Windat(win2,data2)

        winds = ((1,wind1),(3,wind2))

        head = fits.Header()
        head['OBJECT'] = ('Fake 1','Name')
        self.ccd = CCD(winds, 2048, 1024, head)

    def test_ccd_copy(self):
        ccd = self.ccd.copy()
        ccd[1].data[0,0] = 0.
        self.assertEqual(ccd[1].data[0,0]+self.level1,self.ccd[1].data[0,0],
                         'copy of CCD failed')

        ccd = copy.copy(self.ccd)
        ccd[1].data[0,0] = 0.
        self.assertEqual(ccd[1].data[0,0]+self.level1,self.ccd[1].data[0,0],
                         'copy.copy of CCD failed')

        ccd = copy.deepcopy(self.ccd)
        ccd[1].data[0,0] = 0.
        self.assertEqual(ccd[1].data[0,0]+self.level1,self.ccd[1].data[0,0],
                         'copy.deepcopy of CCD failed')

    def test_ccd_iadd(self):
        self.ccd += 5.
        self.assertTrue(self.ccd[1].data[0,0] == self.level1 + 5.
                        and self.ccd[3].data[0,0] == self.level3 + 5.,
                        'in place addition of a constant to a CCD failed')

    def test_ccd_isub(self):
        self.ccd -= 5.
        self.assertTrue(self.ccd[1].data[0,0] == self.level1 - 5.
                        and self.ccd[3].data[0,0] == self.level3 - 5.,
                        'in place subtraction of a constant from a CCD failed')

    def test_ccd_imul(self):
        self.ccd *= 2.
        self.assertTrue(self.ccd[1].data[0,0] == 2.*self.level1
                        and self.ccd[3].data[0,0] == 2.*self.level3,
                        'in place multiplication of a CCD by a constant failed')

    def test_ccd_idiv(self):
        self.ccd /= 2.
        self.assertTrue(self.ccd[1].data[0,0] == self.level1/2.
                        and self.ccd[3].data[0,0] == self.level3/2.,
                        'in place division of a CCD by a constant failed')

    def test_ccd_add(self):
        ccd = self.ccd + 5.
        self.assertTrue(ccd[1].data[0,0] == self.level1 + 5.
                        and ccd[3].data[0,0] == self.level3 + 5.,
                        'addition of a constant to a CCD failed')

if __name__ == '__main__':
    unittest.main()
