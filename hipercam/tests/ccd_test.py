import unittest
import copy
from astropy.io import fits

import numpy as np
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

        head = fits.Header()
        head['OBJECT'] = ('Fake 1','Name')
        self.ccd = CCD(((1,wind1),(3,wind2)),2048,1024,head)
        print(self.ccd)

    def test_ccd_iadd(self):
        ccd = copy.deepcopy(self.ccd)
        ccd += 5.
        self.assertTrue(ccd[1].data[0,0] == 15. and ccd[3].data[0,0] == 25.,
                        'in place addition of constant to CCD failed')

if __name__ == '__main__':
    unittest.main()
