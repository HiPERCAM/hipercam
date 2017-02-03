import unittest
import copy

import numpy as np
from hipercam import Target, HipercamError

class TestTarget(unittest.TestCase):
    """
    Provides unit tests of the Target class
    """

    def setUp(self):
        self.xcen = 51.
        self.ycen = 65.
        self.height = 113.
        self.fwhm1 = 4.5
        self.fwhm2 = 5.5
        self.angle = 40.
        self.beta = 4.

        self.targ = Target(self.xcen,self.ycen,self.height,self.fwhm1,
                           self.fwhm2,self.angle,self.beta)

    def test_target_xcen(self):
        self.assertEqual(self.targ.xcen, self.xcen,
                         'xcen incorrectly stored / retrieved')

    def test_target_ycen(self):
        self.assertEqual(self.targ.ycen, self.ycen,
                         'ycen incorrectly stored / retrieved')

    def test_target_fwhm1(self):
        self.assertEqual(self.targ.fwhm1, self.fwhm1,
                         'fwhm1 incorrectly stored')

    def test_target_fwhm2(self):
        self.assertEqual(self.targ.fwhm2, self.fwhm2,
                         'fwhm2 incorrectly stored')

    def test_target_angle(self):
        self.assertEqual(self.targ.angle, self.angle,
                         'angle incorrectly stored')

    def test_target_beta(self):
        self.assertEqual(self.targ.beta, self.beta,
                         'beta incorrectly stored')

    def test_target_copy(self):
        targ = self.targ.copy()
        self.assertEqual(self.targ.height, targ.height,
                         'target height did not transfer')

    def test_target_call(self):
        # calculates project at half FWHM1 along axis 1.
        # should therefore return height/2.
        x = self.xcen+(self.fwhm1/2)*np.cos(np.radians(self.angle))
        y = self.ycen+(self.fwhm1/2)*np.sin(np.radians(self.angle))
        val = self.targ(x,y)

        self.assertTrue(abs(val-self.targ.height/2.) < 1.e-12*self.targ.height,
                         'profile did not a factor of 2 as expected')




if __name__ == '__main__':
    unittest.main()
