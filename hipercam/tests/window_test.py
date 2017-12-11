import unittest
import copy

import numpy as np
from hipercam import Winhead

class TestWinhead(unittest.TestCase):
    """
    Provides simple tests of every Winhead method and attribute. If any of
    these fail, something is badly wrong.
    """

    def setUp(self):

        # set up a window with named parameters
        self.llx = 3
        self.lly = 4
        self.nx = 100
        self.ny = 200
        self.xbin = 1
        self.ybin = 2
        self.win = Winhead(self.llx,self.lly,self.nx,self.ny,self.xbin,self.ybin)

    def test_window_llx(self):
        self.assertEqual(self.win.llx, self.llx,
                         'incorrect lower-left X pixel')

    def test_window_lly(self):
        self.assertEqual(self.win.lly, self.lly,
                         'incorrect lower-left Y pixel')

    def test_window_nx(self):
        self.assertEqual(self.win.nx, self.nx,
                         'incorrect number of pixels in X')

    def test_window_ny(self):
        self.assertEqual(self.win.ny, self.ny,
                         'incorrect number of pixels in Y')

    def test_window_xbin(self):
        self.assertEqual(self.win.xbin, self.xbin,
                         'incorrect X binning factor')

    def test_window_ybin(self):
        self.assertEqual(self.win.ybin, self.ybin,
                         'incorrect Y binning factor')

    def test_window_clash(self):
        win = self.win.copy()
        win.llx += 1
        self.assertRaises(ValueError, self.win.clash, win)

    def test_window_urx(self):
        self.assertEqual(self.win.urx, self.llx + self.xbin*self.nx - 1,
                         'incorrect upper-right X pixel')

    def test_window_ury(self):
        self.assertEqual(self.win.ury, self.lly + self.ybin*self.ny - 1,
                         'incorrect upper-right Y pixel')

    def test_window_extent(self):
        x1 = self.llx-0.5
        x2 = self.win.urx+0.5
        y1 = self.lly-0.5
        y2 = self.win.ury+0.5
        self.assertEqual(self.win.extent(), (x1,x2,y1,y2),
                         'incorrect window extent')

    def test_window_inside(self):
        win = self.win.copy()
        win.llx -= self.xbin
        win.nx += 2
        win.lly -= self.ybin
        win.ny += 2
        self.assertTrue(self.win.inside(win),
                         'first inside test failed')

        win = self.win.copy()
        win.llx += self.nx*self.xbin
        self.assertFalse(self.win.inside(win),
                         'second inside test failed')

    def test_window_outside(self):
        win = self.win.copy()
        win.llx -= self.xbin
        win.nx += 2
        win.lly -= self.ybin
        win.ny += 2
        self.assertTrue(win.outside(self.win),
                         'first outside test failed')

        win = self.win.copy()
        win.llx += self.nx*self.xbin
        self.assertFalse(win.outside(self.win),
                         'second outside test failed')

    def test_window_xy(self):
        # tiny window to test xy
        win = Winhead(11,21,2,3,2,3)
        X = np.array([[11.5,13.5],[11.5,13.5],[11.5,13.5]])
        Y = np.array([[22.,22.],[25.,25.],[28.,28.]])
        x,y = win.xy()
        self.assertTrue((x == X).all(),
                        'incorrect X array returned')
        self.assertTrue((y == Y).all(),
                        'incorrect Y array returned')

    def test_window_matches(self):
        win = self.win.copy()
        win.llx += 1
        self.assertRaises(ValueError, self.win.matches, win)

    def test_window_copy_eq(self):
        win = self.win.copy()
        self.assertTrue(self.win == win)

    def test_window_ne(self):
        win = self.win.copy()
        self.assertFalse(self.win != win)

if __name__ == '__main__':
    unittest.main()
