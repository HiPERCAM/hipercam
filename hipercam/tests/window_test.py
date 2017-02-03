import unittest
import copy

import numpy as np
from hipercam import Window, HipercamError

class TestWindow(unittest.TestCase):
    """
    Provides simple tests of every Window method and attribute. If any of
    these fail, something is badly wrong.
    """

    def setUp(self):
        # win1 and win2 should clash (they overlap)
        # win1 is inside win3
        # win1c is a copy of win1
        self.win1 = Window(3,4,100,200,1,2)
        self.win1c = self.win1.copy()
        self.win2 = Window(50,50,100,200,1,2)
        self.win3 = Window(1,2,200,300,1,2)

        # win4 is tiny for a quick test of xy; x and y are the
        # expected returns.
        self.win4 = Window(11,21,2,3,2,3)
        self.x = np.array([[11.5,13.5],[11.5,13.5],[11.5,13.5]])
        self.y = np.array([[22.,22.],[25.,25.],[28.,28.]])

    def test_window_llx(self):
        self.assertEqual(self.win1.llx, 3,
                         'incorrect lower-left X pixel')

    def test_window_lly(self):
        self.assertEqual(self.win1.lly, 4,
                         'incorrect lower-left Y pixel')

    def test_window_nx(self):
        self.assertEqual(self.win1.nx, 100,
                         'incorrect number of pixels in X')

    def test_window_ny(self):
        self.assertEqual(self.win1.ny, 200,
                         'incorrect number of pixels in Y')

    def test_window_xbin(self):
        self.assertEqual(self.win1.xbin, 1,
                         'incorrect X binning factor')

    def test_window_ybin(self):
        self.assertEqual(self.win1.ybin, 2,
                         'incorrect Y binning factor')

    def test_window_clash(self):
        self.assertRaises(HipercamError, self.win1.clash, self.win2)

    def test_window_urx(self):
        self.assertEqual(self.win1.urx, 102,
                         'incorrect upper-right X pixel')

    def test_window_ury(self):
        self.assertEqual(self.win1.ury, 403,
                         'incorrect upper-right Y pixel')

    def test_window_extent(self):
        self.assertEqual(self.win1.extent(), (2.5, 102.5, 3.5, 403.5),
                         'incorrect window extent')

    def test_window_inside(self):
        self.assertTrue(self.win1.inside(self.win3),
                         'win1 should be inside win3')
        self.assertFalse(self.win1.inside(self.win2),
                         'win1 should not be inside win2')

    def test_window_outside(self):
        self.assertTrue(self.win3.outside(self.win1),
                         'win3 should be outside win1')
        self.assertFalse(self.win2.outside(self.win1),
                         'win2 should not be outside win1')

    def test_window_xy(self):
        x,y = self.win4.xy()

        self.assertTrue((x == self.x).all(),
                        'incorrect X array returned')
        self.assertTrue((y == self.y).all(),
                        'incorrect Y array returned')

    def test_window_matches(self):
        self.assertRaises(HipercamError, self.win1.matches, self.win2)

    def test_window_copy_eq(self):
        self.assertTrue(self.win1 == self.win1c)

    def test_window_ne(self):
        self.assertFalse(self.win1 != self.win1c)

if __name__ == '__main__':
    unittest.main()
