import unittest

from hipercam import Window, HipercamError

class TestWindow(unittest.TestCase):

    def setUp(self):
        # win1 and win2 should clash (they overlap)
        # win1 is inside win3
        self.win1 = Window(3,4,100,200,1,2)
        self.win2 = Window(50,50,100,200,1,2)
        self.win3 = Window(1,2,200,300,1,2)


    def test_window_llx_value(self):
        self.assertEqual(self.win1.llx, 3,
                         'incorrect lower-left X pixel')

    def test_window_lly_value(self):
        self.assertEqual(self.win1.lly, 4,
                         'incorrect lower-left Y pixel')

    def test_window_nx_value(self):
        self.assertEqual(self.win1.nx, 100,
                         'incorrect number of pixels in X')

    def test_window_ny_value(self):
        self.assertEqual(self.win1.ny, 200,
                         'incorrect number of pixels in Y')

    def test_window_xbin_value(self):
        self.assertEqual(self.win1.xbin, 1,
                         'incorrect X binning factor')

    def test_window_ybin_value(self):
        self.assertEqual(self.win1.ybin, 2,
                         'incorrect Y binning factor')

    def test_window_clash(self):
        self.assertRaises(HipercamError, self.win1.clash, self.win2)

    def test_window_urx_value(self):
        self.assertEqual(self.win1.urx, 102,
                         'incorrect upper-right X pixel')

    def test_window_ury_value(self):
        self.assertEqual(self.win1.ury, 403,
                         'incorrect upper-right Y pixel')

    def test_window_extent_value(self):
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



if __name__ == '__main__':
    unittest.main()
